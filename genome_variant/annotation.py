"""GFF3-driven variant annotation for parsed VCF documents.

Public API:

- :class:`AnnotationFormatError` — malformed GFF3 input (bad CDS rows,
  coordinates or phase, missing Parent, a strand other than ``+``/``-``,
  CDS fragments spanning several sequences or contradictory reading
  frames) or a VCF that already declares or uses the ``GVANN`` key.
- :func:`annotate_vcf` — append one ``GVANN`` INFO item per ALT allele to
  every record of a :class:`~genome_variant.vcf.VcfFile`.

CDS features are grouped by their GFF3 ``Parent`` attribute into
transcripts.  The standard genetic code is used.  ``A``/``C``/``G``/``T``
single-base substitutions inside a CDS are annotated with a coding
consequence, and pure ``A``/``C``/``G``/``T`` insertions and deletions
are annotated as ``FRAMESHIFT`` (length change not a multiple of three,
``HIGH``) or ``INFRAME_INSERTION``/``INFRAME_DELETION`` (``MODERATE``)
when the event lies inside one contiguous CDS fragment.  Every other ALT
(complex replacements, ambiguous or symbolic alleles, indels that cross
CDS fragment boundaries or miss a coding anchor) receives
``UNSUPPORTED``/``MODIFIER`` and records without a CDS hit receive
``NON_CODING``/``MODIFIER``.
"""

from __future__ import annotations

import io
import os
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace

from .sequence_io import SequenceRecord
from .vcf import (
    ReferenceMismatchError,
    VcfFile,
    VcfFormatError,
    VcfHeader,
    VcfRecord,
)

__all__ = [
    "AnnotationFormatError",
    "annotate_vcf",
]

GVANN_ID = "GVANN"

#: GFF3 column 3 type that delimits the annotation ranges.
_CDS_TYPE = "CDS"

#: Standard genetic code, codons read on the coding strand.
_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: Impact each consequence maps to.
_CONSEQUENCE_IMPACT = {
    "START_LOST": "HIGH",
    "STOP_GAINED": "HIGH",
    "STOP_LOST": "HIGH",
    "FRAMESHIFT": "HIGH",
    "SYNONYMOUS": "LOW",
    "MISSENSE": "MODERATE",
    "INFRAME_DELETION": "MODERATE",
    "INFRAME_INSERTION": "MODERATE",
}

_GVANN_HEADER = (
    '##INFO=<ID=GVANN,Number=.,Type=String,Description="Genome variant '
    'annotation: ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|'
    'AA_CHANGE; multiple transcripts are comma separated">'
)

_NO_TRANSCRIPT = ".|.|.|."


class AnnotationFormatError(ValueError):
    """Malformed GFF3 data or a VCF that already carries ``GVANN``."""


@dataclass(frozen=True)
class _CdsFragment:
    """One CDS GFF3 row: 1-based inclusive interval on one sequence."""

    seqid: str
    start: int
    end: int
    strand: str
    phase: int
    line_number: int


@dataclass(frozen=True)
class _CodonSite:
    """Everything known about one coding base without reading the reference."""

    cds_pos: int
    codon_number: int
    position_in_codon: int
    codon_indices: tuple[int, ...]


@dataclass(frozen=True)
class _Transcript:
    """A Parent group assembled into coding-strand codons."""

    parent: str
    strand: str
    seqid: str
    #: Genomic 0-based index -> coding site, for bases in complete codons.
    sites: dict[int, _CodonSite]
    #: Genomic 0-based index -> 1-based coding position, every CDS base.
    coding_positions: dict[int, int]
    #: Original fragments as (1-based start, 1-based end, source line).
    fragments: tuple[tuple[int, int, int], ...]


def _gff_error(source: str, line_number: int, message: str) -> AnnotationFormatError:
    return AnnotationFormatError(f"{source}:{line_number}: {message}")


def _decode_gff3(text: str) -> str:
    """Decode one percent-encoded GFF3 attribute value component."""
    output: list[str] = []
    index = 0
    while index < len(text):
        char = text[index]
        if char == "%":
            if index + 2 >= len(text):
                raise ValueError("truncated percent escape")
            try:
                output.append(chr(int(text[index + 1 : index + 3], 16)))
            except ValueError:
                raise ValueError(
                    f"invalid percent escape {text[index:index + 3]!r}"
                )
            index += 3
        else:
            output.append(char)
            index += 1
    return "".join(output)


def _parse_attributes(
    column: str, source: str, line_number: int
) -> dict[str, list[str]]:
    """Parse GFF3 column 9 into tag -> decoded values in column order."""
    attributes: dict[str, list[str]] = {}
    if column == "":
        return attributes
    for assignment in column.split(";"):
        if not assignment:
            raise _gff_error(source, line_number, "empty GFF3 attribute assignment")
        if "=" not in assignment:
            raise _gff_error(
                source, line_number, f"GFF3 attribute {assignment!r} lacks '='"
            )
        tag, raw_values = assignment.split("=", 1)
        if not tag:
            raise _gff_error(source, line_number, "empty GFF3 attribute tag")
        values: list[str] = []
        for raw_value in raw_values.split(","):
            try:
                values.append(_decode_gff3(raw_value))
            except ValueError as exc:
                raise _gff_error(
                    source, line_number, f"GFF3 attribute {tag!r}: {exc}"
                ) from None
        attributes.setdefault(tag, []).extend(values)
    return attributes


def _parse_cds_fragments(
    stream: Iterable[str], source: str
) -> dict[str, list[_CdsFragment]]:
    """Read a GFF3 stream and collect CDS fragments keyed by Parent."""
    fragments: dict[str, list[_CdsFragment]] = defaultdict(list)
    in_fasta = False
    for line_number, raw_line in enumerate(stream, start=1):
        line = raw_line.rstrip("\r\n")
        if line == "##FASTA":
            # An optional embedded FASTA section follows; it is not used.
            in_fasta = True
            continue
        if in_fasta or not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 9:
            raise _gff_error(
                source,
                line_number,
                f"expected 9 tab-separated GFF3 columns, got {len(columns)}",
            )
        if columns[2] != _CDS_TYPE:
            continue

        seqid = columns[0]
        start_text, end_text = columns[3], columns[4]
        strand, phase_text = columns[6], columns[7]
        attributes_text = columns[8]

        if not seqid:
            raise _gff_error(source, line_number, "CDS feature has an empty seqid")
        try:
            start = int(start_text)
        except ValueError:
            raise _gff_error(
                source, line_number, f"CDS start {start_text!r} is not an integer"
            ) from None
        try:
            end = int(end_text)
        except ValueError:
            raise _gff_error(
                source, line_number, f"CDS end {end_text!r} is not an integer"
            ) from None
        if start < 1 or end < 1:
            raise _gff_error(
                source,
                line_number,
                f"CDS coordinates {start}-{end} are not 1-based positive",
            )
        if start > end:
            raise _gff_error(
                source,
                line_number,
                f"CDS start {start} is greater than its end {end}",
            )
        if strand not in ("+", "-"):
            raise _gff_error(
                source, line_number, f"CDS strand {strand!r} is not '+' or '-'"
            )
        if phase_text not in ("0", "1", "2"):
            raise _gff_error(
                source,
                line_number,
                f"CDS phase {phase_text!r} is not one of '0', '1', '2'",
            )

        attributes = _parse_attributes(attributes_text, source, line_number)
        parents = attributes.get("Parent")
        if not parents:
            raise _gff_error(
                source, line_number, "CDS feature is missing a Parent attribute"
            )
        if any(parent == "" for parent in parents):
            raise _gff_error(
                source, line_number, "CDS feature has an empty Parent value"
            )
        phase = int(phase_text)
        for parent in parents:
            fragments[parent].append(
                _CdsFragment(
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    phase=phase,
                    line_number=line_number,
                )
            )
    return fragments


def _build_layout(
    parent: str, fragments: list[_CdsFragment], source: str
) -> list[int]:
    """Return CDS bases as genomic 0-based indices in coding order.

    Validates that every fragment shares one sequence and one strand and
    that each fragment's declared phase agrees with the assembled frame.
    """
    first = fragments[0]
    for fragment in fragments[1:]:
        if fragment.seqid != first.seqid:
            raise _gff_error(
                source,
                fragment.line_number,
                f"CDS for Parent {parent!r} spans multiple sequences "
                f"({first.seqid!r} and {fragment.seqid!r})",
            )
        if fragment.strand != first.strand:
            raise _gff_error(
                source,
                fragment.line_number,
                f"CDS for Parent {parent!r} mixes strands "
                f"({first.strand!r} and {fragment.strand!r})",
            )

    # Coding order: ascending coordinates on the plus strand, descending
    # on the minus strand (each minus-strand fragment read in reverse).
    ordered = sorted(fragments, key=lambda fragment: fragment.start)
    reverse = first.strand == "-"
    if reverse:
        ordered.reverse()

    layout: list[int] = []
    consumed = 0
    for fragment in ordered:
        # Phase is the number of coding bases missing at this fragment's
        # 5' edge before the next codon boundary: 1 -> 2, 2 -> 1.
        expected_phase = (3 - (consumed % 3)) % 3
        if fragment.phase != expected_phase:
            raise _gff_error(
                source,
                fragment.line_number,
                f"CDS phase {fragment.phase} for Parent {parent!r} contradicts "
                f"the assembled reading frame (expected {expected_phase})",
            )
        length = fragment.end - fragment.start + 1
        indices = range(fragment.start - 1, fragment.end)
        if reverse:
            indices = reversed(indices)
        layout.extend(indices)
        consumed += length

    return layout


def _build_transcript(
    parent: str, fragments: list[_CdsFragment], source: str
) -> _Transcript:
    layout = _build_layout(parent, fragments, source)
    sites: dict[int, _CodonSite] = {}
    for codon_number, start in enumerate(range(0, len(layout) - 2, 3), start=1):
        codon_indices = tuple(layout[start : start + 3])
        for position_in_codon, genomic_index in enumerate(codon_indices):
            sites[genomic_index] = _CodonSite(
                cds_pos=start + position_in_codon + 1,
                codon_number=codon_number,
                position_in_codon=position_in_codon,
                codon_indices=codon_indices,
            )

    coding_positions = {
        genomic_index: cds_pos
        for cds_pos, genomic_index in enumerate(layout, start=1)
    }

    return _Transcript(
        parent=parent,
        strand=fragments[0].strand,
        seqid=fragments[0].seqid,
        sites=sites,
        coding_positions=coding_positions,
        fragments=tuple(
            sorted(
                (fragment.start, fragment.end, fragment.line_number)
                for fragment in fragments
            )
        ),
    )


def _read_gff3(
    source: "str | os.PathLike[str] | io.TextIOBase",
) -> tuple[tuple[_Transcript, ...], str]:
    if isinstance(source, io.IOBase):
        stream: io.TextIOBase = source
        source_name = getattr(stream, "name", None)
        if not isinstance(source_name, str) or not source_name:
            source_name = "<stream>"
        close = False
    elif isinstance(source, (str, os.PathLike)):
        path = os.fspath(source)
        stream = open(path, "r", encoding="utf-8", newline="")
        source_name = path
        close = True
    else:
        raise TypeError("source must be a text path or a text stream")

    try:
        fragments = _parse_cds_fragments(stream, source_name)
        transcripts = [
            _build_transcript(parent, parent_fragments, source_name)
            for parent, parent_fragments in fragments.items()
        ]
        transcripts.sort(key=lambda transcript: transcript.parent)
        return tuple(transcripts), source_name
    finally:
        if close:
            stream.close()


def _build_reference(
    reference: "Mapping[str, str] | Iterable[SequenceRecord]",
) -> tuple[dict[str, str], set[str]]:
    """Map CHROM to upper-case sequence, tracking duplicated identifiers."""
    if isinstance(reference, Mapping):
        return {key: value.upper() for key, value in reference.items()}, set()
    sequences: dict[str, str] = {}
    duplicates: set[str] = set()
    for record in reference:
        if record.identifier in sequences:
            duplicates.add(record.identifier)
        else:
            sequences[record.identifier] = record.sequence.upper()
    return sequences, duplicates


def _is_acgt_snv(ref: str, alt: str) -> bool:
    return (
        len(ref) == 1
        and len(alt) == 1
        and ref in "ACGT"
        and alt in "ACGT"
        and ref != alt
    )


def _is_acgt(text: str) -> bool:
    return bool(text) and all(base in "ACGT" for base in text)


def _trim_alleles(ref: str, alt: str) -> tuple[int, int, str, str]:
    """Remove the longest common prefix, then suffix, of REF and ALT.

    Returns ``(prefix_length, suffix_length, ref_remainder, alt_remainder)``.
    The prefix is stripped first so the shared VCF anchor base stays on
    the left; the suffix is then measured on the remainders, never
    re-using a prefix base.
    """
    limit = min(len(ref), len(alt))
    prefix = 0
    while prefix < limit and ref[prefix] == alt[prefix]:
        prefix += 1
    ref_remainder = ref[prefix:]
    alt_remainder = alt[prefix:]
    suffix = 0
    while (
        suffix < len(ref_remainder)
        and suffix < len(alt_remainder)
        and ref_remainder[-1 - suffix] == alt_remainder[-1 - suffix]
    ):
        suffix += 1
    if suffix:
        ref_remainder = ref_remainder[:-suffix]
        alt_remainder = alt_remainder[:-suffix]
    return prefix, suffix, ref_remainder, alt_remainder


def _indel_kind(ref: str, alt: str) -> tuple[str, int, int] | None:
    """Classify a pure ACGT ALT as ``"deletion"``/``"insertion"``.

    Returns ``(kind, prefix_length, suffix_length)`` or ``None`` for a
    complex replacement (remainders on both sides), a no-op allele or an
    allele carrying non-ACGT symbols.
    """
    if not _is_acgt(ref) or not _is_acgt(alt):
        return None
    prefix, suffix, ref_remainder, alt_remainder = _trim_alleles(ref, alt)
    if ref_remainder and not alt_remainder:
        return "deletion", prefix, suffix
    if alt_remainder and not ref_remainder:
        return "insertion", prefix, suffix
    return None


def _classify(codon_number: int, old_codon: str, new_codon: str) -> tuple[str, str]:
    old_aa = _GENETIC_CODE[old_codon]
    new_aa = _GENETIC_CODE[new_codon]
    if codon_number == 1 and old_aa == "M" and new_aa != "M":
        consequence = "START_LOST"
    elif old_aa != "*" and new_aa == "*":
        consequence = "STOP_GAINED"
    elif old_aa == "*" and new_aa != "*":
        consequence = "STOP_LOST"
    elif old_aa == new_aa:
        consequence = "SYNONYMOUS"
    else:
        consequence = "MISSENSE"
    return consequence, _CONSEQUENCE_IMPACT[consequence]


def _fragment_containing(
    transcript: _Transcript, start: int, end: int
) -> tuple[int, int, int] | None:
    """The unique CDS fragment covering the 1-based inclusive interval.

    ``None`` when no fragment covers it or more than one does (overlapping
    fragments or an interval straddling a fragment boundary).
    """
    matches = [
        fragment
        for fragment in transcript.fragments
        if fragment[0] <= start and end <= fragment[1]
    ]
    return matches[0] if len(matches) == 1 else None


def _annotate_indel_transcript(
    record: VcfRecord,
    alt: str,
    kind: str,
    prefix: int,
    suffix: int,
    transcript: _Transcript,
) -> str | None:
    """Annotate one pure ACGT indel against one transcript.

    Returns the transcript's GVANN item (``FRAMESHIFT``/``HIGH`` or
    ``INFRAME_*``/``MODERATE`` on success, ``UNSUPPORTED``/``MODIFIER``
    when the event touches the CDS but cannot be judged), or ``None``
    when the event never touches this transcript's CDS.
    """
    unsupported = f"{alt}|UNSUPPORTED|MODIFIER|{transcript.parent}|.|.|."

    if kind == "deletion":
        # The deleted bases follow the shared anchor (record POS plus the
        # stripped common prefix) and exclude the trailing common base.
        length = len(record.ref) - prefix - suffix
        start_index = record.pos - 1 + prefix
        positions = [
            transcript.coding_positions.get(index)
            for index in range(start_index, start_index + length)
        ]
        if all(position is None for position in positions):
            return None
        if any(position is None for position in positions):
            # Part of the deleted interval is outside every CDS fragment.
            return unsupported
        # The whole deletion must fall inside one CDS fragment.
        if _fragment_containing(
            transcript, start_index + 1, start_index + length
        ) is None:
            return unsupported
        # ...and the deleted bases must be consecutive in coding order.
        first = min(positions)
        if sorted(positions) != list(range(first, first + length)):
            return unsupported
        cds_pos = first
        consequence = "FRAMESHIFT" if length % 3 else "INFRAME_DELETION"
    else:
        # The inserted bases sit between the last shared-prefix base and
        # the next genomic base (which is the first suffix base when a
        # shared suffix exists, since prefix plus suffix exhaust REF):
        # genomic 0-based indices L and R.
        length = len(alt) - prefix - suffix
        left_index = record.pos - 2 + prefix
        right_index = record.pos - 1 + prefix
        left_position = transcript.coding_positions.get(left_index)
        right_position = transcript.coding_positions.get(right_index)
        if left_position is None and right_position is None:
            return None
        if left_position is None or right_position is None:
            # One coding anchor base is missing (outside the CDS).
            return unsupported
        # Both flanking bases must live in the same CDS fragment.
        if _fragment_containing(
            transcript, left_index + 1, right_index + 1
        ) is None:
            return unsupported
        # They must be adjacent in the assembled coding sequence; the
        # reported position is the downstream one in coding direction,
        # which is strand independent.
        if abs(right_position - left_position) != 1:
            return unsupported
        cds_pos = max(left_position, right_position)
        consequence = "FRAMESHIFT" if length % 3 else "INFRAME_INSERTION"

    impact = "HIGH" if consequence == "FRAMESHIFT" else "MODERATE"
    return f"{alt}|{consequence}|{impact}|{transcript.parent}|{cds_pos}|.|."


def _annotate_alt(
    record: VcfRecord,
    alt: str,
    transcripts: tuple[_Transcript, ...],
    chromosome: str,
) -> str:
    indel = _indel_kind(record.ref, alt)
    if indel is not None:
        kind, prefix, suffix = indel
        items: list[str] = []
        # transcripts are already ordered by Parent, so hits stay lexicographic.
        for transcript in transcripts:
            if transcript.seqid != record.chrom:
                continue
            item = _annotate_indel_transcript(
                record, alt, kind, prefix, suffix, transcript
            )
            if item is not None:
                items.append(item)
        if not items:
            return f"{alt}|NON_CODING|MODIFIER|{_NO_TRANSCRIPT}"
        return ",".join(items)

    if not _is_acgt_snv(record.ref, alt):
        return f"{alt}|UNSUPPORTED|MODIFIER|{_NO_TRANSCRIPT}"

    genomic_index = record.pos - 1
    items = []
    # transcripts are already ordered by Parent, so hits stay lexicographic.
    for transcript in transcripts:
        if transcript.seqid != record.chrom:
            continue
        site = transcript.sites.get(genomic_index)
        if site is None:
            continue

        old_bases: list[str] = []
        for index in site.codon_indices:
            base = chromosome[index]
            if transcript.strand == "-":
                base = base.translate(_COMPLEMENT)
            old_bases.append(base)
        if any(base not in "ACGT" for base in old_bases):
            # Ambiguous coding context: the standard code cannot judge it.
            items.append(
                f"{alt}|UNSUPPORTED|MODIFIER|{transcript.parent}|"
                f"{site.cds_pos}|.|."
            )
            continue

        old_codon = "".join(old_bases)
        alt_coding_base = (
            alt.translate(_COMPLEMENT) if transcript.strand == "-" else alt
        )
        new_bases = list(old_codon)
        new_bases[site.position_in_codon] = alt_coding_base
        new_codon = "".join(new_bases)

        consequence, impact = _classify(site.codon_number, old_codon, new_codon)
        old_aa = _GENETIC_CODE[old_codon]
        new_aa = _GENETIC_CODE[new_codon]
        items.append(
            "|".join(
                (
                    alt,
                    consequence,
                    impact,
                    transcript.parent,
                    str(site.cds_pos),
                    f"{old_codon}>{new_codon}",
                    f"{old_aa}>{new_aa}",
                )
            )
        )

    if not items:
        return f"{alt}|NON_CODING|MODIFIER|{_NO_TRANSCRIPT}"
    return ",".join(items)


def _meta_declares_gvann(line: str) -> bool:
    if not line.startswith("##"):
        return False
    body = line[2:]
    if body == GVANN_ID or body.startswith(GVANN_ID + "="):
        return True
    # Structured meta-information, e.g. ##INFO=<ID=GVANN,...>.
    return f"<ID={GVANN_ID}," in body or f"<ID={GVANN_ID}>" in body


def _info_has_gvann(info: str) -> bool:
    if not info:
        return False
    for item in info.split(";"):
        if item.split("=", 1)[0] == GVANN_ID:
            return True
    return False


def _check_existing_gvann(document: VcfFile) -> None:
    for line in document.header.meta_lines:
        if _meta_declares_gvann(line):
            raise AnnotationFormatError(
                f"{GVANN_ID} is already declared in the VCF meta-information"
            )
    for record in document.records:
        if _info_has_gvann(record.info):
            source = document.source
            prefix = f"{source}:{record.line_number}: " if source and record.line_number else ""
            raise AnnotationFormatError(
                f"{prefix}{GVANN_ID} is already present in a VCF record"
            )


def annotate_vcf(
    document: VcfFile,
    reference: "Mapping[str, str] | Iterable[SequenceRecord]",
    features: "str | os.PathLike[str] | io.TextIOBase",
) -> VcfFile:
    """Annotate every ALT allele of *document* with a ``GVANN`` INFO item.

    *reference* is a mapping of CHROM to upper-case sequence or an
    iterable of :class:`~genome_variant.sequence_io.SequenceRecord`;
    *features* is a GFF3 text path or text stream whose CDS features
    define the annotated ranges, grouped into transcripts by ``Parent``.

    Record order, ALT order, sample columns and every existing field are
    preserved; one unique ``##INFO`` declaration for ``GVANN`` is
    appended to the meta-information and a ``GVANN`` item is appended to
    the end of each record's INFO column (starting the column when INFO
    was missing).  Coding SNV consequences are ``START_LOST``,
    ``STOP_GAINED``, ``STOP_LOST``, ``SYNONYMOUS`` and ``MISSENSE`` with
    impacts ``HIGH``, ``HIGH``, ``HIGH``, ``LOW`` and ``MODERATE``.  A
    pure ``A``/``C``/``G``/``T`` insertion or deletion whose longest
    common REF/ALT prefix and suffix trim to a single event is judged per
    Parent: wholly within one CDS fragment, with the deleted bases (or
    the insertion's two flanking bases) consecutive in the assembled
    coding sequence, it is ``FRAMESHIFT``/``HIGH`` when the length change
    is not a multiple of three and otherwise
    ``INFRAME_DELETION``/``INFRAME_INSERTION`` with ``MODERATE``;
    ``CDS_POS`` is the first affected coding position (the downstream
    base for an insertion) and ``CODON_CHANGE``/``AA_CHANGE`` are ``.``.
    An indel that touches a CDS across fragment boundaries, over
    non-coding bases or without a coding anchor on both sides gets
    ``UNSUPPORTED``/``MODIFIER`` for that transcript.  Transcripts are
    listed by lexicographic Parent.  Variants outside every CDS get
    ``NON_CODING``/``MODIFIER``; complex replacements (remainders on both
    sides) and ambiguous, symbolic, break-end, spanning or missing ALTs
    get ``UNSUPPORTED``/``MODIFIER``.

    Every record's REF is checked against the reference first.  Raises
    :class:`~genome_variant.vcf.ReferenceMismatchError` for an
    absent/duplicated CHROM or a REF mismatch,
    :class:`~genome_variant.vcf.VcfFormatError` for an out-of-bounds
    record and :class:`AnnotationFormatError` for malformed GFF3 data or
    a document that already declares or uses ``GVANN``.
    """
    transcripts, features_source = _read_gff3(features)
    _check_existing_gvann(document)
    sequences, duplicates = _build_reference(reference)

    # A CDS interval cannot run past the end of the sequence it names.
    for transcript in transcripts:
        chromosome = sequences.get(transcript.seqid)
        if chromosome is None:
            continue
        length = len(chromosome)
        for start, end, line_number in transcript.fragments:
            if end > length:
                raise AnnotationFormatError(
                    f"{features_source}:{line_number}: CDS coordinates "
                    f"{start}-{end} run past the end of sequence "
                    f"{transcript.seqid!r} (length {length})"
                )

    annotated_records: list[VcfRecord] = []
    for record in document.records:
        if record.chrom in duplicates:
            raise ReferenceMismatchError(
                f"CHROM {record.chrom} POS {record.pos}: CHROM occurs more than "
                "once in the reference"
            )
        try:
            chromosome = sequences[record.chrom]
        except KeyError:
            raise ReferenceMismatchError(
                f"CHROM {record.chrom} POS {record.pos}: no matching record in "
                "the reference"
            ) from None

        end = record.pos - 1 + len(record.ref)
        if record.pos < 1 or end > len(chromosome):
            prefix = (
                f"{document.source}:{record.line_number}: "
                if document.source and record.line_number
                else ""
            )
            raise VcfFormatError(
                f"{prefix}CHROM {record.chrom} POS {record.pos}: REF of length "
                f"{len(record.ref)} is out of bounds for a reference of length "
                f"{len(chromosome)}"
            )
        expected = chromosome[record.pos - 1 : end]
        if expected != record.ref:
            raise ReferenceMismatchError(
                f"CHROM {record.chrom} POS {record.pos}: REF does not match the "
                f"reference: expected {expected!r}, found {record.ref!r}"
            )

        items = [
            _annotate_alt(record, alt, transcripts, chromosome)
            for alt in record.alt
        ]
        gvann = f"{GVANN_ID}=" + ",".join(items)
        # The lone "." is the VCF missing-value placeholder, i.e. no INFO:
        # GVANN then starts the column; otherwise it is appended last.
        if record.info and record.info != ".":
            info = f"{record.info};{gvann}"
        else:
            info = gvann
        annotated_records.append(replace(record, info=info))

    header = VcfHeader(
        meta_lines=document.header.meta_lines + (_GVANN_HEADER,),
        samples=document.header.samples,
    )
    return VcfFile(header, tuple(annotated_records), document.source)
