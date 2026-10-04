"""GFF3-driven variant annotation for parsed VCF documents.

Public API:

- :class:`AnnotationFormatError` — malformed GFF3 input (bad CDS rows,
  coordinates or phase, missing Parent, a strand other than ``+``/``-``,
  CDS fragments spanning several sequences or contradictory reading
  frames) or a VCF that already declares or uses the ``GVANN`` key.
- :func:`annotate_vcf` — append one ``GVANN`` INFO item per ALT allele to
  every record of a :class:`~genome_variant.vcf.VcfFile`.

CDS features are grouped by their GFF3 ``Parent`` attribute into
transcripts.  Attribute values are split on unescaped commas and each
value's ``%HH`` escapes are then restored as UTF-8 bytes (unescaped
non-ASCII text is kept as-is), so ``Parent=tx%E5%9F%BA%E5%9B%A0`` and
``Parent=tx基因`` name the same transcript.  In the ``GVANN`` TRANSCRIPT
field, characters that would break the VCF INFO or GVANN structure
(``%``, ``,``, ``;``, ``|``, ``=`` and ASCII whitespace or control
characters) are percent-encoded by UTF-8 byte with upper-case hex;
other non-ASCII characters are emitted directly.  The standard genetic
code is used.  Single-base ``A``/``C``/
``G``/``T`` substitutions inside a CDS are classified by their codon
change.  Ordinary ``A``/``C``/``G``/``T`` insertions and deletions are
graded by their length change: a non-multiple of three is ``FRAMESHIFT``
(``HIGH``) and an in-frame gain or loss of codons is
``INFRAME_INSERTION``/``INFRAME_DELETION`` (``MODERATE``), provided the
event lies inside a single CDS fragment of the transcript with the
coding anchors (and, for deletions, every deleted base) intact; events
that touch a CDS across fragments, through non-coding bases or without a
coding anchor on both sides receive ``UNSUPPORTED``/``MODIFIER`` for that
transcript.  Every other ALT (complex replacements, ambiguous or
symbolic alleles) receives ``UNSUPPORTED``/``MODIFIER`` and records
without a CDS hit receive ``NON_CODING``/``MODIFIER``.
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
    #: Original fragments as (1-based start, 1-based end, source line).
    fragments: tuple[tuple[int, int, int], ...]
    #: Genomic 0-based index -> 1-based coding position, every CDS base.
    coding_positions: dict[int, int]


def _gff_error(source: str, line_number: int, message: str) -> AnnotationFormatError:
    return AnnotationFormatError(f"{source}:{line_number}: {message}")


_HEX_DIGITS = frozenset("0123456789abcdefABCDEF")


def _decode_gff3(text: str) -> str:
    """Decode one percent-encoded GFF3 attribute value component.

    Escapes are UTF-8 bytes, not Unicode code points: the ``%HH``
    sequences (consecutive or scattered) contribute their bytes and the
    unescaped text — already decoded UTF-8 — contributes its own
    encoding, then the assembled bytes are decoded as UTF-8.  Hex digits
    are case-insensitive.  Truncated escapes, non-hex escapes and escape
    bytes that do not form valid UTF-8 raise :class:`ValueError`.
    """
    buffer = bytearray()
    index = 0
    while index < len(text):
        char = text[index]
        if char == "%":
            if index + 2 >= len(text):
                raise ValueError("truncated percent escape")
            digits = text[index + 1 : index + 3]
            if any(digit not in _HEX_DIGITS for digit in digits):
                raise ValueError(
                    f"invalid percent escape {text[index:index + 3]!r}"
                )
            buffer.append(int(digits, 16))
            index += 3
        else:
            buffer.extend(char.encode("utf-8"))
            index += 1
    try:
        return bytes(buffer).decode("utf-8")
    except UnicodeDecodeError:
        raise ValueError(
            "percent escapes do not form valid UTF-8"
        ) from None


#: Characters that would break the VCF INFO column or the GVANN
#: ``|``/``,``/``=`` structure, and therefore cannot appear literally in
#: a TRANSCRIPT field.
_GVANN_ESCAPED = frozenset("%,;|=")


def _encode_transcript(name: str) -> str:
    """Render a transcript name safe for the GVANN TRANSCRIPT field.

    ``%``, ``,``, ``;``, ``|``, ``=`` and ASCII whitespace or control
    characters are percent-encoded by UTF-8 byte with upper-case hex;
    every other character, including non-ASCII ones, is emitted as-is.
    """
    output: list[str] = []
    for char in name:
        code = ord(char)
        if char in _GVANN_ESCAPED or code <= 0x20 or code == 0x7F:
            output.extend(f"%{byte:02X}" for byte in char.encode("utf-8"))
        else:
            output.append(char)
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
    coding_positions: dict[int, int] = {}
    for cds_pos, genomic_index in enumerate(layout, start=1):
        coding_positions[genomic_index] = cds_pos
    for codon_number, start in enumerate(range(0, len(layout) - 2, 3), start=1):
        codon_indices = tuple(layout[start : start + 3])
        for position_in_codon, genomic_index in enumerate(codon_indices):
            sites[genomic_index] = _CodonSite(
                cds_pos=start + position_in_codon + 1,
                codon_number=codon_number,
                position_in_codon=position_in_codon,
                codon_indices=codon_indices,
            )
    return _Transcript(
        parent=parent,
        strand=fragments[0].strand,
        seqid=fragments[0].seqid,
        sites=sites,
        fragments=tuple(
            sorted(
                (fragment.start, fragment.end, fragment.line_number)
                for fragment in fragments
            )
        ),
        coding_positions=coding_positions,
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


def _is_plain_allele(allele: str) -> bool:
    """Whether *allele* is an ordinary sequence allele without symbols.

    Mirrors the VCF reader's notion: symbolic (``<DEL>``), break-end
    (brackets or a ``.`` join), spanning (``*``) and missing (``.``)
    alleles are not plain.
    """
    return (
        bool(allele)
        and allele not in (".", "*")
        and not any(char in "<>[]*." for char in allele)
    )


def _is_acgt_text(text: str) -> bool:
    return bool(text) and all(char in "ACGT" for char in text)


def _trim_alleles(ref: str, alt: str) -> tuple[int, int]:
    """Measure REF/ALT's longest common prefix and suffix.

    Returns ``(prefix, suffix)`` shared base counts; the suffix
    comparison never reuses a base consumed by the prefix.
    """
    prefix = 0
    while prefix < len(ref) and prefix < len(alt) and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while (
        suffix < len(ref) - prefix
        and suffix < len(alt) - prefix
        and ref[len(ref) - 1 - suffix] == alt[len(alt) - 1 - suffix]
    ):
        suffix += 1
    return prefix, suffix


def _containing_fragment(
    transcript: _Transcript, genomic_index: int
) -> int | None:
    """Index of the CDS fragment covering *genomic_index*, if any."""
    for number, (start, end, _line_number) in enumerate(transcript.fragments):
        if start - 1 <= genomic_index <= end - 1:
            return number
    return None


def _indel_transcript(
    is_deletion: bool,
    length: int,
    left_anchor: int,
    right_anchor: int,
    deleted_indices: tuple[int, ...],
    chromosome_length: int,
    transcript: _Transcript,
) -> tuple[str, str, int] | None:
    """Grade one pure indel against one transcript.

    Returns ``(consequence, impact, cds_pos)`` — consequence is
    ``UNSUPPORTED`` when the event touches the CDS but cannot be judged —
    or ``None`` when the event does not touch the transcript.

    Anchors are the genomic bases immediately flanking the event, found
    from coordinates alone (a minimized indel carries only one flank in
    the allele); an index outside the chromosome has no anchor at all.
    """
    coding = transcript.coding_positions

    def coding_at(index: int) -> int | None:
        if 0 <= index < chromosome_length:
            return coding.get(index)
        return None

    left_cds = coding_at(left_anchor)
    right_cds = coding_at(right_anchor)
    deleted_cds = [
        coding[index]
        for index in deleted_indices
        if index in coding
    ]

    if not deleted_cds and left_cds is None and right_cds is None:
        return None

    if is_deletion:
        cds_pos = min(
            deleted_cds,
            default=(
                left_cds if left_cds is not None else right_cds
            ),
        )

        # The deleted bases must all be coding and consecutive in the
        # assembled coding sequence (coding positions work on either
        # strand), and the deleted stretch together with both flanking
        # anchors must lie inside one single CDS fragment.  A missing or
        # non-coding anchor (fragment edge, intron, chromosome end) makes
        # the event unjudgeable for this transcript.
        anchor_positions = [
            value
            for value in (left_cds, right_cds)
            if value is not None
        ]
        all_positions = anchor_positions + deleted_cds
        expected_count = len(deleted_indices) + 2
        contiguous = (
            len(all_positions) == expected_count
            and max(all_positions) - min(all_positions) + 1
            == expected_count
        )
        fragments = {
            _containing_fragment(transcript, index)
            for index in (left_anchor, right_anchor, *deleted_indices)
        }
        eligible = (
            len(deleted_cds) == len(deleted_indices)
            and left_cds is not None
            and right_cds is not None
            and contiguous
            and len(fragments) == 1
            and next(iter(fragments)) is not None
        )
        if not eligible:
            return ("UNSUPPORTED", "MODIFIER", cds_pos)
        consequence = (
            "FRAMESHIFT" if length % 3 else "INFRAME_DELETION"
        )
        return (consequence, _CONSEQUENCE_IMPACT[consequence], cds_pos)

    # Insertion: the boundary must lie between two coding bases adjacent
    # in coding direction (coding order descends with genomic coordinates
    # on the minus strand), and both flanks must belong to one CDS
    # fragment.  CDS_POS names the first coding base to the right of the
    # boundary in coding direction, the larger of the two positions.
    cds_pos = (
        max(left_cds, right_cds)
        if left_cds is not None and right_cds is not None
        else (left_cds if left_cds is not None else right_cds)
    )
    adjacent = (
        left_cds is not None
        and right_cds is not None
        and abs(right_cds - left_cds) == 1
    )
    left_fragment = _containing_fragment(transcript, left_anchor)
    right_fragment = _containing_fragment(transcript, right_anchor)
    eligible = (
        adjacent
        and left_fragment is not None
        and left_fragment == right_fragment
    )
    if not eligible:
        return ("UNSUPPORTED", "MODIFIER", cds_pos)
    consequence = (
        "FRAMESHIFT" if length % 3 else "INFRAME_INSERTION"
    )
    return (consequence, _CONSEQUENCE_IMPACT[consequence], cds_pos)


def _annotate_indel(
    record: VcfRecord,
    alt: str,
    transcripts: tuple[_Transcript, ...],
    chromosome: str,
) -> str | None:
    """Annotate one pure ACGT indel ALT, or return ``None`` if not one.

    ``None`` covers symbolic, break-end, spanning and missing ALTs,
    ambiguous alleles and complex replacements (remainders on both
    sides), which stay globally ``UNSUPPORTED``.
    """
    if not (
        _is_plain_allele(record.ref)
        and _is_plain_allele(alt)
        and _is_acgt_text(record.ref)
        and _is_acgt_text(alt)
    ):
        return None

    prefix, suffix = _trim_alleles(record.ref, alt)
    ref_remainder = record.ref[prefix : len(record.ref) - suffix]
    alt_remainder = alt[prefix : len(alt) - suffix]
    if ref_remainder and alt_remainder:
        # Bases changed on both sides: a complex replacement, not an indel.
        return None
    if not ref_remainder and not alt_remainder:
        return None
    is_deletion = bool(ref_remainder)
    length = len(ref_remainder) if is_deletion else len(alt_remainder)

    # Flanks come from genomic coordinates, not from the allele's
    # retained bases, so minimized and right-anchored representations of
    # one event are judged identically.  Both indices may point outside
    # the chromosome at its edges, which means there is no anchor.
    if is_deletion:
        first = record.pos - 1 + prefix
        deleted_indices = tuple(range(first, first + length))
        left_anchor = first - 1
        right_anchor = first + length
    else:
        deleted_indices = ()
        gap = record.pos - 1 + prefix
        left_anchor = gap - 1
        right_anchor = gap

    items: list[str] = []
    # transcripts are already ordered by Parent, so hits stay lexicographic.
    for transcript in transcripts:
        if transcript.seqid != record.chrom:
            continue
        result = _indel_transcript(
            is_deletion,
            length,
            left_anchor,
            right_anchor,
            deleted_indices,
            len(chromosome),
            transcript,
        )
        if result is None:
            continue
        consequence, impact, cds_pos = result
        items.append(
            f"{alt}|{consequence}|{impact}|"
            f"{_encode_transcript(transcript.parent)}|{cds_pos}|.|."
        )

    if not items:
        return f"{alt}|NON_CODING|MODIFIER|{_NO_TRANSCRIPT}"
    return ",".join(items)


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


def _annotate_alt(
    record: VcfRecord,
    alt: str,
    transcripts: tuple[_Transcript, ...],
    chromosome: str,
) -> str:
    if not _is_acgt_snv(record.ref, alt):
        indel = _annotate_indel(record, alt, transcripts, chromosome)
        return indel if indel is not None else (
            f"{alt}|UNSUPPORTED|MODIFIER|{_NO_TRANSCRIPT}"
        )

    genomic_index = record.pos - 1
    items: list[str] = []
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
                f"{alt}|UNSUPPORTED|MODIFIER|"
                f"{_encode_transcript(transcript.parent)}|"
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
                    _encode_transcript(transcript.parent),
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
    define the annotated ranges, grouped into transcripts by ``Parent``
    after percent-decoding attribute values as UTF-8 bytes.

    Record order, ALT order, sample columns and every existing field are
    preserved; one unique ``##INFO`` declaration for ``GVANN`` is
    appended to the meta-information and a ``GVANN`` item is appended to
    the end of each record's INFO column (starting the column when INFO
    was missing).  Coding consequences are ``START_LOST``,
    ``STOP_GAINED``, ``STOP_LOST``, ``SYNONYMOUS`` and ``MISSENSE`` with
    impacts ``HIGH``, ``HIGH``, ``HIGH``, ``LOW`` and ``MODERATE``;
    transcripts are listed by lexicographic Parent.  Pure ``ACGT``
    insertions and deletions are graded per transcript after trimming the
    longest common prefix and suffix of REF and ALT: a length change not
    divisible by three is ``FRAMESHIFT``/``HIGH`` and an in-frame change
    is ``INFRAME_INSERTION`` or ``INFRAME_DELETION``/``MODERATE``, when
    the event lies inside one CDS fragment with intact coding anchors;
    events that touch a CDS across fragments, through non-coding bases or
    without a coding anchor on both sides get ``UNSUPPORTED``/``MODIFIER``
    for that transcript.  SNVs and pure indels outside every CDS get
    ``NON_CODING``/``MODIFIER``; complex replacements and non-SNV,
    ambiguous, symbolic, break-end, spanning or missing ALTs get
    ``UNSUPPORTED``/``MODIFIER``.

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
