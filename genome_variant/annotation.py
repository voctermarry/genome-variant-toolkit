"""Variant consequence annotation of VCF records against a GFF3 annotation.

Public API:

- :class:`AnnotationFormatError` — malformed GFF3 (bad CDS fields,
  coordinates, phase, strand or parentage; fragments spanning sequences
  or an inconsistent reading frame) or a VCF that already declares or
  uses the ``GVANN`` INFO key.
- :func:`annotate_vcf` — annotate every record of a parsed
  :class:`~genome_variant.vcf.VcfFile` and return a new document.

GFF3 ``CDS`` features sharing a ``Parent`` form one transcript.  Their
fragments are concatenated in transcription order (ascending genomic
coordinates on the ``+`` strand, descending on the ``-`` strand) and the
GFF3 ``phase`` column fixes the reading frame: phase 0/1/2 means 0/2/1
bases precede the first codon base at the 5' end of the feature, and
every later fragment's phase must match the cumulative frame.  Only
single-base ``A``/``C``/``G``/``T`` substitutions of a reference base
inside a coding region are translated with the standard genetic code;
every other allele is reported as ``UNSUPPORTED`` (or ``NON_CODING``
when it misses every CDS) while the record itself is preserved apart
from the appended ``GVANN`` INFO item.
"""

from __future__ import annotations

import io
import os
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace

from .sequence_io import SequenceRecord
from .vcf import (
    ReferenceMismatchError,
    VcfFile,
    VcfHeader,
    VcfRecord,
    _build_reference,
    _is_plain_allele,
    _record_error,
)

__all__ = ["AnnotationFormatError", "annotate_vcf"]

_INFO_KEY = "GVANN"
_META_LINE = (
    '##INFO=<ID=GVANN,Number=.,Type=String,Description="Genomic Variant '
    'Annotation: ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|AA_CHANGE">'
)

_IMPACT = {
    "START_LOST": "HIGH",
    "STOP_GAINED": "HIGH",
    "STOP_LOST": "HIGH",
    "SYNONYMOUS": "LOW",
    "MISSENSE": "MODERATE",
    "NON_CODING": "MODIFIER",
    "UNSUPPORTED": "MODIFIER",
}
_DNA_BASES = frozenset("ACGT")

# Standard genetic code, encoded-strand codons; '*' is a stop.
_CODON_TABLE: dict[str, str] = {}


def _build_codon_table() -> None:
    blocks = (
        ("TTT", "F", "TTC", "F", "TTA", "L", "TTG", "L"),
        ("CTT", "L", "CTC", "L", "CTA", "L", "CTG", "L"),
        ("ATT", "I", "ATC", "I", "ATA", "I", "ATG", "M"),
        ("GTT", "V", "GTC", "V", "GTA", "V", "GTG", "V"),
        ("TCT", "S", "TCC", "S", "TCA", "S", "TCG", "S"),
        ("CCT", "P", "CCC", "P", "CCA", "P", "CCG", "P"),
        ("ACT", "T", "ACC", "T", "ACA", "T", "ACG", "T"),
        ("GCT", "A", "GCC", "A", "GCA", "A", "GCG", "A"),
        ("TAT", "Y", "TAC", "Y", "TAA", "*", "TAG", "*"),
        ("CAT", "H", "CAC", "H", "CAA", "Q", "CAG", "Q"),
        ("AAT", "N", "AAC", "N", "AAA", "K", "AAG", "K"),
        ("GAT", "D", "GAC", "D", "GAA", "E", "GAG", "E"),
        ("TGT", "C", "TGC", "C", "TGA", "*", "TGG", "W"),
        ("CGT", "R", "CGC", "R", "CGA", "R", "CGG", "R"),
        ("AGT", "S", "AGC", "S", "AGA", "R", "AGG", "R"),
        ("GGT", "G", "GGC", "G", "GGA", "G", "GGG", "G"),
    )
    for block in blocks:
        for index in range(0, len(block), 2):
            _CODON_TABLE[block[index]] = block[index + 1]


_build_codon_table()

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class AnnotationFormatError(ValueError):
    """Malformed GFF3 annotation or an illegal GVANN usage in the VCF."""


@dataclass(frozen=True)
class _CdsFragment:
    seqid: str
    start: int  # 1-based inclusive
    end: int
    phase: int
    strand: str
    line_number: int


@dataclass(frozen=True)
class _Transcript:
    name: str
    seqid: str
    strand: str
    cds: str  # all CDS bases, encoded strand, 5' -> 3'
    leader: int  # untranslated bases at the 5' end fixed by the first phase
    # 1-based genomic coordinate (forward reference) of every CDS base,
    # in 5' -> 3' transcription order.
    genome_positions: tuple[int, ...]

    @property
    def coding_positions(self) -> tuple[int, ...]:
        return self.genome_positions[self.leader :]


def _feature_error(source: str, line_number: int, message: str) -> AnnotationFormatError:
    return AnnotationFormatError(f"{source}:{line_number}: {message}")


def _record_annotation_error(
    source: str, record: VcfRecord, message: str
) -> AnnotationFormatError:
    if source and record.line_number:
        return AnnotationFormatError(f"{source}:{record.line_number}: {message}")
    if record.line_number:
        return AnnotationFormatError(f"line {record.line_number}: {message}")
    return AnnotationFormatError(message)


def _decode_attribute_value(value: str) -> str:
    """Decode one GFF3 percent-encoded attribute value."""
    decoded: list[str] = []
    index = 0
    while index < len(value):
        char = value[index]
        if char == "%" and index + 2 < len(value):
            try:
                decoded.append(chr(int(value[index + 1 : index + 3], 16)))
            except ValueError:
                decoded.append(value[index : index + 3])
            index += 3
        else:
            decoded.append(char)
            index += 1
    return "".join(decoded)


def _parent_attribute(attr_text: str) -> str | None:
    if not attr_text or attr_text == ".":
        return None
    for item in attr_text.split(";"):
        if item.startswith("Parent="):
            return item[len("Parent=") :]
    return None


def _parse_gff(
    stream: Iterable[str], source_name: str
) -> dict[str, list[_CdsFragment]]:
    """Collect CDS fragments keyed by their Parent transcript."""
    fragments: dict[str, list[_CdsFragment]] = defaultdict(list)
    for number, raw_line in enumerate(stream, start=1):
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            # Only CDS features participate in annotation; a CDS row
            # without nine columns is malformed, other lines are ignored.
            if len(fields) >= 3 and fields[2] == "CDS":
                raise _feature_error(
                    source_name,
                    number,
                    f"expected 9 tab-separated GFF3 columns, got {len(fields)}",
                )
            continue
        if fields[2] != "CDS":
            continue

        seqid = fields[0]
        start_text, end_text = fields[3], fields[4]
        strand, phase_text = fields[6], fields[7]

        if not seqid:
            raise _feature_error(source_name, number, "CDS feature has an empty seqid")
        if strand not in ("+", "-"):
            raise _feature_error(
                source_name,
                number,
                f"CDS feature strand must be '+' or '-', got {strand!r}",
            )
        if not start_text.isdigit() or not end_text.isdigit():
            raise _feature_error(
                source_name,
                number,
                "CDS feature coordinates must be positive integers, got "
                f"{start_text!r}-{end_text!r}",
            )
        start, end = int(start_text), int(end_text)
        if start < 1 or end < start:
            raise _feature_error(
                source_name, number, f"CDS feature coordinates {start}-{end} are illegal"
            )
        if phase_text not in ("0", "1", "2"):
            raise _feature_error(
                source_name,
                number,
                f"CDS feature phase must be 0, 1 or 2, got {phase_text!r}",
            )

        parent_value = _parent_attribute(fields[8])
        if parent_value is None or parent_value == "":
            raise _feature_error(
                source_name, number, "CDS feature is missing a Parent attribute"
            )
        if "," in parent_value:
            raise _feature_error(
                source_name,
                number,
                f"CDS feature Parent {parent_value!r} lists more than one parent",
            )
        parent = _decode_attribute_value(parent_value)
        if not parent:
            raise _feature_error(
                source_name, number, "CDS feature has an empty Parent attribute"
            )

        fragments[parent].append(
            _CdsFragment(seqid, start, end, int(phase_text), strand, number)
        )
    return fragments


def _complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)


def _build_transcript(
    name: str,
    fragments: Sequence[_CdsFragment],
    sequences: Mapping[str, str],
    duplicates: set[str],
    source_name: str,
) -> _Transcript:
    ascending = sorted(fragments, key=lambda fragment: (fragment.start, fragment.end))
    first = ascending[0]
    seqid, strand = first.seqid, first.strand

    for fragment in ascending:
        if fragment.seqid != seqid:
            raise _feature_error(
                source_name,
                fragment.line_number,
                f"CDS fragments of transcript {name!r} span multiple sequences "
                f"({seqid!r} and {fragment.seqid!r})",
            )
        if fragment.strand != strand:
            raise _feature_error(
                source_name,
                fragment.line_number,
                f"CDS fragments of transcript {name!r} disagree on the strand",
            )
    for earlier, later in zip(ascending, ascending[1:]):
        if later.start <= earlier.end:
            raise _feature_error(
                source_name,
                later.line_number,
                f"CDS fragments of transcript {name!r} overlap at position "
                f"{later.start}",
            )

    if seqid in duplicates:
        raise ReferenceMismatchError(
            f"CHROM {seqid}: CHROM occurs more than once in the reference"
        )
    try:
        chromosome = sequences[seqid]
    except KeyError:
        raise ReferenceMismatchError(
            f"CHROM {seqid}: CDS feature references a sequence absent from the "
            "reference"
        ) from None

    for fragment in ascending:
        if fragment.end > len(chromosome):
            raise _feature_error(
                source_name,
                fragment.line_number,
                f"CDS feature {fragment.start}-{fragment.end} runs past the end "
                f"of sequence {seqid!r} (length {len(chromosome)})",
            )

    # 5' -> 3' traversal of the fragments.
    traversal = ascending if strand == "+" else list(reversed(ascending))

    pieces: list[str] = []
    positions: list[int] = []
    for fragment in traversal:
        segment = chromosome[fragment.start - 1 : fragment.end].upper()
        fragment_positions = list(range(fragment.start, fragment.end + 1))
        if strand == "-":
            segment = _complement(segment)[::-1]
            fragment_positions.reverse()
        pieces.append(segment)
        positions.extend(fragment_positions)
    cds = "".join(pieces)

    # GFF3 phase p means (3 - p) mod 3 bases precede the first codon base
    # at the 5' end of the 5'-most fragment.
    leader = (3 - traversal[0].phase) % 3
    if leader >= len(cds) or (len(cds) - leader) % 3 != 0:
        raise _feature_error(
            source_name,
            traversal[0].line_number,
            f"CDS fragments of transcript {name!r} do not form an in-frame "
            f"coding sequence (length {len(cds)}, phase {traversal[0].phase})",
        )

    # Every later fragment's phase must agree with the cumulative frame.
    consumed = (traversal[0].end - traversal[0].start + 1) - leader
    for fragment in traversal[1:]:
        # GFF3 phase: p bases before the next codon boundary are encoded
        # as (3 - p) mod 3, so p == (-consumed) mod 3.
        expected_phase = (-consumed) % 3
        if fragment.phase != expected_phase:
            raise _feature_error(
                source_name,
                fragment.line_number,
                f"CDS fragment phase {fragment.phase} of transcript {name!r} is "
                f"inconsistent with the reading frame (expected {expected_phase})",
            )
        consumed += fragment.end - fragment.start + 1

    return _Transcript(name, seqid, strand, cds, leader, tuple(positions))


def _classify(codon: str, variant_codon: str, codon_index: int) -> str:
    amino = _CODON_TABLE[codon]
    variant_amino = _CODON_TABLE[variant_codon]
    if amino == variant_amino:
        return "SYNONYMOUS"
    if amino == "*":
        return "STOP_LOST"
    if variant_amino == "*":
        return "STOP_GAINED"
    if amino == "M" and codon_index == 1:
        return "START_LOST"
    return "MISSENSE"


def _fallback_item(allele: str, consequence: str) -> str:
    return "|".join((allele, consequence, _IMPACT[consequence], "", "", "", ""))


def _annotate_alt(
    allele: str,
    ref: str,
    pos: int,
    transcripts: Sequence[_Transcript],
) -> list[str]:
    single_substitution = (
        _is_plain_allele(allele)
        and len(ref) == 1
        and len(allele) == 1
        and ref in _DNA_BASES
        and allele in _DNA_BASES
        and allele != ref
    )

    hits: list[tuple[_Transcript, int]] = []
    if single_substitution:
        for transcript in transcripts:  # Parent lexicographic order
            for cds_position, genome_position in enumerate(
                transcript.coding_positions, start=1
            ):
                if genome_position == pos:
                    hits.append((transcript, cds_position))
                    break

    if not hits:
        return [_fallback_item(allele, "NON_CODING" if single_substitution else "UNSUPPORTED")]

    decoded_hits: list[tuple[_Transcript, int, str, str]] = []
    for transcript, cds_position in hits:
        codon_index = (cds_position - 1) // 3
        within = (cds_position - 1) % 3
        codon = transcript.cds[
            transcript.leader
            + codon_index * 3 : transcript.leader
            + codon_index * 3
            + 3
        ]
        variant_bases = list(codon)
        variant_base = allele if transcript.strand == "+" else _complement(allele)
        variant_bases[within] = variant_base
        decoded_hits.append(
            (transcript, cds_position, codon, "".join(variant_bases))
        )

    # An ambiguous base inside any reference codon cannot be translated
    # with the standard code; report the allele once as unsupported
    # rather than guessing a consequence per transcript.
    if any(
        codon not in _CODON_TABLE or variant_codon not in _CODON_TABLE
        for _, _, codon, variant_codon in decoded_hits
    ):
        return [_fallback_item(allele, "UNSUPPORTED")]

    items: list[str] = []
    for transcript, cds_position, codon, variant_codon in decoded_hits:
        codon_index = (cds_position - 1) // 3
        consequence = _classify(codon, variant_codon, codon_index + 1)
        amino = _CODON_TABLE[codon]
        variant_amino = _CODON_TABLE[variant_codon]
        aa_change = f"{amino}{codon_index + 1}{variant_amino}"
        items.append(
            "|".join(
                (
                    allele,
                    consequence,
                    _IMPACT[consequence],
                    transcript.name,
                    str(cds_position),
                    f"{codon}>{variant_codon}",
                    aa_change,
                )
            )
        )
    return items


def _info_already_used(info: str) -> bool:
    if info in ("", "."):
        return False
    return any(item.split("=", 1)[0] == _INFO_KEY for item in info.split(";"))


def _meta_declares_gvann(meta_lines: Sequence[str]) -> int | None:
    """Return the 1-based source line of a conflicting meta line, if any."""
    for index, line in enumerate(meta_lines, start=1):
        if not line.startswith("##") or "=" not in line:
            continue
        key, _, value = line[2:].partition("=")
        if key == _INFO_KEY:
            return index
        if key == "INFO":
            head = value.strip()
            if head.startswith("<") and "ID=" in head:
                ident = head.split("ID=", 1)[1]
                ident = ident.split(",", 1)[0].split(">", 1)[0]
                if ident == _INFO_KEY:
                    return index
    return None


def _append_info(info: str, annotation: str) -> str:
    field = f"{_INFO_KEY}={annotation}"
    return field if info in ("", ".") else f"{info};{field}"


def annotate_vcf(
    document: VcfFile,
    reference: "Mapping[str, str] | Iterable[SequenceRecord]",
    features: "str | os.PathLike[str] | io.TextIOBase",
) -> VcfFile:
    """Annotate every record of *document* from GFF3 CDS features.

    *reference* is a mapping of CHROM to upper-case reference sequence or
    an iterable of :class:`~genome_variant.sequence_io.SequenceRecord`.
    *features* is a GFF3 text path or text stream; only ``CDS`` features
    are read.

    Each record keeps its fixed columns, ALT order, sample columns and
    existing INFO fields; one ``GVANN`` item is appended to the INFO
    column, with entries ordered by ALT allele and then by transcript
    Parent name.  A unique ``##INFO`` declaration for ``GVANN`` is added
    last to the meta-information.

    Raises :class:`AnnotationFormatError` for malformed GFF3 (bad
    columns, coordinates, phase, strand or Parent, fragments spanning
    sequences, overlapping fragments or an inconsistent reading frame)
    and when the VCF already declares or uses ``GVANN``.  Absent or
    duplicated reference sequences and REF disagreements raise
    :class:`~genome_variant.vcf.ReferenceMismatchError`, and out-of-bounds
    records raise :class:`~genome_variant.vcf.VcfFormatError`.
    """
    if isinstance(reference, Mapping):
        sequences = dict(reference)
        duplicates: set[str] = set()
    else:
        sequences, duplicates = _build_reference(reference)

    conflicting = _meta_declares_gvann(document.header.meta_lines)
    if conflicting is not None:
        prefix = f"{document.source}:" if document.source else ""
        raise AnnotationFormatError(
            f"{prefix}{conflicting}: meta-information already declares {_INFO_KEY}"
        )

    if isinstance(features, io.IOBase):
        stream: io.TextIOBase = features
        features_source = getattr(stream, "name", None)
        if not isinstance(features_source, str) or not features_source:
            features_source = "<stream>"
        close = False
    elif isinstance(features, (str, os.PathLike)):
        path = os.fspath(features)
        stream = open(path, "r", encoding="utf-8", newline="")
        features_source = path
        close = True
    else:
        raise TypeError("features must be a GFF3 text path or a text stream")

    try:
        fragments = _parse_gff(stream, features_source)
    finally:
        if close:
            stream.close()

    transcripts = tuple(
        _build_transcript(name, fragments[name], sequences, duplicates, features_source)
        for name in sorted(fragments)
    )
    by_sequence: dict[str, list[_Transcript]] = defaultdict(list)
    for transcript in transcripts:
        by_sequence[transcript.seqid].append(transcript)

    annotated_records: list[VcfRecord] = []
    for record in document.records:
        if _info_already_used(record.info):
            raise _record_annotation_error(
                document.source, record, f"record already carries a {_INFO_KEY} INFO item"
            )

        if record.chrom in duplicates:
            raise ReferenceMismatchError(
                f"CHROM {record.chrom} POS {record.pos}: CHROM occurs more than once "
                "in the reference"
            )
        try:
            chromosome = sequences[record.chrom]
        except KeyError:
            raise ReferenceMismatchError(
                f"CHROM {record.chrom} POS {record.pos}: no matching record in the reference"
            ) from None

        end = record.pos - 1 + len(record.ref)
        if record.pos < 1 or end > len(chromosome):
            raise _record_error(
                document.source,
                record,
                f"CHROM {record.chrom} POS {record.pos}: REF of length "
                f"{len(record.ref)} is out of bounds for a reference of length "
                f"{len(chromosome)}",
            )
        expected = chromosome[record.pos - 1 : end].upper()
        if expected != record.ref:
            raise ReferenceMismatchError(
                f"CHROM {record.chrom} POS {record.pos}: REF does not match the "
                f"reference: expected {expected!r}, found {record.ref!r}"
            )

        entries: list[str] = []
        for allele in record.alt:
            entries.extend(
                _annotate_alt(allele, record.ref, record.pos, by_sequence.get(record.chrom, ()))
            )

        annotated_records.append(
            replace(record, info=_append_info(record.info, ",".join(entries)))
        )

    header = VcfHeader(document.header.meta_lines + (_META_LINE,), document.header.samples)
    return VcfFile(header, tuple(annotated_records), document.source)
