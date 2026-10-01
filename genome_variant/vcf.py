"""Reference-aware VCF reading, normalization and writing.

Public API:

- :class:`VcfHeader` — the ``##`` meta-information lines (kept verbatim
  and in order) together with the sample names declared by the single
  ``#CHROM`` column header.
- :class:`VcfRecord` — one VCF data row: the eight fixed columns, an
  optional ``FORMAT`` column and the per-sample text columns.
- :class:`VcfFile` — a parsed VCF document pairing a header with its
  records.
- :func:`read_vcf` — parse VCF text from a path or text stream.
- :func:`normalize_record` / :func:`normalize_vcf` — trim common
  allele suffixes/prefixes and left-align indels against a reference.
- :func:`render_vcf` / :func:`write_vcf` — serialize a VCF document
  with tabs and ``"\\n"``.

Ordinary sequence alleles only contain the IUPAC DNA symbols
``ACGTRYSWKMBDHVN`` (lower-case input is upper-cased).  Records carrying
symbolic (``<...>``), break-end (brackets or a ``.`` join), spanning
(``*``) or missing (``.``) ALT alleles are validated against the
reference but never trimmed or shifted.
"""

from __future__ import annotations

import io
import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace

from .sequence_io import SequenceRecord, _VALID_BASES

__all__ = [
    "VcfHeader",
    "VcfRecord",
    "VcfFile",
    "VcfFormatError",
    "ReferenceMismatchError",
    "read_vcf",
    "normalize_record",
    "normalize_vcf",
    "render_vcf",
    "write_vcf",
]

_FIXED_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
_MIN_COLUMNS = 8


class VcfFormatError(ValueError):
    """Malformed VCF input (missing/duplicate header, bad columns or fields)."""


class ReferenceMismatchError(ValueError):
    """A record's CHROM is absent/duplicated in the reference or REF disagrees."""


@dataclass(frozen=True)
class VcfHeader:
    """The VCF header.

    ``meta_lines`` holds the ``##`` meta-information lines verbatim,
    without their line terminators and in input order.  ``samples``
    holds the sample names taken from the single ``#CHROM`` header line
    (an empty tuple for sample-less VCF).
    """

    meta_lines: tuple[str, ...] = ()
    samples: tuple[str, ...] = ()

    @property
    def column_names(self) -> tuple[str, ...]:
        """All column names as declared by the ``#CHROM`` line."""
        if self.samples:
            return _FIXED_COLUMNS + ("FORMAT",) + self.samples
        return _FIXED_COLUMNS


@dataclass(frozen=True)
class VcfRecord:
    """One VCF data row.

    The fixed columns are ``chrom``, ``pos`` (1-based), ``id``, ``ref``,
    ``alt`` (a tuple of ALT allele strings), ``qual``, ``filter`` and
    ``info``.  For a VCF with samples, ``format_text`` holds the FORMAT
    string and ``sample_text`` one verbatim column per sample; both are
    empty for sample-less VCF.  ``line_number`` is the 1-based source
    line of a parsed record (``0`` for records built directly).
    """

    chrom: str
    pos: int
    id: str
    ref: str
    alt: tuple[str, ...]
    qual: str
    filter: str
    info: str
    format_text: str | None = None
    sample_text: tuple[str, ...] = ()
    line_number: int = 0

    @property
    def alt_text(self) -> str:
        """The ALT column as it appears on disk."""
        return ",".join(self.alt)


@dataclass(frozen=True)
class VcfFile:
    """A parsed VCF document: a header and its records in input order.

    ``source`` is the path or stream label the document was read from
    and is used to qualify error messages during normalization.
    """

    header: VcfHeader
    records: tuple[VcfRecord, ...] = ()
    source: str = ""


def _line_error(source: str, line_number: int, message: str) -> VcfFormatError:
    return VcfFormatError(f"{source}:{line_number}: {message}")


def _is_plain_allele(allele: str) -> bool:
    """Whether *allele* is an ordinary sequence allele of real bases.

    The gVCF spanning value ``*``, the missing value ``.``, symbolic
    alleles (``<DEL>``) and break-end alleles (``]chr:123]N``, ``N.``)
    are not plain.
    """
    if not allele or allele == "*" or allele == ".":
        return False
    if "<" in allele or ">" in allele:
        return False
    if "[" in allele or "]" in allele or "." in allele or "*" in allele:
        return False
    return True


def _validate_allele(allele: str) -> str:
    """Upper-case and validate one ordinary sequence allele."""
    upper = allele.upper()
    for offset, symbol in enumerate(upper):
        if symbol not in _VALID_BASES:
            raise ValueError(
                f"invalid allele symbol {symbol!r} at offset {offset} of {allele!r}"
            )
    return upper


def _parse_pos(pos_text: str, source: str, number: int) -> int:
    if not pos_text or any(char < "0" or char > "9" for char in pos_text):
        raise _line_error(
            source, number, f"POS {pos_text!r} is not a positive integer"
        )
    pos = int(pos_text)
    if pos <= 0:
        raise _line_error(
            source, number, f"POS {pos_text!r} is not a positive integer"
        )
    return pos


def _parse_record_line(
    line: str, number: int, source: str, samples: Sequence[str]
) -> VcfRecord:
    fields = line.split("\t")
    if len(fields) < _MIN_COLUMNS:
        raise _line_error(
            source,
            number,
            f"expected at least 8 tab-separated columns, got {len(fields)}",
        )
    if samples:
        expected = _MIN_COLUMNS + 1 + len(samples)
        if len(fields) != expected:
            raise _line_error(
                source,
                number,
                f"expected {expected} columns to match the #CHROM header, got {len(fields)}",
            )
        format_text = fields[8]
        if not format_text:
            raise _line_error(source, number, "empty FORMAT column")
        sample_text = tuple(fields[9:])
    else:
        if len(fields) != _MIN_COLUMNS:
            raise _line_error(
                source,
                number,
                f"expected 8 columns to match the sample-less #CHROM header, "
                f"got {len(fields)}",
            )
        format_text = None
        sample_text = ()

    chrom, pos_text, ident, ref_text, alt_text, qual, flt, info = fields[:8]

    if not chrom:
        raise _line_error(source, number, "empty CHROM column")
    pos = _parse_pos(pos_text, source, number)

    if not ref_text:
        raise _line_error(source, number, "empty REF allele")
    try:
        ref = _validate_allele(ref_text)
    except ValueError as exc:
        raise _line_error(source, number, str(exc)) from None

    if not alt_text:
        raise _line_error(source, number, "empty ALT column")
    alts: list[str] = []
    for allele in alt_text.split(","):
        if allele == "":
            raise _line_error(source, number, "empty ALT allele")
        if not _is_plain_allele(allele):
            alts.append(allele)
            continue
        try:
            alts.append(_validate_allele(allele))
        except ValueError as exc:
            raise _line_error(source, number, str(exc)) from None

    return VcfRecord(
        chrom=chrom,
        pos=pos,
        id=ident,
        ref=ref,
        alt=tuple(alts),
        qual=qual,
        filter=flt,
        info=info,
        format_text=format_text,
        sample_text=sample_text,
        line_number=number,
    )


def _parse_header_line(
    line: str, number: int, source: str
) -> tuple[str, ...]:
    """Validate the ``#CHROM`` line and return its sample names."""
    columns = line.split("\t")
    if tuple(columns[:8]) != _FIXED_COLUMNS:
        raise _line_error(
            source,
            number,
            "the first eight #CHROM header columns must be exactly "
            "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO",
        )
    tail = columns[8:]
    if not tail:
        return ()
    if tail[0] != "FORMAT":
        raise _line_error(
            source, number, "the ninth column, when present, must be FORMAT"
        )
    samples = tuple(tail[1:])
    if not samples:
        raise _line_error(source, number, "FORMAT column requires at least one sample")
    if any(not sample for sample in samples):
        raise _line_error(source, number, "empty sample name in #CHROM header")
    return samples


def read_vcf(
    source: "str | os.PathLike[str] | io.TextIOBase",
) -> VcfFile:
    """Read and validate a whole VCF document from *source*.

    *source* is either a filesystem path (``str`` or :class:`os.PathLike`)
    or an already-open text stream.  ``##`` meta-information lines are
    kept verbatim and in order; exactly one legal ``#CHROM`` line is
    required and every data row must match its columns.  Records are
    returned in input order.

    Raises :class:`VcfFormatError` for malformed input and :class:`OSError`
    when opening a path fails.
    """
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
        meta_lines: list[str] = []
        samples: tuple[str, ...] | None = None
        records: list[VcfRecord] = []

        for number, raw_line in enumerate(stream, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith("##"):
                if samples is not None:
                    raise _line_error(
                        source_name,
                        number,
                        "meta-information line after the #CHROM header",
                    )
                meta_lines.append(line)
            elif line.startswith("#"):
                if samples is not None:
                    raise _line_error(
                        source_name, number, "duplicate #CHROM column header"
                    )
                samples = _parse_header_line(line, number, source_name)
            else:
                if samples is None:
                    raise _line_error(
                        source_name, number, "data row before the #CHROM column header"
                    )
                records.append(
                    _parse_record_line(line, number, source_name, samples)
                )

        if samples is None:
            raise VcfFormatError(
                f"{source_name}:1: missing #CHROM column header"
            )
        return VcfFile(
            VcfHeader(tuple(meta_lines), samples), tuple(records), source_name
        )
    finally:
        if close:
            stream.close()


def _build_reference(
    reference: Iterable[SequenceRecord],
) -> tuple[dict[str, str], set[str]]:
    """Map CHROM to sequence, tracking identifiers appearing more than once."""
    sequences: dict[str, str] = {}
    duplicates: set[str] = set()
    for record in reference:
        if record.identifier in sequences:
            duplicates.add(record.identifier)
        else:
            sequences[record.identifier] = record.sequence
    return sequences, duplicates


def _record_error(source: str, record: VcfRecord, message: str) -> VcfFormatError:
    if source and record.line_number:
        return VcfFormatError(f"{source}:{record.line_number}: {message}")
    if record.line_number:
        return VcfFormatError(f"line {record.line_number}: {message}")
    return VcfFormatError(message)


def _resolve_chromosome(
    record: VcfRecord,
    sequences: Mapping[str, str],
    duplicates: set[str],
) -> str:
    if record.chrom in duplicates:
        raise ReferenceMismatchError(
            f"CHROM {record.chrom} POS {record.pos}: CHROM occurs more than once "
            "in the reference"
        )
    try:
        return sequences[record.chrom]
    except KeyError:
        raise ReferenceMismatchError(
            f"CHROM {record.chrom} POS {record.pos}: no matching record in the reference"
        ) from None


def _minimize(
    pos: int, ref: str, alts: Sequence[str]
) -> tuple[int, str, tuple[str, ...]]:
    """Remove the common suffix and prefix shared by all alleles.

    At least one base is always retained per allele and *pos* (1-based)
    is advanced by the number of stripped leading bases.
    """
    # Common suffix: compare each ALT against REF up to the shorter
    # length, always leaving at least one base per allele.
    suffix = 0
    while True:
        index = len(ref) - 1 - suffix
        if index <= 0:
            break
        ref_base = ref[index]
        ok = True
        for alt in alts:
            alt_index = len(alt) - 1 - suffix
            if alt_index <= 0 or alt[alt_index] != ref_base:
                ok = False
                break
        if not ok:
            break
        suffix += 1
    if suffix:
        ref = ref[: len(ref) - suffix]
        alts = tuple(alt[: len(alt) - suffix] for alt in alts)

    # Common prefix: at least one base must remain in every allele.
    prefix = 0
    shortest = min(len(ref), min(len(alt) for alt in alts))
    while prefix < shortest - 1 and all(
        alt[prefix] == ref[prefix] for alt in alts
    ):
        prefix += 1
    if prefix:
        ref = ref[prefix:]
        alts = tuple(alt[prefix:] for alt in alts)
        pos += prefix

    return pos, ref, tuple(alts)


def _lengths_differ(ref: str, alts: Sequence[str]) -> bool:
    return any(len(alt) != len(ref) for alt in alts)


def _normalize_against(
    record: VcfRecord,
    sequences: Mapping[str, str],
    duplicates: set[str],
    source: str,
) -> VcfRecord:
    chromosome = _resolve_chromosome(record, sequences, duplicates)

    end = record.pos - 1 + len(record.ref)
    if record.pos < 1 or end > len(chromosome):
        raise _record_error(
            source,
            record,
            f"CHROM {record.chrom} POS {record.pos}: REF of length {len(record.ref)} "
            f"is out of bounds for a reference of length {len(chromosome)}",
        )

    expected = chromosome[record.pos - 1 : end].upper()
    if expected != record.ref:
        raise ReferenceMismatchError(
            f"CHROM {record.chrom} POS {record.pos}: REF does not match the "
            f"reference: expected {expected!r}, found {record.ref!r}"
        )

    if not all(_is_plain_allele(allele) for allele in record.alt):
        # Symbolic, break-end, spanning or missing ALT: verify only.
        return record

    pos, ref, alts = _minimize(record.pos, record.ref, record.alt)

    # Left-align: prepend the preceding reference base to every allele
    # and re-minimize.  A move only happens when re-minimization removes
    # a (different) right-hand base, i.e. the extended window represents
    # the same haplotypes; once a step leaves the alleles unchanged the
    # minimum POS has been reached.
    while pos > 1 and _lengths_differ(ref, alts):
        preceding = chromosome[pos - 2].upper()
        next_pos, next_ref, next_alts = _minimize(
            pos - 1,
            preceding + ref,
            tuple(preceding + alt for alt in alts),
        )
        if (next_pos, next_ref, next_alts) == (pos, ref, alts):
            break
        pos, ref, alts = next_pos, next_ref, next_alts

    return replace(record, pos=pos, ref=ref, alt=alts)


def normalize_record(
    record: VcfRecord,
    reference: "Mapping[str, str] | Iterable[SequenceRecord]",
) -> VcfRecord:
    """Normalize one record against *reference*.

    *reference* is either a mapping of CHROM to upper-case reference
    sequence or an iterable of :class:`~genome_variant.sequence_io.SequenceRecord`
    (duplicate identifiers raise :class:`ReferenceMismatchError`).

    The REF fragment at POS is always verified.  Records whose ALT
    contains a symbolic, break-end, spanning (``*``) or missing (``.``)
    allele are returned unchanged.  Otherwise all sequence alleles are
    trimmed of their shared suffix and prefix (one base is always kept)
    and, when allele lengths differ, left-shifted across repeated
    reference bases to the smallest POS representing the same
    haplotypes, re-minimizing after every shift.  ID, QUAL, FILTER,
    INFO, FORMAT, sample text and ALT order are preserved.

    Raises :class:`ReferenceMismatchError` for an absent/duplicate CHROM
    or a REF that disagrees with the reference, and
    :class:`VcfFormatError` for a record that runs past the reference.
    """
    if isinstance(reference, Mapping):
        sequences = dict(reference)
        duplicates: set[str] = set()
    else:
        sequences, duplicates = _build_reference(reference)
    return _normalize_against(record, sequences, duplicates, "")


def normalize_vcf(
    document: VcfFile,
    reference: "Mapping[str, str] | Iterable[SequenceRecord]",
) -> VcfFile:
    """Normalize every record of *document* against *reference*.

    Record order is preserved; see :func:`normalize_record` for the
    per-record rules.  Every record is validated before the normalized
    document is returned, so an error never comes with partial output.
    """
    if isinstance(reference, Mapping):
        sequences = dict(reference)
        duplicates: set[str] = set()
    else:
        sequences, duplicates = _build_reference(reference)

    normalized = tuple(
        _normalize_against(record, sequences, duplicates, document.source)
        for record in document.records
    )
    return VcfFile(document.header, normalized, document.source)


def _record_line(record: VcfRecord) -> str:
    fields = [
        record.chrom,
        str(record.pos),
        record.id,
        record.ref,
        record.alt_text,
        record.qual,
        record.filter,
        record.info,
    ]
    if record.format_text is not None:
        fields.append(record.format_text)
        fields.extend(record.sample_text)
    return "\t".join(fields)


def render_vcf(document: VcfFile) -> str:
    """Serialize *document* to VCF text.

    Output uses tabs and ``"\\n"`` line endings and always ends with
    exactly one trailing newline; a record-less VCF still writes the
    complete header (``##`` lines followed by the ``#CHROM`` line).
    """
    lines = list(document.header.meta_lines)
    lines.append("\t".join(document.header.column_names))
    lines.extend(_record_line(record) for record in document.records)
    return "\n".join(lines) + "\n"


def write_vcf(
    document: VcfFile,
    output: "str | os.PathLike[str] | io.TextIOBase",
) -> None:
    """Write *document* to *output* as VCF text.

    Tabs and ``"\\n"`` are used throughout with exactly one trailing
    newline; an empty-record VCF still writes the complete header.
    """
    text = render_vcf(document)
    if isinstance(output, io.IOBase):
        output.write(text)
    elif isinstance(output, (str, os.PathLike)):
        with open(os.fspath(output), "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
    else:
        raise TypeError("output must be a text path or a text stream")
