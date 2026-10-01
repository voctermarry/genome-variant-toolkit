"""Reference-aware VCF reading, normalization and writing.

Public API:

- :class:`VcfHeader` — the ``##`` meta-information lines (verbatim, in
  order) and the fields of the single ``#CHROM`` column header line.
- :class:`VcfRecord` — one VCF record: the eight fixed columns plus an
  optional FORMAT column and per-sample texts.
- :func:`read_vcf` — read a header and all records from VCF text (a path
  or a text stream), preserving input order.
- :func:`normalize_vcf` — validate records against a reference and
  normalize the sequence-type alleles (trim common affixes, left-shift
  length-changing alleles to the minimal POS).
- :func:`write_vcf` — write a header and records as tab-separated VCF
  text with ``"\\n"`` line endings.

Plain sequence REF/ALT alleles are upper-cased and limited to the IUPAC
DNA symbols ``ACGTRYSWKMBDHVN``.  Symbolic alleles (``<DEL>`` …),
breakend alleles (``]chr:pos]`` …), ``*`` and ``.`` are kept verbatim;
records carrying them are validated against the reference but never
trimmed or left-shifted.

Errors:

- :class:`VcfFormatError` — malformed VCF structure or alleles, and
  records that fall outside the reference sequence; messages carry the
  source name and a 1-based line number.
- :class:`ReferenceMismatchError` — a CHROM missing from or duplicated
  in the reference, or a REF allele inconsistent with the reference;
  messages carry CHROM, POS and the determinable expected/actual values.
"""

from __future__ import annotations

import io
import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace

__all__ = [
    "VcfFormatError",
    "ReferenceMismatchError",
    "VcfHeader",
    "VcfRecord",
    "read_vcf",
    "normalize_vcf",
    "write_vcf",
]

_VALID_BASES = frozenset("ACGTRYSWKMBDHVN")
_FIXED_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


class VcfFormatError(ValueError):
    """Malformed VCF input (header, columns, POS, REF/ALT, out-of-range)."""


class ReferenceMismatchError(ValueError):
    """A record cannot be reconciled with the reference sequence.

    ``chrom`` and ``pos`` locate the offending record; ``expected`` and
    ``actual`` hold the reference segment and the record's REF allele
    when both are determinable (``None`` otherwise).
    """

    def __init__(
        self,
        message: str,
        *,
        chrom: str | None = None,
        pos: int | None = None,
        expected: str | None = None,
        actual: str | None = None,
    ) -> None:
        super().__init__(message)
        self.chrom = chrom
        self.pos = pos
        self.expected = expected
        self.actual = actual


@dataclass(frozen=True)
class VcfHeader:
    """The VCF header.

    ``meta_lines`` holds the ``##`` lines verbatim (including the ``##``
    prefix) in their original order.  ``columns`` holds the fields of
    the single ``#CHROM`` column header line: the eight fixed columns,
    optionally followed by ``FORMAT`` and one column per sample.
    """

    meta_lines: tuple[str, ...]
    columns: tuple[str, ...]

    @property
    def sample_names(self) -> tuple[str, ...]:
        """Sample column names (empty when the header has no FORMAT)."""
        if len(self.columns) > 9:
            return self.columns[9:]
        return ()


@dataclass(frozen=True)
class VcfRecord:
    """One VCF record.

    ``ref`` and sequence-type entries of ``alts`` are upper-case IUPAC
    sequences; symbolic, breakend, ``*`` and ``.`` alleles are kept
    verbatim.  ``qual``, ``filter`` and ``info`` are the raw column
    texts.  ``format`` is ``None`` when the header has no FORMAT column.
    ``line_number`` is the 1-based source line (``0`` when unknown) and
    is never written out.
    """

    chrom: str
    pos: int
    id: str = "."
    ref: str = ""
    alts: tuple[str, ...] = ()
    qual: str = "."
    filter: str = "."
    info: str = "."
    format: str | None = None
    samples: tuple[str, ...] = ()
    line_number: int = 0


def _line_error(source: str, line_number: int, message: str) -> VcfFormatError:
    return VcfFormatError(f"{source}:{line_number}: {message}")


def _is_sequence_allele(allele: str) -> bool:
    return bool(allele) and all(symbol in _VALID_BASES for symbol in allele)


def _parse_alt_allele(allele: str, source: str, line_number: int) -> str:
    if not allele:
        raise _line_error(source, line_number, "empty ALT allele")
    upper = allele.upper()
    if _is_sequence_allele(upper):
        return upper
    if allele in ("*", "."):
        return allele
    if len(allele) >= 2 and allele.startswith("<") and allele.endswith(">"):
        return allele
    if "[" in allele or "]" in allele:
        return allele
    raise _line_error(source, line_number, f"invalid ALT allele {allele!r}")


def _check_columns(columns: tuple[str, ...], source: str, line_number: int) -> None:
    if len(columns) < len(_FIXED_COLUMNS):
        raise _line_error(
            source,
            line_number,
            f"#CHROM column header has {len(columns)} fields; "
            f"expected at least {len(_FIXED_COLUMNS)}",
        )
    if columns[: len(_FIXED_COLUMNS)] != _FIXED_COLUMNS:
        raise _line_error(
            source,
            line_number,
            "#CHROM column header must start with the fixed columns "
            + " ".join(_FIXED_COLUMNS),
        )
    if len(columns) > len(_FIXED_COLUMNS) and columns[8] != "FORMAT":
        raise _line_error(
            source,
            line_number,
            f"#CHROM column header field 9 must be 'FORMAT', got {columns[8]!r}",
        )


def _parse_record(
    line: str, source: str, line_number: int, columns: tuple[str, ...]
) -> VcfRecord:
    fields = line.split("\t")
    if len(fields) != len(columns):
        raise _line_error(
            source,
            line_number,
            f"record has {len(fields)} columns but the #CHROM column header "
            f"declares {len(columns)}",
        )
    pos_text = fields[1]
    if not pos_text.isascii() or not pos_text.isdigit() or int(pos_text) < 1:
        raise _line_error(
            source,
            line_number,
            f"invalid POS {pos_text!r}: expected a positive integer",
        )
    pos = int(pos_text)

    ref = fields[3].upper()
    if not ref:
        raise _line_error(source, line_number, "empty REF allele")
    if not _is_sequence_allele(ref):
        raise _line_error(source, line_number, f"invalid REF allele {fields[3]!r}")

    if not fields[4]:
        raise _line_error(source, line_number, "empty ALT")
    alts = tuple(
        _parse_alt_allele(allele, source, line_number)
        for allele in fields[4].split(",")
    )

    has_format = len(columns) > len(_FIXED_COLUMNS)
    return VcfRecord(
        chrom=fields[0],
        pos=pos,
        id=fields[2],
        ref=ref,
        alts=alts,
        qual=fields[5],
        filter=fields[6],
        info=fields[7],
        format=fields[8] if has_format else None,
        samples=tuple(fields[9:]) if has_format else (),
        line_number=line_number,
    )


def _parse_vcf(stream: Iterable[str], source: str) -> tuple[VcfHeader, list[VcfRecord]]:
    meta_lines: list[str] = []
    columns: tuple[str, ...] | None = None
    records: list[VcfRecord] = []
    last_line = 0

    for number, raw_line in enumerate(stream, start=1):
        last_line = number
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("##"):
            if columns is not None:
                raise _line_error(
                    source,
                    number,
                    "meta-information line after the #CHROM column header",
                )
            meta_lines.append(line)
            continue
        if line.startswith("#"):
            if columns is not None:
                raise _line_error(source, number, "duplicate #CHROM column header")
            columns = tuple(line.split("\t"))
            _check_columns(columns, source, number)
            continue
        if columns is None:
            raise _line_error(
                source, number, "record before the #CHROM column header"
            )
        records.append(_parse_record(line, source, number, columns))

    if columns is None:
        raise _line_error(source, last_line + 1, "missing #CHROM column header")
    return VcfHeader(tuple(meta_lines), columns), records


def read_vcf(
    source: "str | os.PathLike[str] | io.TextIOBase",
) -> tuple[VcfHeader, list[VcfRecord]]:
    """Read VCF text from *source* (a text path or a text stream).

    Returns the :class:`VcfHeader` and all :class:`VcfRecord` objects in
    input order.  The ``##`` meta-information lines are preserved
    verbatim and in order; exactly one valid ``#CHROM`` column header is
    required and every record must match its column count.  Malformed
    input raises :class:`VcfFormatError` with the source name and a
    1-based line number; errors while opening a path propagate as
    :class:`OSError`.
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
        return _parse_vcf(stream, source_name)
    finally:
        if close:
            stream.close()


def _trim_alleles(alleles: list[str], pos: int) -> tuple[list[str], int]:
    """Remove the longest common suffix then prefix shared by *alleles*.

    Every allele always keeps at least one base; *pos* is advanced by
    the number of trimmed prefix bases.
    """
    keep = min(len(allele) for allele in alleles) - 1
    suffix = 0
    first = alleles[0]
    while suffix < keep and all(
        allele[-1 - suffix] == first[-1 - suffix] for allele in alleles[1:]
    ):
        suffix += 1
    if suffix:
        alleles = [allele[: len(allele) - suffix] for allele in alleles]

    keep = min(len(allele) for allele in alleles) - 1
    prefix = 0
    while prefix < keep and all(
        allele[prefix] == alleles[0][prefix] for allele in alleles[1:]
    ):
        prefix += 1
    if prefix:
        alleles = [allele[prefix:] for allele in alleles]
        pos += prefix
    return alleles, pos


def _left_shift(
    alleles: list[str], pos: int, reference: str
) -> tuple[list[str], int]:
    """Shift length-changing alleles to the minimal equivalent POS.

    A shift by one base is possible while every allele ends with the
    same base (the last base of the REF region, which equals the
    reference there); the preceding reference base is then prepended to
    every allele and the last base dropped, which leaves every
    haplotype unchanged.  Alleles are re-trimmed after each move.
    """
    while pos > 1:
        last_base = alleles[0][-1]
        if any(allele[-1] != last_base for allele in alleles[1:]):
            break
        previous = reference[pos - 2]
        alleles = [previous + allele[:-1] for allele in alleles]
        pos -= 1
        alleles, pos = _trim_alleles(alleles, pos)
    return alleles, pos


def _normalize_record(
    record: VcfRecord,
    sequences: Mapping[str, str],
    duplicates: Mapping[str, int],
    source: str,
) -> VcfRecord:
    chrom = record.chrom
    if chrom in duplicates:
        raise ReferenceMismatchError(
            f"{chrom}:{record.pos}: chromosome {chrom!r} occurs "
            f"{duplicates[chrom]} times in the reference",
            chrom=chrom,
            pos=record.pos,
        )
    reference = sequences.get(chrom)
    if reference is None:
        raise ReferenceMismatchError(
            f"{chrom}:{record.pos}: chromosome {chrom!r} is not present "
            "in the reference",
            chrom=chrom,
            pos=record.pos,
        )

    start = record.pos - 1
    end = start + len(record.ref)
    if end > len(reference):
        location = (
            f"{source}:{record.line_number}" if record.line_number else source
        )
        raise VcfFormatError(
            f"{location}: record at {chrom}:{record.pos} with REF of length "
            f"{len(record.ref)} extends beyond reference sequence "
            f"{chrom!r} of length {len(reference)}"
        )
    expected = reference[start:end]
    if expected != record.ref:
        raise ReferenceMismatchError(
            f"{chrom}:{record.pos}: REF {record.ref!r} does not match the "
            f"reference {expected!r}",
            chrom=chrom,
            pos=record.pos,
            expected=expected,
            actual=record.ref,
        )

    if not record.alts or any(
        not _is_sequence_allele(allele) for allele in record.alts
    ):
        # Symbolic, breakend, '*' and '.' alleles are never trimmed or
        # left-shifted; the REF validation above still applies.
        return record

    alleles, pos = _trim_alleles([record.ref, *record.alts], record.pos)
    if len({len(allele) for allele in alleles}) > 1:
        alleles, pos = _left_shift(alleles, pos, reference)
    if pos == record.pos and alleles[0] == record.ref and alleles[1:] == list(record.alts):
        return record
    return replace(record, pos=pos, ref=alleles[0], alts=tuple(alleles[1:]))


def normalize_vcf(
    records: Iterable[VcfRecord],
    reference: "Mapping[str, str] | Iterable[object]",
    source: str = "<stream>",
) -> list[VcfRecord]:
    """Validate and normalize *records* against *reference*.

    *reference* is either a mapping of chromosome name to sequence or an
    iterable of objects with ``identifier`` and ``sequence`` attributes
    (such as :class:`~genome_variant.sequence_io.SequenceRecord`).
    Records are processed in input order and the normalized records are
    returned as a list in the same order; ALT order and the ID, QUAL,
    FILTER, INFO, FORMAT and sample texts are preserved.

    Every record's REF is checked against the reference at POS.  Records
    whose alleles are all plain sequences are trimmed (longest common
    suffix, then prefix, always keeping one base per allele) and, when
    the alleles differ in length, left-shifted to the minimal POS
    representing the same haplotype.  Records carrying a symbolic,
    breakend, ``*`` or ``.`` ALT are validated but never modified.

    *source* names the VCF origin in out-of-range :class:`VcfFormatError`
    messages (together with each record's 1-based line number).
    """
    if isinstance(reference, Mapping):
        items = list(reference.items())
    else:
        items = [
            (getattr(entry, "identifier"), getattr(entry, "sequence"))
            for entry in reference
        ]
    sequences: dict[str, str] = {}
    duplicates: dict[str, int] = {}
    for name, sequence in items:
        if name in sequences:
            duplicates[name] = duplicates.get(name, 1) + 1
        else:
            sequences[name] = sequence

    return [
        _normalize_record(record, sequences, duplicates, source)
        for record in records
    ]


def _format_record(record: VcfRecord, header: VcfHeader) -> str:
    if (
        not isinstance(record.pos, int)
        or isinstance(record.pos, bool)
        or record.pos < 1
    ):
        raise ValueError(
            f"record at {record.chrom!r}: POS must be a positive integer"
        )
    if not record.ref:
        raise ValueError(f"record at {record.chrom!r}: empty REF allele")
    if not record.alts:
        raise ValueError(f"record at {record.chrom!r}: no ALT allele")
    fields = [
        record.chrom,
        str(record.pos),
        record.id,
        record.ref,
        ",".join(record.alts),
        record.qual,
        record.filter,
        record.info,
    ]
    if record.format is not None:
        fields.append(record.format)
        fields.extend(record.samples)
    if len(fields) != len(header.columns):
        raise ValueError(
            f"record at {record.chrom!r}:{record.pos} has {len(fields)} columns "
            f"but the header declares {len(header.columns)}"
        )
    return "\t".join(fields)


def write_vcf(
    header: VcfHeader,
    records: Iterable[VcfRecord],
    output: "str | os.PathLike[str] | io.TextIOBase",
) -> None:
    """Write *header* and *records* to *output* as VCF text.

    *output* is a filesystem path (``str`` or :class:`os.PathLike`) or an
    already-open text stream.  The ``##`` meta-information lines are
    written verbatim, followed by the ``#CHROM`` column header and one
    tab-separated line per record.  Output uses ``"\\n"`` line endings
    and ends with exactly one trailing newline, even when there are no
    records.
    """
    if (
        len(header.columns) < len(_FIXED_COLUMNS)
        or header.columns[: len(_FIXED_COLUMNS)] != _FIXED_COLUMNS
        or (len(header.columns) > len(_FIXED_COLUMNS) and header.columns[8] != "FORMAT")
    ):
        raise ValueError(
            "header columns must be the eight fixed VCF columns, optionally "
            "followed by FORMAT and one column per sample"
        )
    if isinstance(output, io.IOBase):
        stream: io.TextIOBase = output
        close = False
    elif isinstance(output, (str, os.PathLike)):
        stream = open(os.fspath(output), "w", encoding="utf-8", newline="")
        close = True
    else:
        raise TypeError("output must be a text path or a text stream")

    try:
        for meta_line in header.meta_lines:
            stream.write(meta_line + "\n")
        stream.write("\t".join(header.columns) + "\n")
        for record in records:
            stream.write(_format_record(record, header) + "\n")
    finally:
        if close:
            stream.close()
