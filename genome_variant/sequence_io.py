"""FASTA/FASTQ reading, validation and normalization.

Public surface:

* :class:`SequenceRecord` -- one sequence with identifier, description and
  optional Phred quality scores.
* :func:`read_sequences` -- lazily read records from a text path or stream.
* :func:`write_sequences` -- write records in FASTA or FASTQ form.
* :class:`SequenceFormatError`, :class:`SequenceValidationError`.

Sequences are normalized to upper case and may only contain the IUPAC DNA
symbols ``ACGTRYSWKMBDHVN``. FASTQ quality strings are decoded with the
Phred+33 convention (ASCII 33-126).
"""

from __future__ import annotations

import itertools
import os
from dataclasses import dataclass
from typing import IO, Iterator, Literal, Union

__all__ = [
    "SequenceRecord",
    "SequenceFormatError",
    "SequenceValidationError",
    "read_sequences",
    "write_sequences",
    "FASTA",
    "FASTQ",
    "AUTO",
]

FASTA = "fasta"
FASTQ = "fastq"
AUTO = "auto"

_FORMAT_ALIASES = {
    "fasta": FASTA,
    "fastq": FASTQ,
    "auto": AUTO,
}

#: IUPAC nucleotide symbols allowed after upper-casing.
ALLOWED_BASES = frozenset("ACGTRYSWKMBDHVN")

#: Phred+33 quality scores use ASCII 33 through 126.
QUAL_MIN = 33
QUAL_MAX = 126

PathOrStream = Union[str, os.PathLike[str], IO[str]]
Format = Literal["fasta", "fastq", "auto"]


class SequenceFormatError(ValueError):
    """Raised when input cannot be parsed as the declared sequence format."""


class SequenceValidationError(ValueError):
    """Raised when a sequence contains a symbol outside the allowed alphabet."""


@dataclass(frozen=True)
class SequenceRecord:
    """A single sequence record.

    ``quality`` is ``None`` for FASTA records and a tuple of integer Phred
    scores for FASTQ records.
    """

    id: str
    sequence: str
    description: str = ""
    quality: tuple[int, ...] | None = None


class _Builder:
    """Mutable accumulator for one record while it is being parsed."""

    __slots__ = (
        "identifier",
        "description",
        "raw_title",
        "header_line",
        "seq_parts",
        "qual_parts",
        "seq_length",
    )

    def __init__(self, line: str, lineno: int, marker: str, source: str) -> None:
        rest = line[1:]
        stripped = rest.lstrip()
        if not stripped:
            raise SequenceFormatError(
                f"{source}:{lineno}: empty record identifier in "
                f"{marker!r} header"
            )
        parts = stripped.split(maxsplit=1)
        self.identifier = parts[0]
        self.description = parts[1] if len(parts) == 2 else ""
        # Raw title text (everything after the marker, byte-for-byte), used
        # to validate a repeated FASTQ "+" title.
        self.raw_title = rest
        self.header_line = lineno
        self.seq_parts: list[str] = []
        self.qual_parts: list[str] = []
        self.seq_length = 0

    def append_sequence(self, chunk: str, lineno: int, source: str) -> None:
        offset = self.seq_length
        for index, symbol in enumerate(chunk):
            if symbol not in ALLOWED_BASES:
                raise SequenceValidationError(
                    f"{source}:{lineno}: invalid symbol {symbol!r} in record "
                    f"{self.identifier!r} at position {offset + index + 1}"
                )
        self.seq_parts.append(chunk)
        self.seq_length += len(chunk)

    def sequence(self) -> str:
        return "".join(self.seq_parts)

    def to_record(self) -> SequenceRecord:
        return SequenceRecord(
            id=self.identifier,
            sequence="".join(self.seq_parts),
            description=self.description,
            quality=None,
        )


def _clean_line(raw: str) -> str:
    """Remove a single trailing newline, leaving every other byte intact."""
    return raw.rstrip("\r\n")


def _is_blank(line: str) -> bool:
    return line.strip() == ""


def _read_fasta(lines: Iterator[tuple[int, str]], source: str) -> Iterator[SequenceRecord]:
    builder: _Builder | None = None
    for lineno, line in lines:
        if builder is None:
            if _is_blank(line):
                continue
            if not line.startswith(">"):
                raise SequenceFormatError(
                    f"{source}:{lineno}: expected FASTA header starting with '>'"
                )
            builder = _Builder(line, lineno, ">", source)
        else:
            if _is_blank(line):
                continue
            if line.startswith(">"):
                yield _finish_fasta(builder, source)
                builder = _Builder(line, lineno, ">", source)
            else:
                # Drop inline whitespace and normalize to upper case.
                builder.append_sequence("".join(line.split()).upper(), lineno, source)
    if builder is None:
        raise SequenceFormatError(
            f"{source}:1: empty input: no FASTA records found"
        )
    yield _finish_fasta(builder, source)


def _finish_fasta(builder: _Builder, source: str) -> SequenceRecord:
    if builder.seq_length == 0:
        raise SequenceFormatError(
            f"{source}:{builder.header_line}: record {builder.identifier!r} "
            f"has an empty sequence"
        )
    return builder.to_record()


def _read_fastq(lines: Iterator[tuple[int, str]], source: str) -> Iterator[SequenceRecord]:
    # State of a small machine: header -> seq -> qual -> header ...
    builder: _Builder | None = None
    expecting = "header"
    qual_length = 0
    saw_record = False

    for lineno, line in lines:
        if expecting == "header":
            if _is_blank(line):
                continue
            if not line.startswith("@"):
                raise SequenceFormatError(
                    f"{source}:{lineno}: expected FASTQ header starting with '@'"
                )
            builder = _Builder(line, lineno, "@", source)
            expecting = "seq"
        elif expecting == "seq":
            if _is_blank(line):
                continue
            # '+' is not a valid base, so a line starting with it is
            # unambiguously the separator even for multi-line sequences.
            if line.startswith("+"):
                repeated = line[1:]
                if repeated and repeated != builder.raw_title:
                    raise SequenceFormatError(
                        f"{source}:{lineno}: '+' separator title does not match "
                        f"the header of record {builder.identifier!r}"
                    )
                expecting = "qual"
                qual_length = 0
                if builder.seq_length == 0:
                    raise SequenceFormatError(
                        f"{source}:{builder.header_line}: record "
                        f"{builder.identifier!r} has an empty sequence"
                    )
            else:
                builder.append_sequence("".join(line.split()).upper(), lineno, source)
        else:  # expecting quality characters
            if _is_blank(line):
                # A blank quality line adds nothing; whitespace-only lines
                # carry spaces (ASCII 32) and fail the range check below.
                continue
            for symbol in line:
                if not QUAL_MIN <= ord(symbol) <= QUAL_MAX:
                    raise SequenceFormatError(
                        f"{source}:{lineno}: invalid quality character "
                        f"{symbol!r} (outside ASCII {QUAL_MIN}-{QUAL_MAX}) "
                        f"in record {builder.identifier!r}"
                    )
            builder.qual_parts.append(line)
            qual_length += len(line)
            if qual_length > builder.seq_length:
                raise SequenceFormatError(
                    f"{source}:{lineno}: quality length {qual_length} exceeds "
                    f"sequence length {builder.seq_length} in record "
                    f"{builder.identifier!r}"
                )
            if qual_length == builder.seq_length:
                yield _finish_fastq(builder, qual_length, source)
                builder = None
                expecting = "header"
                saw_record = True

    if builder is not None:
        if expecting == "seq":
            raise SequenceFormatError(
                f"{source}:{builder.header_line}: incomplete FASTQ record "
                f"{builder.identifier!r}: missing '+' separator"
            )
        raise SequenceFormatError(
            f"{source}:{builder.header_line}: incomplete FASTQ record "
            f"{builder.identifier!r}: quality length {qual_length} does not "
            f"match sequence length {builder.seq_length}"
        )
    if not saw_record:
        raise SequenceFormatError(
            f"{source}:1: empty input: no FASTQ records found"
        )


def _finish_fastq(builder: _Builder, qual_length: int, source: str) -> SequenceRecord:
    if builder.seq_length == 0:
        raise SequenceFormatError(
            f"{source}:{builder.header_line}: record {builder.identifier!r} "
            f"has an empty sequence"
        )
    if qual_length != builder.seq_length:  # pragma: no cover - guarded at EOF
        raise SequenceFormatError(
            f"{source}:{builder.header_line}: quality length {qual_length} does "
            f"not match sequence length {builder.seq_length} in record "
            f"{builder.identifier!r}"
        )
    quality = tuple(ord(symbol) - QUAL_MIN for symbol in "".join(builder.qual_parts))
    record = builder.to_record()
    return SequenceRecord(record.id, record.sequence, record.description, quality)


def _iter_sequences(
    stream: IO[str], fmt: str, source: str, close_stream: bool
) -> Iterator[SequenceRecord]:
    try:
        numbered: Iterator[tuple[int, str]] = (
            (lineno, _clean_line(raw)) for lineno, raw in enumerate(stream, 1)
        )
        if fmt == AUTO:
            first: tuple[int, str] | None = None
            for candidate in numbered:
                if not _is_blank(candidate[1]):
                    first = candidate
                    break
            if first is None:
                raise SequenceFormatError(
                    f"{source}: empty input: no sequence records found"
                )
            lineno, line = first
            if line.startswith(">"):
                fmt = FASTA
            elif line.startswith("@"):
                fmt = FASTQ
            else:
                raise SequenceFormatError(
                    f"{source}:{lineno}: cannot determine sequence format: "
                    f"first non-empty line must start with '>' or '@'"
                )
            numbered = itertools.chain([first], numbered)

        if fmt == FASTA:
            yield from _read_fasta(numbered, source)
        elif fmt == FASTQ:
            yield from _read_fastq(numbered, source)
        else:  # pragma: no cover - validated by the public entry point
            raise ValueError(f"unknown format: {fmt!r}")
    finally:
        if close_stream:
            stream.close()


def read_sequences(
    source: PathOrStream,
    format: Format = AUTO,
    *,
    source_name: str | None = None,
) -> Iterator[SequenceRecord]:
    """Lazily yield :class:`SequenceRecord` objects from ``source``.

    ``source`` is either a text file path (opened for reading; an ``OSError``
    from opening propagates to the caller) or an already open text stream.
    ``format`` is ``"fasta"``, ``"fastq"`` or ``"auto"`` (the default), which
    decides from the first non-empty line. Parsing and I/O errors raised while
    iterating are never swallowed.
    """
    fmt = _FORMAT_ALIASES.get(format)
    if fmt is None:
        raise ValueError(
            f"unknown format {format!r}: expected 'fasta', 'fastq' or 'auto'"
        )

    if isinstance(source, (str, os.PathLike)):
        path = os.fspath(source)
        # Open eagerly so a failure surfaces from the call itself, rather
        # than from the first ``next()``.
        stream = open(path, "r", encoding="utf-8", newline="")
        name = path if source_name is None else source_name
        return _iter_sequences(stream, fmt, name, True)

    name = source_name
    if name is None:
        attr = getattr(source, "name", None)
        name = attr if isinstance(attr, str) and attr else "<input>"
    return _iter_sequences(source, fmt, name, False)


def _record_title(record: SequenceRecord, marker: str) -> str:
    if not record.id:
        raise ValueError("record identifier must not be empty")
    if record.description:
        return f"{marker}{record.id} {record.description}"
    return f"{marker}{record.id}"


def write_sequences(
    records: Iterator[SequenceRecord] | list[SequenceRecord],
    output: PathOrStream,
    format: Format = FASTA,
    *,
    line_width: int = 80,
) -> None:
    """Write ``records`` to ``output`` in FASTA or FASTQ form.

    FASTA sequences are wrapped at ``line_width`` columns (a positive
    integer); the option does not affect FASTQ, which always uses the fixed
    four-line layout. Writing a record without quality scores as FASTQ raises
    :class:`ValueError`. Lines end with ``"\\n"`` and the output ends with
    exactly one trailing newline whenever it contains records.
    """
    if format not in (FASTA, FASTQ):
        raise ValueError(
            f"unknown output format {format!r}: expected 'fasta' or 'fastq'"
        )
    if isinstance(line_width, bool) or not isinstance(line_width, int):
        raise ValueError("line_width must be a positive integer")
    if line_width <= 0:
        raise ValueError("line_width must be a positive integer")

    close_stream = False
    if isinstance(output, (str, os.PathLike)):
        stream = open(os.fspath(output), "w", encoding="utf-8", newline="\n")
        close_stream = True
    else:
        stream = output

    try:
        for record in records:
            if format == FASTA:
                stream.write(_record_title(record, ">") + "\n")
                sequence = record.sequence
                if not sequence:
                    raise ValueError(
                        f"record {record.id!r} has an empty sequence"
                    )
                for start in range(0, len(sequence), line_width):
                    stream.write(sequence[start : start + line_width] + "\n")
            else:
                if record.quality is None:
                    raise ValueError(
                        f"cannot write record {record.id!r} as FASTQ: no "
                        f"quality scores available"
                    )
                if len(record.quality) != len(record.sequence):
                    raise ValueError(
                        f"cannot write record {record.id!r} as FASTQ: quality "
                        f"length {len(record.quality)} does not match sequence "
                        f"length {len(record.sequence)}"
                    )
                encoded = []
                for score in record.quality:
                    if not 0 <= score <= QUAL_MAX - QUAL_MIN:
                        raise ValueError(
                            f"quality score {score} in record {record.id!r} is "
                            f"outside the Phred+33 range 0-{QUAL_MAX - QUAL_MIN}"
                        )
                    encoded.append(chr(score + QUAL_MIN))
                stream.write(_record_title(record, "@") + "\n")
                stream.write(record.sequence + "\n")
                stream.write("+\n")
                stream.write("".join(encoded) + "\n")
    finally:
        if close_stream:
            stream.close()
