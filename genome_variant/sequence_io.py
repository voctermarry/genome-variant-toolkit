"""FASTA/FASTQ reading and normalization.

Public API:

- :class:`SequenceRecord` — one sequence with identifier, description and
  optional per-base quality values.
- :func:`read_sequences` — lazily read records from a text path or text
  stream; the format is ``"fasta"``, ``"fastq"`` or ``"auto"``.
- :func:`write_sequences` — write records in FASTA or FASTQ form.

Sequences only contain the IUPAC DNA symbols ``ACGTRYSWKMBDHVN``.  FASTQ
quality values use the Phred+33 encoding with ASCII characters 33–126.
"""

from __future__ import annotations

import io
import os
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

__all__ = [
    "SequenceRecord",
    "SequenceFormatError",
    "SequenceValidationError",
    "read_sequences",
    "write_sequences",
]

_VALID_BASES = frozenset("ACGTRYSWKMBDHVN")
_INPUT_FORMATS = ("fasta", "fastq", "auto")
_OUTPUT_FORMATS = ("fasta", "fastq")


class SequenceFormatError(ValueError):
    """Malformed FASTA/FASTQ input (structure, identifiers, lengths)."""


class SequenceValidationError(ValueError):
    """A sequence contains a symbol outside the allowed IUPAC set."""


@dataclass(frozen=True)
class SequenceRecord:
    """A single sequence record.

    ``identifier`` is the first whitespace-delimited field of the title
    line; ``description`` holds the remainder (``""`` when absent).
    ``quality`` is ``None`` for FASTA records and a tuple of Phred scores
    for FASTQ records.
    """

    identifier: str
    sequence: str
    description: str = ""
    quality: tuple[int, ...] | None = None


def _clean_sequence(raw: str) -> str:
    """Strip in-line whitespace and upper-case a raw sequence fragment."""
    return "".join(raw.split()).upper()


def _validate_sequence(identifier: str, sequence: str) -> None:
    for position, symbol in enumerate(sequence, start=1):
        if symbol not in _VALID_BASES:
            raise SequenceValidationError(
                f"record {identifier!r}: invalid symbol {symbol!r} at position {position}"
            )


def _split_title(title: str, source: str, line_number: int) -> tuple[str, str]:
    fields = title.split(None, 1)
    if not fields:
        raise SequenceFormatError(f"{source}:{line_number}: empty sequence identifier")
    identifier = fields[0]
    description = fields[1] if len(fields) == 2 else ""
    return identifier, description


def _line_error(source: str, line_number: int, message: str) -> SequenceFormatError:
    return SequenceFormatError(f"{source}:{line_number}: {message}")


def _fasta_records(
    lines: Iterator[tuple[int, str]], source: str
) -> Iterator[SequenceRecord]:
    header_line: str | None = None
    header_number = 0
    seq_parts: list[str] = []

    for number, raw_line in lines:
        line = raw_line.rstrip("\r\n")
        if line.startswith(">"):
            if header_line is not None:
                yield _finish_fasta(header_line, header_number, seq_parts, source)
            header_line = line
            header_number = number
            seq_parts = []
        elif line.strip():
            if header_line is None:
                raise _line_error(
                    source,
                    number,
                    "content before first title line (expected a line starting with '>')",
                )
            seq_parts.append(line)
        # Blank lines outside or inside a record are ignored.

    if header_line is not None:
        yield _finish_fasta(header_line, header_number, seq_parts, source)


def _finish_fasta(
    header_line: str, header_number: int, parts: list[str], source: str
) -> SequenceRecord:
    identifier, description = _split_title(header_line[1:], source, header_number)
    sequence = _clean_sequence("".join(parts))
    if not sequence:
        raise _line_error(source, header_number, "empty sequence")
    _validate_sequence(identifier, sequence)
    return SequenceRecord(identifier, sequence, description)


def _check_quality_segment(segment: str, source: str, line_number: int) -> None:
    for char in segment:
        value = ord(char)
        if value < 33 or value > 126:
            raise _line_error(
                source,
                line_number,
                f"quality character {char!r} is outside the ASCII 33-126 range",
            )


def _fastq_records(
    lines: Iterator[tuple[int, str]], source: str
) -> Iterator[SequenceRecord]:
    # State machine over logical lines.  Sequence and quality may each
    # span several physical lines; record boundaries are determined by the
    # '+' separator and the number of quality characters matching the
    # sequence length.
    header_line: str | None = None
    header_number = 0
    seq_parts: list[str] = []
    qual_segments: list[str] = []
    in_quality = False

    for number, raw_line in lines:
        line = raw_line.rstrip("\r\n")

        if header_line is None:
            if not line:
                raise _line_error(
                    source,
                    number,
                    "incomplete FASTQ record (expected a title line starting with '@')",
                )
            if not line.startswith("@"):
                raise _line_error(
                    source, number, "expected a FASTQ title line starting with '@'"
                )
            header_line = line
            header_number = number
            seq_parts = []
            qual_segments = []
            in_quality = False
            continue

        if not in_quality:
            if line.startswith("+"):
                plus_title = line[1:]
                if plus_title and plus_title != header_line[1:]:
                    raise _line_error(
                        source,
                        number,
                        "FASTQ separator title does not match the record title",
                    )
                if not _clean_sequence("".join(seq_parts)):
                    raise _line_error(source, header_number, "empty sequence")
                in_quality = True
            else:
                if not line.strip():
                    raise _line_error(source, number, "empty sequence")
                seq_parts.append(line)
            continue

        if not line:
            raise _line_error(source, number, "incomplete FASTQ record (missing quality values)")
        _check_quality_segment(line, source, number)
        qual_segments.append(line)

        sequence = _clean_sequence("".join(seq_parts))
        quality_text = "".join(qual_segments)
        if len(quality_text) > len(sequence):
            raise _line_error(
                source,
                number,
                f"quality length {len(quality_text)} does not match sequence length {len(sequence)}",
            )
        if len(quality_text) == len(sequence):
            yield _finish_fastq(
                header_line, header_number, sequence, quality_text, source
            )
            header_line = None

    if header_line is not None:
        raise _line_error(source, header_number, "incomplete FASTQ record")


def _finish_fastq(
    header_line: str,
    header_number: int,
    sequence: str,
    quality_text: str,
    source: str,
) -> SequenceRecord:
    identifier, description = _split_title(header_line[1:], source, header_number)
    if not sequence:
        raise _line_error(source, header_number, "empty sequence")
    _validate_sequence(identifier, sequence)
    if len(quality_text) != len(sequence):
        raise _line_error(
            source,
            header_number,
            f"quality length {len(quality_text)} does not match sequence length {len(sequence)}",
        )
    quality = tuple(ord(char) - 33 for char in quality_text)
    return SequenceRecord(identifier, sequence, description, quality)


def _enumerate_lines(stream: Iterable[str]) -> Iterator[tuple[int, str]]:
    for number, line in enumerate(stream, start=1):
        yield number, line


def read_sequences(
    source: "str | os.PathLike[str] | io.TextIOBase",
    format: str = "auto",
) -> Iterator[SequenceRecord]:
    """Lazily read sequence records from *source*.

    *source* is either a filesystem path (``str`` or :class:`os.PathLike`)
    or an already-open text stream.  *format* is ``"fasta"``, ``"fastq"``
    or ``"auto"``; with ``"auto"`` the format is detected from the first
    non-empty line (``>`` → FASTA, ``@`` → FASTQ).

    Records are produced in input order.  Errors raised while opening a
    path propagate as :class:`OSError`; malformed input raises
    :class:`SequenceFormatError` or :class:`SequenceValidationError`
    during iteration and are never swallowed.
    """
    if format not in _INPUT_FORMATS:
        raise ValueError(
            f"unsupported format {format!r}; expected one of {', '.join(_INPUT_FORMATS)}"
        )

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
        yield from _read_from_stream(stream, source_name, format)
    finally:
        if close:
            stream.close()


def _read_from_stream(
    stream: Iterable[str], source_name: str, format: str
) -> Iterator[SequenceRecord]:
    numbered = _enumerate_lines(stream)

    first_number = 0
    first_line: str | None = None
    for number, raw_line in numbered:
        candidate = raw_line.rstrip("\r\n")
        if candidate.strip():
            first_number = number
            first_line = candidate
            break

    if first_line is None:
        if format == "auto":
            raise SequenceFormatError(
                f"{source_name}: empty input; cannot determine format"
            )
        raise SequenceFormatError(f"{source_name}: empty input")

    if format == "auto":
        if first_line.startswith(">"):
            detected = "fasta"
        elif first_line.startswith("@"):
            detected = "fastq"
        else:
            raise _line_error(
                source_name,
                first_number,
                "cannot determine format: first non-empty line must start with '>' or '@'",
            )
    else:
        detected = format

    def replay() -> Iterator[tuple[int, str]]:
        yield first_number, first_line + "\n"
        yield from numbered

    # The per-format parsers reject a first line that lacks their title
    # marker themselves (reported as content before / unexpected title).
    if detected == "fasta":
        yield from _fasta_records(replay(), source_name)
    else:
        yield from _fastq_records(replay(), source_name)


def write_sequences(
    records: Iterable[SequenceRecord],
    output: "str | os.PathLike[str] | io.TextIOBase",
    format: str = "fasta",
    line_width: int = 60,
) -> None:
    """Write *records* to *output* in FASTA or FASTQ form.

    FASTA sequences are wrapped at *line_width* characters (a positive
    integer); FASTQ always uses the fixed four-line layout regardless of
    *line_width*.  Records without quality values cannot be written as
    FASTQ.  Output uses ``"\\n"`` line endings and ends with exactly one
    trailing newline.
    """
    if format not in _OUTPUT_FORMATS:
        raise ValueError(f"unsupported format {format!r}; expected 'fasta' or 'fastq'")
    if not isinstance(line_width, int) or isinstance(line_width, bool) or line_width <= 0:
        raise ValueError("line_width must be a positive integer")

    if isinstance(output, io.IOBase):
        stream: io.TextIOBase = output
        close = False
    elif isinstance(output, (str, os.PathLike)):
        stream = open(os.fspath(output), "w", encoding="utf-8", newline="")
        close = True
    else:
        raise TypeError("output must be a text path or a text stream")

    try:
        for record in records:
            _validate_sequence(record.identifier, record.sequence)
            if format == "fasta":
                stream.write(_format_fasta_record(record, line_width))
            else:
                stream.write(_format_fastq_record(record))
    finally:
        if close:
            stream.close()


def _format_title(record: SequenceRecord) -> str:
    if record.description:
        return record.identifier + " " + record.description
    return record.identifier


def _format_fasta_record(record: SequenceRecord, line_width: int) -> str:
    parts = [">" + _format_title(record) + "\n"]
    sequence = record.sequence
    for start in range(0, len(sequence), line_width):
        parts.append(sequence[start : start + line_width])
        parts.append("\n")
    return "".join(parts)


def _format_fastq_record(record: SequenceRecord) -> str:
    if record.quality is None:
        raise ValueError(
            f"record {record.identifier!r} has no quality values; cannot write it as FASTQ"
        )
    if len(record.quality) != len(record.sequence):
        raise ValueError(
            f"record {record.identifier!r}: quality length does not match sequence length"
        )
    quality_chars = []
    for value in record.quality:
        if not 0 <= value <= 93:
            raise ValueError(
                f"record {record.identifier!r}: quality score {value} is outside the "
                "Phred+33 range 0-93"
            )
        quality_chars.append(chr(value + 33))
    title = _format_title(record)
    return f"@{title}\n{record.sequence}\n+{title}\n{''.join(quality_chars)}\n"
