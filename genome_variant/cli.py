"""Command line entry point for genome-variant-toolkit."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable, Iterator, Sequence

from . import __version__
from .sequence_io import (
    SequenceFormatError,
    SequenceRecord,
    SequenceValidationError,
    read_sequences,
    write_sequences,
)


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r} is not a positive integer")
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"{value!r} is not a positive integer")
    return parsed


def _prefix(first: SequenceRecord, rest: Iterable[SequenceRecord]) -> Iterator[SequenceRecord]:
    yield first
    yield from rest


class _ArgumentParser(argparse.ArgumentParser):
    """Argument parser that reports errors as one line on stderr."""

    def error(self, message: str) -> None:
        sys.stderr.write(f"{self.prog}: error: {message}\n")
        raise SystemExit(2)


def _build_parser() -> argparse.ArgumentParser:
    parser = _ArgumentParser(
        prog="genome-variant-toolkit",
        description="Genomic variant analysis toolkit from reads to annotated calls",
    )
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("version", help="print the current version")

    normalize = sub.add_parser(
        "normalize-sequences",
        help="read FASTA/FASTQ sequences and write them in a normalized form",
    )
    normalize.add_argument("input", help="input file, or '-' for standard input")
    normalize.add_argument(
        "--input-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="input format (default: auto-detect from the first title line)",
    )
    normalize.add_argument(
        "--output-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="output format (default: keep the detected input format)",
    )
    normalize.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )
    normalize.add_argument(
        "--line-width",
        type=_positive_int,
        default=60,
        help="maximum FASTA sequence line width in characters (default: 60); "
        "FASTQ layout is always four lines per record",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "version":
        print(__version__)
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "normalize-sequences":
        return _run_normalize(args, parser)

    parser.print_help()  # pragma: no cover - every subcommand is handled above
    return 0


def _run_normalize(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.input == "-":
        input_stream = sys.stdin
        close_input = False
    else:
        try:
            input_stream = open(args.input, "r", encoding="utf-8", newline="")
        except OSError as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 1
        close_input = True

    try:
        records: Iterable[SequenceRecord] = read_sequences(
            input_stream, format=args.input_format
        )
        output_format = args.output_format

        if output_format == "auto":
            # read_sequences resolves "auto" and validates lazily; peek at
            # the first record to learn the input format.
            try:
                first = next(records)
            except (SequenceFormatError, SequenceValidationError) as exc:
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 2
            output_format = "fastq" if first.quality is not None else "fasta"
            records = _prefix(first, records)

        if args.output == "-":
            # Keep output byte-stable across platforms: no newline
            # translation on standard output.
            if hasattr(sys.stdout, "reconfigure"):
                sys.stdout.reconfigure(newline="")
            output_stream = sys.stdout
            close_output = False
        else:
            try:
                output_stream = open(args.output, "w", encoding="utf-8", newline="")
            except OSError as exc:
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 1
            close_output = True

        try:
            try:
                write_sequences(
                    records,
                    output_stream,
                    format=output_format,
                    line_width=args.line_width,
                )
            except (SequenceFormatError, SequenceValidationError, ValueError) as exc:
                # ValueError covers usage problems such as writing FASTA
                # records without quality values as FASTQ.
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 2
            except OSError as exc:
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 1
        finally:
            if close_output:
                try:
                    output_stream.close()
                except OSError as exc:
                    print(f"{parser.prog}: {exc}", file=sys.stderr)
                    return 1
    finally:
        if close_input:
            input_stream.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
