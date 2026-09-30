"""Command line entry point for genome-variant-toolkit."""

from __future__ import annotations

import argparse
import sys
from itertools import chain

from . import __version__
from .sequence_io import (
    AUTO,
    FASTA,
    FASTQ,
    SequenceFormatError,
    SequenceValidationError,
    read_sequences,
    write_sequences,
)


class _UsageError(Exception):
    """Raised instead of ``SystemExit`` so argument errors stay single-line."""


class _ArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:  # pragma: no cover - exercised via CLI
        raise _UsageError(f"{self.prog}: error: {message}")


def _positive_int(value: str) -> int:
    try:
        width = int(value)
    except ValueError:
        width = 0
    if width <= 0:
        raise argparse.ArgumentTypeError(
            f"line-width must be a positive integer, got {value!r}"
        )
    return width


def _build_parser() -> argparse.ArgumentParser:
    parser = _ArgumentParser(
        prog="genome-variant-toolkit",
        description=(
            "Genomic variant analysis toolkit from reads to annotated calls"
        ),
    )
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("version", help="print the current version")

    normalize = sub.add_parser(
        "normalize-sequences",
        help="read FASTA/FASTQ sequences and write them back in a canonical form",
    )
    normalize.add_argument("input", help="input file, or '-' for standard input")
    normalize.add_argument(
        "--input-format",
        choices=("fasta", "fastq", "auto"),
        default=AUTO,
        help="input format (default: auto-detect from the first non-empty line)",
    )
    normalize.add_argument(
        "--output-format",
        choices=("fasta", "fastq"),
        default=None,
        help="output format (default: keep the detected input format)",
    )
    normalize.add_argument(
        "--output",
        default=None,
        help="output file (default: standard output); '-' also means stdout",
    )
    normalize.add_argument(
        "--line-width",
        type=_positive_int,
        default=80,
        help="maximum FASTA sequence line width; does not affect FASTQ",
    )
    return parser


def _run_normalize(args: argparse.Namespace) -> int:
    source_name = "<stdin>" if args.input == "-" else args.input
    records = read_sequences(
        sys.stdin if args.input == "-" else args.input,
        format=args.input_format,
        source_name=source_name,
    )

    # Pull the first record eagerly: this surfaces empty input and format
    # errors before any output file is touched, and lets us pick the output
    # format when it was left implicit.
    try:
        first = next(records)
    except StopIteration:  # pragma: no cover - readers reject empty input
        raise SequenceFormatError(f"{source_name}:1: empty input")

    if args.output_format is not None:
        output_format = args.output_format
    else:
        output_format = FASTQ if first.quality is not None else FASTA

    if output_format == FASTQ and first.quality is None:
        raise ValueError(
            f"{source_name}: cannot write FASTQ output from FASTA input: "
            f"record {first.id!r} has no quality scores"
        )

    if not args.output or args.output == "-":
        # Never translate newlines: output is always LF, even on Windows.
        try:
            sys.stdout.reconfigure(newline="\n")
        except (AttributeError, ValueError):  # pragma: no cover - non-standard stream
            pass
        output = sys.stdout
    else:
        output = args.output

    write_sequences(
        chain([first], records),
        output,
        format=output_format,
        line_width=args.line_width,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    try:
        args = parser.parse_args(argv)
    except _UsageError as error:
        print(str(error), file=sys.stderr)
        return 2

    if args.command == "version":
        print(__version__)
        return 0

    if args.command == "normalize-sequences":
        try:
            return _run_normalize(args)
        except (SequenceFormatError, SequenceValidationError, ValueError) as error:
            print(str(error), file=sys.stderr)
            return 2
        except OSError as error:
            print(str(error) or error.__class__.__name__, file=sys.stderr)
            return 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
