"""Command line entry point for genome-variant-toolkit."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections.abc import Iterable, Iterator, Sequence

from . import __version__
from .kmer import build_kmer_index
from .quality import ReadQualityError, filter_reads
from .sequence_io import (
    SequenceFormatError,
    SequenceRecord,
    SequenceValidationError,
    read_sequences,
    write_sequences,
)

# Read/parse/processing errors that make the invocation invalid usage.
# ReadQualityError is a ValueError subclass but is listed explicitly.
_READ_ERRORS = (
    SequenceFormatError,
    SequenceValidationError,
    ReadQualityError,
    ValueError,
)


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r} is not a positive integer")
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"{value!r} is not a positive integer")
    return parsed


def _quality_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not an integer in the range 0-93"
        )
    if not 0 <= parsed <= 93:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not an integer in the range 0-93"
        )
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

    filter_reads_cmd = sub.add_parser(
        "filter-reads",
        help="trim low-quality ends and filter FASTQ reads",
    )
    filter_reads_cmd.add_argument("input", help="FASTQ input file, or '-' for standard input")
    filter_reads_cmd.add_argument(
        "--min-end-quality",
        type=_quality_int,
        default=20,
        metavar="N",
        help="trim bases with Phred quality strictly below N from both ends "
        "(default: 20)",
    )
    filter_reads_cmd.add_argument(
        "--min-mean-quality",
        type=_quality_int,
        default=20,
        metavar="N",
        help="drop reads whose mean Phred quality is below N after trimming "
        "(default: 20)",
    )
    filter_reads_cmd.add_argument(
        "--min-length",
        type=_positive_int,
        default=30,
        metavar="N",
        help="drop reads shorter than N bases after trimming (default: 30)",
    )
    filter_reads_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    kmer_index_cmd = sub.add_parser(
        "kmer-index",
        help="build a k-mer index and write it as JSON Lines",
    )
    kmer_index_cmd.add_argument("input", help="input file, or '-' for standard input")
    kmer_index_cmd.add_argument(
        "--input-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="input format (default: auto-detect from the first title line)",
    )
    kmer_index_cmd.add_argument(
        "--k",
        type=_positive_int,
        required=True,
        metavar="N",
        help="k-mer length (a positive integer)",
    )
    kmer_index_cmd.add_argument(
        "--no-canonical",
        action="store_true",
        help="index k-mers on the forward strand only (default: canonicalize "
        "each k-mer against its reverse complement)",
    )
    kmer_index_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
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

    if args.command == "filter-reads":
        return _run_filter_reads(args, parser)

    if args.command == "kmer-index":
        return _run_kmer_index(args, parser)

    if args.command == "kmer-index":
        return _run_kmer_index(args, parser)

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


def _run_filter_reads(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
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
            input_stream, format="fastq"
        )
        filtered = filter_reads(
            records,
            min_end_quality=args.min_end_quality,
            min_mean_quality=args.min_mean_quality,
            min_length=args.min_length,
        )

        if args.output == "-":
            # Keep output byte-stable across platforms: no newline
            # translation on standard output.
            if hasattr(sys.stdout, "reconfigure"):
                sys.stdout.reconfigure(newline="")
            output_stream = sys.stdout
            close_output = False
            temporary_path = None
        else:
            output_stream, temporary_path = _open_temporary_output(
                args.output, parser
            )
            if output_stream is None:
                return 1
            close_output = True

        try:
            try:
                write_sequences(filtered, output_stream, format="fastq")
            except _READ_ERRORS as exc:
                if temporary_path is not None:
                    _discard_temporary(output_stream, temporary_path)
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 2
            except OSError as exc:
                if temporary_path is not None:
                    _discard_temporary(output_stream, temporary_path)
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 1

            if temporary_path is not None:
                try:
                    output_stream.flush()
                    output_stream.close()
                    os.replace(temporary_path, args.output)
                except OSError as exc:
                    _discard_temporary(output_stream, temporary_path)
                    print(f"{parser.prog}: {exc}", file=sys.stderr)
                    return 1
        finally:
            # On failure the temporary file was already closed and removed
            # by _discard_temporary; on success it was closed explicitly
            # before os.replace. This only closes a still-open temp stream.
            if close_output and not output_stream.closed:
                try:
                    output_stream.close()
                except OSError:
                    pass
    finally:
        if close_input:
            input_stream.close()

    return 0


def _run_kmer_index(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
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
        # Build the whole index before writing anything: all reading and
        # validation happens here, so a failure never produces partial
        # output or touches the output file.
        try:
            index = build_kmer_index(
                records, args.k, canonical=not args.no_canonical
            )
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 1

        lines = []
        for kmer, occurrences in index.items():
            lines.append(
                json.dumps(
                    {
                        "kmer": kmer,
                        "count": len(occurrences),
                        "occurrences": [
                            {
                                "record": occurrence.record,
                                "id": occurrence.id,
                                "position": occurrence.position,
                                "strand": occurrence.strand,
                            }
                            for occurrence in occurrences
                        ],
                    },
                    ensure_ascii=False,
                    separators=(",", ":"),
                )
            )
        text = "".join(line + "\n" for line in lines)

        if args.output == "-":
            # Keep output byte-stable across platforms: no newline
            # translation on standard output.
            if hasattr(sys.stdout, "reconfigure"):
                sys.stdout.reconfigure(newline="")
            try:
                sys.stdout.write(text)
                sys.stdout.flush()
            except OSError as exc:
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 1
            return 0

        output_stream, temporary_path = _open_temporary_output(
            args.output, parser, prefix=".kmer-index-"
        )
        if output_stream is None:
            return 1
        try:
            try:
                output_stream.write(text)
                output_stream.flush()
                output_stream.close()
                os.replace(temporary_path, args.output)
            except OSError as exc:
                _discard_temporary(output_stream, temporary_path)
                print(f"{parser.prog}: {exc}", file=sys.stderr)
                return 1
        finally:
            if not output_stream.closed:
                try:
                    output_stream.close()
                except OSError:
                    pass
    finally:
        if close_input:
            input_stream.close()

    return 0


def _discard_temporary(stream: object, path: str) -> None:
    """Close and remove a temporary output file after a failed run."""
    close = getattr(stream, "close", None)
    if callable(close):
        try:
            close()
        except OSError:
            pass
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
    except OSError:
        pass


def _open_temporary_output(
    destination: str,
    parser: argparse.ArgumentParser,
    prefix: str = ".filter-reads-",
) -> tuple[object, str] | tuple[None, None]:
    """Create a temporary text file next to *destination*.

    Returns ``(stream, path)``; on failure prints one error line and
    returns ``(None, None)``.
    """
    output_dir = os.path.dirname(os.path.abspath(destination))
    try:
        fd, path = tempfile.mkstemp(
            dir=output_dir,
            prefix=prefix,
            suffix=".tmp",
        )
        # mkstemp creates files with mode 0600; match the mode a plain
        # open() would give the final output instead.
        current_umask = os.umask(0)
        os.umask(current_umask)
        os.fchmod(fd, 0o666 & ~current_umask)
    except OSError as exc:
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return None, None
    try:
        return os.fdopen(fd, "w", encoding="utf-8", newline=""), path
    except OSError as exc:
        try:
            os.close(fd)
        finally:
            try:
                os.unlink(path)
            except OSError:
                pass
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return None, None


if __name__ == "__main__":
    sys.exit(main())
