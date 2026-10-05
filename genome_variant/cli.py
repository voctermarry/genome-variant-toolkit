"""Command line entry point for genome-variant-toolkit."""

from __future__ import annotations

import argparse
import io
import json
import math
import os
import sys
import tempfile
from collections.abc import Iterable, Iterator, Sequence

from . import __version__
from .alignment import PairwiseAlignment, align_pair
from .annotation import AnnotationFormatError, annotate_vcf
from .batch import BatchManifestError, batch_call_variants
from .calling import VariantCallingError, call_variants
from .comparison import (
    VariantComparisonError,
    compare_variants,
    render_comparison,
)
from .coverage import coverage_report
from .deduplication import deduplicate_reads
from .kmer import build_kmer_index
from .mapping import MappingReferenceError, map_reads
from .quality import ReadQualityError, filter_reads
from .sequence_io import (
    SequenceFormatError,
    SequenceRecord,
    SequenceValidationError,
    read_sequences,
    write_sequences,
)
from .summary import (
    ManifestError,
    render_summaries,
    summarize_variants,
)
from .vcf import (
    ReferenceMismatchError,
    VcfFormatError,
    normalize_vcf,
    read_vcf,
    render_vcf,
)

# Read/parse/processing errors that make the invocation invalid usage.
# ReadQualityError is a ValueError subclass but is listed explicitly.
_READ_ERRORS = (
    SequenceFormatError,
    SequenceValidationError,
    ReadQualityError,
    VcfFormatError,
    ReferenceMismatchError,
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


def _nonnegative_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not a non-negative integer"
        )
    if parsed < 0:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not a non-negative integer"
        )
    return parsed


def _fraction(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not a number in the range 0-1"
        )
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not a finite number in the range 0-1"
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

    dedup_cmd = sub.add_parser(
        "deduplicate-reads",
        help="deduplicate FASTA/FASTQ reads before mapping and variant calling",
    )
    dedup_cmd.add_argument("input", help="input file, or '-' for standard input")
    dedup_cmd.add_argument(
        "--input-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="input format (default: auto-detect from the first title line)",
    )
    dedup_cmd.add_argument(
        "--canonical",
        action="store_true",
        help="also treat a sequence and its IUPAC reverse complement as "
        "duplicates (default: only identical sequences are grouped)",
    )
    dedup_cmd.add_argument(
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

    align_pair_cmd = sub.add_parser(
        "align-pair",
        help="align two single-record FASTA/FASTQ inputs pairwise",
    )
    align_pair_cmd.add_argument(
        "reference", help="reference input file, or '-' for standard input"
    )
    align_pair_cmd.add_argument(
        "query", help="query input file, or '-' for standard input"
    )
    align_pair_cmd.add_argument(
        "--reference-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="reference input format (default: auto-detect)",
    )
    align_pair_cmd.add_argument(
        "--query-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="query input format (default: auto-detect)",
    )
    align_pair_cmd.add_argument(
        "--mode",
        choices=("global", "local"),
        default="global",
        help="alignment mode (default: global)",
    )
    align_pair_cmd.add_argument(
        "--match",
        "--match-score",
        dest="match",
        type=_positive_int,
        default=2,
        metavar="N",
        help="match score, a positive integer (default: 2)",
    )
    align_pair_cmd.add_argument(
        "--mismatch",
        "--mismatch-penalty",
        dest="mismatch",
        type=_nonnegative_int,
        default=3,
        metavar="N",
        help="mismatch penalty, a non-negative integer (default: 3)",
    )
    align_pair_cmd.add_argument(
        "--gap-open",
        type=_nonnegative_int,
        default=5,
        metavar="N",
        help="gap opening penalty, a non-negative integer (default: 5)",
    )
    align_pair_cmd.add_argument(
        "--gap-extend",
        type=_nonnegative_int,
        default=2,
        metavar="N",
        help="gap extension penalty, a non-negative integer (default: 2)",
    )
    align_pair_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    normalize_vcf_cmd = sub.add_parser(
        "normalize-vcf",
        help="normalize VCF variants against a reference FASTA",
    )
    normalize_vcf_cmd.add_argument(
        "input", help="VCF input file, or '-' for standard input"
    )
    normalize_vcf_cmd.add_argument(
        "--reference",
        required=True,
        help="reference FASTA file (required); '-' denotes standard input",
    )
    normalize_vcf_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    annotate_vcf_cmd = sub.add_parser(
        "annotate-vcf",
        help="annotate VCF SNVs and short ACGT indels against a "
        "reference FASTA and GFF3 features",
    )
    annotate_vcf_cmd.add_argument(
        "input", help="VCF input file, or '-' for standard input"
    )
    annotate_vcf_cmd.add_argument(
        "--reference",
        required=True,
        help="reference FASTA file (required); '-' denotes standard input",
    )
    annotate_vcf_cmd.add_argument(
        "--features",
        required=True,
        help="GFF3 feature file (required); '-' denotes standard input",
    )
    annotate_vcf_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    map_reads_cmd = sub.add_parser(
        "map-reads",
        help="map FASTA/FASTQ reads against a multi-record FASTA reference",
    )
    map_reads_cmd.add_argument(
        "reference", help="reference FASTA file, or '-' for standard input"
    )
    map_reads_cmd.add_argument(
        "reads", help="reads FASTA/FASTQ file, or '-' for standard input"
    )
    map_reads_cmd.add_argument(
        "--reads-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="reads input format (default: auto-detect)",
    )
    map_reads_cmd.add_argument(
        "--match",
        "--match-score",
        dest="match",
        type=_positive_int,
        default=2,
        metavar="N",
        help="match score, a positive integer (default: 2)",
    )
    map_reads_cmd.add_argument(
        "--mismatch",
        "--mismatch-penalty",
        dest="mismatch",
        type=_nonnegative_int,
        default=3,
        metavar="N",
        help="mismatch penalty, a non-negative integer (default: 3)",
    )
    map_reads_cmd.add_argument(
        "--gap-open",
        type=_nonnegative_int,
        default=5,
        metavar="N",
        help="gap opening penalty, a non-negative integer (default: 5)",
    )
    map_reads_cmd.add_argument(
        "--gap-extend",
        type=_nonnegative_int,
        default=2,
        metavar="N",
        help="gap extension penalty, a non-negative integer (default: 2)",
    )
    map_reads_cmd.add_argument(
        "--min-score",
        type=_positive_int,
        default=1,
        metavar="N",
        help="minimum local-alignment score to consider mapped (default: 1)",
    )
    map_reads_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    call_variants_cmd = sub.add_parser(
        "call-variants",
        help="call single-sample SNVs from FASTQ reads against a FASTA reference",
    )
    call_variants_cmd.add_argument(
        "reference", help="reference FASTA file, or '-' for standard input"
    )
    call_variants_cmd.add_argument(
        "reads", help="reads FASTQ file, or '-' for standard input"
    )
    call_variants_cmd.add_argument(
        "--reads-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="reads input format (default: auto-detect); reads without "
        "quality values are a calling error",
    )
    call_variants_cmd.add_argument(
        "--min-base-quality",
        type=_quality_int,
        default=20,
        metavar="N",
        help="admit only bases with Phred quality at least N (default: 20)",
    )
    call_variants_cmd.add_argument(
        "--min-alt-count",
        type=_positive_int,
        default=2,
        metavar="N",
        help="minimum number of ALT observations to call a variant "
        "(default: 2)",
    )
    call_variants_cmd.add_argument(
        "--min-alt-fraction",
        type=_fraction,
        default=0.2,
        metavar="F",
        help="minimum ALT fraction AC/DP required to call a variant "
        "(default: 0.2)",
    )
    call_variants_cmd.add_argument(
        "--homozygous-fraction",
        type=_fraction,
        default=0.8,
        metavar="F",
        help="ALT fraction at or above which the genotype is 1/1 "
        "(default: 0.8); must not be below --min-alt-fraction",
    )
    call_variants_cmd.add_argument(
        "--sample-name",
        default="SAMPLE",
        help="sample name written to the #CHROM header (default: SAMPLE)",
    )
    call_variants_cmd.add_argument(
        "--call-indels",
        action="store_true",
        help="also call short insertions and deletions from the winning "
        "alignments (default: SNVs only)",
    )
    call_variants_cmd.add_argument(
        "--max-indel-length",
        type=_positive_int,
        default=None,
        metavar="N",
        help="longest insertion or deletion to call, a positive integer "
        "(default: 50); requires --call-indels",
    )
    call_variants_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    batch_call_cmd = sub.add_parser(
        "batch-call-variants",
        help="call variants for every sample of a JSON Lines reads manifest",
    )
    batch_call_cmd.add_argument(
        "manifest", help="JSON Lines manifest file, or '-' for standard input"
    )
    batch_call_cmd.add_argument(
        "--reference",
        required=True,
        help="reference FASTA file (required); '-' denotes standard input",
    )
    batch_call_cmd.add_argument(
        "--min-base-quality",
        type=_quality_int,
        default=20,
        metavar="N",
        help="admit only bases with Phred quality at least N (default: 20)",
    )
    batch_call_cmd.add_argument(
        "--min-alt-count",
        type=_positive_int,
        default=2,
        metavar="N",
        help="minimum number of ALT observations to call a variant "
        "(default: 2)",
    )
    batch_call_cmd.add_argument(
        "--min-alt-fraction",
        type=_fraction,
        default=0.2,
        metavar="F",
        help="minimum ALT fraction AC/DP required to call a variant "
        "(default: 0.2)",
    )
    batch_call_cmd.add_argument(
        "--homozygous-fraction",
        type=_fraction,
        default=0.8,
        metavar="F",
        help="ALT fraction at or above which the genotype is 1/1 "
        "(default: 0.8); must not be below --min-alt-fraction",
    )
    batch_call_cmd.add_argument(
        "--call-indels",
        action="store_true",
        help="also call short insertions and deletions from the winning "
        "alignments (default: SNVs only)",
    )
    batch_call_cmd.add_argument(
        "--max-indel-length",
        type=_positive_int,
        default=None,
        metavar="N",
        help="longest insertion or deletion to call, a positive integer "
        "(default: 50); requires --call-indels",
    )
    batch_call_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    coverage_report_cmd = sub.add_parser(
        "coverage-report",
        help="report per-reference coverage and concordance from mapped reads",
    )
    coverage_report_cmd.add_argument(
        "reference", help="reference FASTA file, or '-' for standard input"
    )
    coverage_report_cmd.add_argument(
        "reads", help="reads FASTA/FASTQ file, or '-' for standard input"
    )
    coverage_report_cmd.add_argument(
        "--reads-format",
        choices=("fasta", "fastq", "auto"),
        default="auto",
        help="reads input format (default: auto-detect)",
    )
    coverage_report_cmd.add_argument(
        "--match",
        "--match-score",
        dest="match",
        type=_positive_int,
        default=2,
        metavar="N",
        help="match score, a positive integer (default: 2)",
    )
    coverage_report_cmd.add_argument(
        "--mismatch",
        "--mismatch-penalty",
        dest="mismatch",
        type=_nonnegative_int,
        default=3,
        metavar="N",
        help="mismatch penalty, a non-negative integer (default: 3)",
    )
    coverage_report_cmd.add_argument(
        "--gap-open",
        type=_nonnegative_int,
        default=5,
        metavar="N",
        help="gap opening penalty, a non-negative integer (default: 5)",
    )
    coverage_report_cmd.add_argument(
        "--gap-extend",
        type=_nonnegative_int,
        default=2,
        metavar="N",
        help="gap extension penalty, a non-negative integer (default: 2)",
    )
    coverage_report_cmd.add_argument(
        "--min-score",
        type=_positive_int,
        default=1,
        metavar="N",
        help="minimum local-alignment score to consider mapped (default: 1)",
    )
    coverage_report_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    summarize_cmd = sub.add_parser(
        "summarize-variants",
        help="summarize single-sample VCFs listed by a JSON Lines manifest",
    )
    summarize_cmd.add_argument(
        "manifest", help="JSON Lines manifest file, or '-' for standard input"
    )
    summarize_cmd.add_argument(
        "--reference",
        required=True,
        help="reference FASTA file (required); '-' denotes standard input",
    )
    summarize_cmd.add_argument(
        "--output",
        default="-",
        help="output file, or '-' for standard output (default: standard output)",
    )

    compare_cmd = sub.add_parser(
        "compare-variants",
        help="compare the call sets of two single-sample VCFs against a "
        "reference FASTA",
    )
    compare_cmd.add_argument(
        "baseline", help="baseline VCF file, or '-' for standard input"
    )
    compare_cmd.add_argument(
        "candidate", help="candidate VCF file, or '-' for standard input"
    )
    compare_cmd.add_argument(
        "--reference",
        required=True,
        help="reference FASTA file (required); '-' denotes standard input",
    )
    compare_cmd.add_argument(
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

    if args.command == "deduplicate-reads":
        return _run_deduplicate_reads(args, parser)

    if args.command == "kmer-index":
        return _run_kmer_index(args, parser)

    if args.command == "align-pair":
        return _run_align_pair(args, parser)

    if args.command == "normalize-vcf":
        return _run_normalize_vcf(args, parser)

    if args.command == "annotate-vcf":
        return _run_annotate_vcf(args, parser)

    if args.command == "map-reads":
        return _run_map_reads(args, parser)

    if args.command == "call-variants":
        return _run_call_variants(args, parser)

    if args.command == "batch-call-variants":
        return _run_batch_call_variants(args, parser)

    if args.command == "coverage-report":
        return _run_coverage_report(args, parser)

    if args.command == "summarize-variants":
        return _run_summarize_variants(args, parser)

    if args.command == "compare-variants":
        return _run_compare_variants(args, parser)

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
            return _fail(parser, exc)
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
            except OSError as exc:
                return _fail(parser, exc)
            output_format = "fastq" if first.quality is not None else "fasta"
            records = _prefix(first, records)

        # Render the whole result before touching the output path, so a
        # late read/validation failure never truncates an existing target
        # or creates the final file before the run has succeeded.
        try:
            text = _render_sequences(
                records, output_format, line_width=args.line_width
            )
        except (SequenceFormatError, SequenceValidationError, ValueError) as exc:
            # ValueError covers usage problems such as writing FASTA
            # records without quality values as FASTQ.
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        return _emit_text(
            text, args.output, parser, prefix=".normalize-sequences-"
        )
    finally:
        if close_input:
            input_stream.close()


def _run_filter_reads(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.input == "-":
        input_stream = sys.stdin
        close_input = False
    else:
        try:
            input_stream = open(args.input, "r", encoding="utf-8", newline="")
        except OSError as exc:
            return _fail(parser, exc)
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

        # Read, trim/filter and render the complete result before the
        # output path is touched, so a late malformed record fails the
        # run without truncating an existing target.
        try:
            text = _render_sequences(filtered, "fastq")
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        return _emit_text(text, args.output, parser, prefix=".filter-reads-")
    finally:
        if close_input:
            input_stream.close()


def _run_deduplicate_reads(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> int:
    if args.input == "-":
        input_stream = sys.stdin
        close_input = False
    else:
        try:
            input_stream = open(args.input, "r", encoding="utf-8", newline="")
        except OSError as exc:
            return _fail(parser, exc)
        close_input = True

    try:
        records: Iterable[SequenceRecord] = read_sequences(
            input_stream, format=args.input_format
        )
        # Resolve the output format from the first record: output always
        # mirrors the input format.  Peeking here also makes sure an empty
        # input is reported as a format error (existing reader semantics)
        # rather than as empty success.
        try:
            first = next(records)
        except (SequenceFormatError, SequenceValidationError) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)
        output_format = "fastq" if first.quality is not None else "fasta"
        records = _prefix(first, records)

        # Deduplicate (validating every record) before any output is
        # prepared, so a failure never creates partial output or touches
        # an existing output file.
        try:
            unique_records = deduplicate_reads(
                records, canonical=args.canonical
            )
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        try:
            text = _render_sequences(unique_records, output_format)
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        return _emit_text(
            text, args.output, parser, prefix=".deduplicate-reads-"
        )
    finally:
        if close_input:
            input_stream.close()


def _run_kmer_index(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.input == "-":
        input_stream = sys.stdin
        close_input = False
    else:
        try:
            input_stream = open(args.input, "r", encoding="utf-8", newline="")
        except OSError as exc:
            return _fail(parser, exc)
        close_input = True

    try:
        records: Iterable[SequenceRecord] = read_sequences(
            input_stream, format=args.input_format
        )
        # Build the whole index before preparing any output: all reading
        # and validation happens here, so a failure never produces a
        # partial output or touches the output file.
        try:
            index = build_kmer_index(
                records, args.k, canonical=not args.no_canonical
            )
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

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

        return _emit_text(text, args.output, parser, prefix=".kmer-index-")
    finally:
        if close_input:
            input_stream.close()


def _read_single_record(
    path: str,
    fmt: str,
    parser: argparse.ArgumentParser,
    label: str,
) -> tuple[SequenceRecord | None, object, bool, int]:
    """Read exactly one record from *path*.

    Returns ``(record, stream, close, exit_code)``; *record* is ``None``
    when the invocation failed and *exit_code* is then ``1`` or ``2``.
    ``"-"`` reads from standard input.
    """
    if path == "-":
        stream = sys.stdin
        close = False
    else:
        try:
            stream = open(path, "r", encoding="utf-8", newline="")
        except OSError as exc:
            print(f"{parser.prog}: {label}: {exc}", file=sys.stderr)
            return None, None, False, 1
        close = True

    try:
        records = read_sequences(stream, format=fmt)
        try:
            first = next(records)
        except StopIteration:
            raise SequenceFormatError(f"{label}: empty input")
        try:
            next(records)
        except StopIteration:
            return first, stream, close, 0
        raise SequenceFormatError(f"{label}: expected exactly one record, got more")
    except (SequenceFormatError, SequenceValidationError, ValueError) as exc:
        if close:
            stream.close()
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return None, None, False, 2
    except OSError as exc:
        if close:
            stream.close()
        print(f"{parser.prog}: {label}: {exc}", file=sys.stderr)
        return None, None, False, 1


def _alignment_to_json(result: PairwiseAlignment) -> str:
    return json.dumps(
        {
            "reference": result.reference,
            "query": result.query,
            "mode": result.mode,
            "score": result.score,
            "reference_start": result.reference_start,
            "reference_end": result.reference_end,
            "query_start": result.query_start,
            "query_end": result.query_end,
            "cigar": result.cigar,
            "aligned_reference": result.aligned_reference,
            "aligned_query": result.aligned_query,
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _run_align_pair(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.reference == "-" and args.query == "-":
        print(
            f"{parser.prog}: reference and query cannot both be read from standard input",
            file=sys.stderr,
        )
        return 2

    reference_record, ref_stream, close_ref, ref_code = _read_single_record(
        args.reference, args.reference_format, parser, "reference"
    )
    if reference_record is None:
        return ref_code

    query_record, query_stream, close_query, query_code = _read_single_record(
        args.query, args.query_format, parser, "query"
    )
    if query_record is None:
        if close_ref:
            ref_stream.close()
        return query_code

    try:
        try:
            result = align_pair(
                reference_record,
                query_record,
                mode=args.mode,
                match_score=args.match,
                mismatch_penalty=args.mismatch,
                gap_open=args.gap_open,
                gap_extend=args.gap_extend,
            )
        except (SequenceValidationError, ValueError) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        text = _alignment_to_json(result) + "\n"

        return _emit_text(text, args.output, parser, prefix=".align-pair-")
    finally:
        if close_ref:
            ref_stream.close()
        if close_query:
            query_stream.close()


def _open_input(path: str, parser: argparse.ArgumentParser) -> tuple[object, bool, int]:
    """Open *path* for reading, supporting '-' for standard input.

    Returns ``(stream, close, exit_code)``; on failure *stream* is
    ``None`` and *exit_code* is ``1``.
    """
    if path == "-":
        return sys.stdin, False, 0
    try:
        return open(path, "r", encoding="utf-8", newline=""), True, 0
    except OSError as exc:
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return None, False, 1


def _run_normalize_vcf(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.input == "-" and args.reference == "-":
        print(
            f"{parser.prog}: VCF input and reference cannot both be read "
            "from standard input",
            file=sys.stderr,
        )
        return 2

    vcf_stream, close_vcf, code = _open_input(args.input, parser)
    if vcf_stream is None:
        return code
    ref_stream, close_ref, ref_code = _open_input(args.reference, parser)
    if ref_stream is None:
        if close_vcf:
            vcf_stream.close()
        return ref_code

    try:
        try:
            document = read_vcf(vcf_stream)
        except VcfFormatError as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        try:
            reference_records = list(read_sequences(ref_stream, format="fasta"))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 1

        # Normalization validates every record (including REF checks)
        # and the full VCF is rendered before any output is prepared, so
        # a failure leaves an existing target untouched.
        try:
            normalized = normalize_vcf(document, reference_records)
        except (VcfFormatError, ReferenceMismatchError) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2

        text = render_vcf(normalized)

        return _emit_text(text, args.output, parser, prefix=".normalize-vcf-")
    finally:
        if close_vcf:
            vcf_stream.close()
        if close_ref:
            ref_stream.close()


def _run_annotate_vcf(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    stdin_inputs = [
        name
        for name, value in (
            ("input", args.input),
            ("reference", args.reference),
            ("features", args.features),
        )
        if value == "-"
    ]
    if len(stdin_inputs) > 1:
        print(
            f"{parser.prog}: at most one of input, --reference and --features "
            "can be read from standard input",
            file=sys.stderr,
        )
        return 2

    vcf_stream, close_vcf, code = _open_input(args.input, parser)
    if vcf_stream is None:
        return code
    ref_stream, close_ref, ref_code = _open_input(args.reference, parser)
    if ref_stream is None:
        if close_vcf:
            vcf_stream.close()
        return ref_code
    feat_stream, close_feat, feat_code = _open_input(args.features, parser)
    if feat_stream is None:
        if close_vcf:
            vcf_stream.close()
        if close_ref:
            ref_stream.close()
        return feat_code

    try:
        # Parse and validate every input completely before preparing the
        # output, so a failure never creates partial output or replaces an
        # existing target.
        try:
            document = read_vcf(vcf_stream)
        except VcfFormatError as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        try:
            reference_records = list(read_sequences(ref_stream, format="fasta"))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 1

        try:
            annotated = annotate_vcf(document, reference_records, feat_stream)
        except (
            AnnotationFormatError,
            VcfFormatError,
            ReferenceMismatchError,
        ) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: features: {exc}", file=sys.stderr)
            return 1

        text = render_vcf(annotated)

        return _emit_text(text, args.output, parser, prefix=".annotate-vcf-")
    finally:
        if close_vcf:
            vcf_stream.close()
        if close_ref:
            ref_stream.close()
        if close_feat:
            feat_stream.close()


def _mapping_to_json(mapping) -> str:
    return json.dumps(
        mapping.to_dict(),
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _run_map_reads(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.reference == "-" and args.reads == "-":
        print(
            f"{parser.prog}: reference and reads cannot both be read from standard input",
            file=sys.stderr,
        )
        return 2

    ref_stream, close_ref, ref_code = _open_input(args.reference, parser)
    if ref_stream is None:
        return ref_code
    read_stream, close_reads, reads_code = _open_input(args.reads, parser)
    if read_stream is None:
        if close_ref:
            ref_stream.close()
        return reads_code

    try:
        # Read the whole reference and validate every read, compute all
        # mappings and render the complete result before any output is
        # prepared: a failure never creates a partial output or touches an
        # existing output file.
        try:
            reference_records = list(read_sequences(ref_stream, format="fasta"))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 1

        try:
            read_records = list(read_sequences(read_stream, format=args.reads_format))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 1

        try:
            mappings = map_reads(
                reference_records,
                read_records,
                match_score=args.match,
                mismatch_penalty=args.mismatch,
                gap_open=args.gap_open,
                gap_extend=args.gap_extend,
                min_score=args.min_score,
            )
            lines = [_mapping_to_json(mapping) for mapping in mappings]
        except (MappingReferenceError, SequenceValidationError, ValueError) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2

        text = "".join(line + "\n" for line in lines)

        return _emit_text(text, args.output, parser, prefix=".map-reads-")
    finally:
        if close_ref:
            ref_stream.close()
        if close_reads:
            read_stream.close()


def _run_call_variants(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.reference == "-" and args.reads == "-":
        print(
            f"{parser.prog}: reference and reads cannot both be read from standard input",
            file=sys.stderr,
        )
        return 2

    if args.homozygous_fraction < args.min_alt_fraction:
        print(
            f"{parser.prog}: --homozygous-fraction must not be smaller than "
            "--min-alt-fraction",
            file=sys.stderr,
        )
        return 2

    if args.max_indel_length is not None and not args.call_indels:
        print(
            f"{parser.prog}: --max-indel-length requires --call-indels",
            file=sys.stderr,
        )
        return 2

    ref_stream, close_ref, ref_code = _open_input(args.reference, parser)
    if ref_stream is None:
        return ref_code
    read_stream, close_reads, reads_code = _open_input(args.reads, parser)
    if read_stream is None:
        if close_ref:
            ref_stream.close()
        return reads_code

    try:
        # Validate thresholds before consuming either input; reading and
        # calling must finish (and the VCF be rendered) completely before
        # the output is prepared, so a failure preserves an existing
        # target.  The reads stream is consumed lazily by call_variants in
        # a single pass: a malformed record fails the call when that
        # record is reached, before any output is produced.
        try:
            reference_records = list(read_sequences(ref_stream, format="fasta"))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 1

        reads = read_sequences(read_stream, format=args.reads_format)
        try:
            document = call_variants(
                reference_records,
                reads,
                min_base_quality=args.min_base_quality,
                min_alt_count=args.min_alt_count,
                min_alt_fraction=args.min_alt_fraction,
                homozygous_fraction=args.homozygous_fraction,
                sample_name=args.sample_name,
                call_indels=args.call_indels,
                max_indel_length=(
                    50 if args.max_indel_length is None else args.max_indel_length
                ),
            )
            text = render_vcf(document)
        except (
            VariantCallingError,
            ReadQualityError,
            SequenceValidationError,
            ValueError,
        ) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            return _fail(parser, exc)

        return _emit_text(text, args.output, parser, prefix=".call-variants-")
    finally:
        if close_ref:
            ref_stream.close()
        if close_reads:
            read_stream.close()


def _run_batch_call_variants(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> int:
    if args.manifest == "-" and args.reference == "-":
        print(
            f"{parser.prog}: manifest and reference cannot both be read "
            "from standard input",
            file=sys.stderr,
        )
        return 2

    if args.homozygous_fraction < args.min_alt_fraction:
        print(
            f"{parser.prog}: --homozygous-fraction must not be smaller than "
            "--min-alt-fraction",
            file=sys.stderr,
        )
        return 2

    if args.max_indel_length is not None and not args.call_indels:
        print(
            f"{parser.prog}: --max-indel-length requires --call-indels",
            file=sys.stderr,
        )
        return 2

    # Pass paths through unchanged so a file manifest's relative reads
    # paths resolve against the manifest's directory; "-" becomes the
    # standard input stream (relative paths then use the working dir).
    manifest_source = sys.stdin if args.manifest == "-" else args.manifest
    reference_source = sys.stdin if args.reference == "-" else args.reference

    # Read and validate the manifest, every sample's reads and the
    # reference completely, call every sample and build the full merged
    # result before preparing the output: a failure never produces
    # partial output or replaces an existing target.
    try:
        summaries = batch_call_variants(
            manifest_source,
            reference_source,
            min_base_quality=args.min_base_quality,
            min_alt_count=args.min_alt_count,
            min_alt_fraction=args.min_alt_fraction,
            homozygous_fraction=args.homozygous_fraction,
            call_indels=args.call_indels,
            max_indel_length=(
                50 if args.max_indel_length is None else args.max_indel_length
            ),
        )
    except (
        BatchManifestError,
        SequenceFormatError,
        SequenceValidationError,
        ReadQualityError,
        VariantCallingError,
        ValueError,
    ) as exc:
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        return _fail(parser, exc)

    text = render_summaries(summaries)

    return _emit_text(
        text, args.output, parser, prefix=".batch-call-variants-"
    )


def _run_coverage_report(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.reference == "-" and args.reads == "-":
        print(
            f"{parser.prog}: reference and reads cannot both be read from standard input",
            file=sys.stderr,
        )
        return 2

    ref_stream, close_ref, ref_code = _open_input(args.reference, parser)
    if ref_stream is None:
        return ref_code
    read_stream, close_reads, reads_code = _open_input(args.reads, parser)
    if read_stream is None:
        if close_ref:
            ref_stream.close()
        return reads_code

    try:
        # Read the whole reference, validate every read and render the
        # complete report before any output is prepared: a failure never
        # creates a partial output or touches an existing output file.
        try:
            reference_records = list(read_sequences(ref_stream, format="fasta"))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: reference: {exc}", file=sys.stderr)
            return 1

        try:
            read_records = list(read_sequences(read_stream, format=args.reads_format))
        except _READ_ERRORS as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2
        except OSError as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 1

        try:
            report = coverage_report(
                reference_records,
                read_records,
                match_score=args.match,
                mismatch_penalty=args.mismatch,
                gap_open=args.gap_open,
                gap_extend=args.gap_extend,
                min_score=args.min_score,
            )
            lines = [
                json.dumps(
                    entry.to_dict(),
                    ensure_ascii=False,
                    separators=(",", ":"),
                )
                for entry in report
            ]
        except (MappingReferenceError, SequenceValidationError, ValueError) as exc:
            print(f"{parser.prog}: {exc}", file=sys.stderr)
            return 2

        text = "".join(line + "\n" for line in lines)

        return _emit_text(text, args.output, parser, prefix=".coverage-report-")
    finally:
        if close_ref:
            ref_stream.close()
        if close_reads:
            read_stream.close()


def _run_summarize_variants(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.manifest == "-" and args.reference == "-":
        print(
            f"{parser.prog}: manifest and reference cannot both be read "
            "from standard input",
            file=sys.stderr,
        )
        return 2

    # Pass paths through unchanged so a file manifest's relative VCF
    # paths resolve against the manifest's directory; "-" becomes the
    # standard input stream (relative paths then use the working dir).
    manifest_source = sys.stdin if args.manifest == "-" else args.manifest
    reference_source = sys.stdin if args.reference == "-" else args.reference

    # Read and validate the manifest, every VCF and the reference
    # completely, and build the full result, before preparing the output:
    # a failure never produces partial output or replaces an existing
    # target.
    try:
        summaries = summarize_variants(manifest_source, reference_source)
    except (
        ManifestError,
        VcfFormatError,
        ReferenceMismatchError,
        SequenceFormatError,
        SequenceValidationError,
    ) as exc:
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        return _fail(parser, exc)

    text = render_summaries(summaries)

    return _emit_text(
        text, args.output, parser, prefix=".summarize-variants-"
    )


def _run_compare_variants(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    stdin_inputs = [
        name
        for name, value in (
            ("baseline", args.baseline),
            ("candidate", args.candidate),
            ("--reference", args.reference),
        )
        if value == "-"
    ]
    if len(stdin_inputs) > 1:
        print(
            f"{parser.prog}: at most one of baseline, candidate and --reference "
            "can be read from standard input",
            file=sys.stderr,
        )
        return 2

    # Pass paths through unchanged; "-" becomes the standard input stream.
    baseline_source = sys.stdin if args.baseline == "-" else args.baseline
    candidate_source = sys.stdin if args.candidate == "-" else args.candidate
    reference_source = sys.stdin if args.reference == "-" else args.reference

    # Read and validate both VCFs and the reference completely and build
    # the full comparison before preparing the output: a failure never
    # produces partial output or replaces an existing target.
    try:
        comparison = compare_variants(
            baseline_source, candidate_source, reference_source
        )
    except (
        VariantComparisonError,
        VcfFormatError,
        ReferenceMismatchError,
        SequenceFormatError,
        SequenceValidationError,
    ) as exc:
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        return _fail(parser, exc)

    text = render_comparison(comparison)

    return _emit_text(text, args.output, parser, prefix=".compare-variants-")


def _fail(parser: argparse.ArgumentParser, exc: BaseException) -> int:
    """Print one error line for *exc* and return the I/O failure code."""
    print(f"{parser.prog}: {exc}", file=sys.stderr)
    return 1


def _write_stdout(text: str, parser: argparse.ArgumentParser) -> int:
    """Write fully materialized *text* to standard output.

    Standard output never uses a temporary file.  Newline translation is
    disabled so output is byte-stable across platforms.
    """
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="")
    try:
        sys.stdout.write(text)
        sys.stdout.flush()
    except OSError as exc:
        return _fail(parser, exc)
    return 0


def _discard_temporary(stream: object, path: str) -> None:
    """Close and remove a temporary output file after a failed run.

    Secondary errors during cleanup are swallowed so they can never mask
    the original failure's return code or error message.
    """
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


def _commit_output(
    text: str,
    destination: str,
    parser: argparse.ArgumentParser,
    prefix: str,
) -> int:
    """Atomically publish already computed *text* at *destination*.

    The transaction covers every output entry point: a temporary file is
    prepared in the destination's directory only after all inputs have
    been read and the full result has been computed, and the destination
    is replaced exactly once after the whole UTF-8 payload has been
    written, flushed and closed.  Any failure (temporary file creation,
    write, flush, close or the final replacement) leaves an existing
    destination untouched and removes the temporary file; errors raised
    during that cleanup never override the reported failure.
    """
    output_dir = os.path.dirname(os.path.abspath(destination))
    try:
        fd, temporary_path = tempfile.mkstemp(
            dir=output_dir,
            prefix=prefix,
            suffix=".tmp",
        )
    except OSError as exc:
        return _fail(parser, exc)
    try:
        # mkstemp creates files with mode 0600; match the mode a plain
        # open() would give the final output instead.
        current_umask = os.umask(0)
        os.umask(current_umask)
        os.fchmod(fd, 0o666 & ~current_umask)
    except OSError as exc:
        try:
            os.close(fd)
        finally:
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
        return _fail(parser, exc)

    try:
        stream = os.fdopen(fd, "w", encoding="utf-8", newline="")
    except OSError as exc:
        try:
            os.close(fd)
        finally:
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
        return _fail(parser, exc)

    try:
        try:
            stream.write(text)
            stream.flush()
            stream.close()
        except OSError as exc:
            _discard_temporary(stream, temporary_path)
            return _fail(parser, exc)
        try:
            os.replace(temporary_path, destination)
        except OSError as exc:
            _discard_temporary(stream, temporary_path)
            return _fail(parser, exc)
    finally:
        # Normally the stream was closed above (explicitly before the
        # replacement, or by _discard_temporary on failure); this only
        # closes a stream that is somehow still open.
        if not stream.closed:
            try:
                stream.close()
            except OSError:
                pass
    return 0


def _emit_text(
    text: str,
    destination: str,
    parser: argparse.ArgumentParser,
    prefix: str,
) -> int:
    """Publish fully materialized *text* per the shared output rules.

    ``"-"`` means standard output (no temporary file); any other path is
    committed atomically in its containing directory.
    """
    if destination == "-":
        return _write_stdout(text, parser)
    return _commit_output(text, destination, parser, prefix)


def _render_sequences(
    records: Iterable[SequenceRecord],
    output_format: str,
    line_width: int | None = None,
) -> str:
    """Render sequence *records* to their complete UTF-8 output text.

    Rendering into a buffer means all remaining read/validation errors
    surface before any output file is touched, and the bytes are
    identical whether they end up on standard output or in a file.
    """
    buffer = io.StringIO()
    if line_width is None:
        write_sequences(records, buffer, format=output_format)
    else:
        write_sequences(
            records, buffer, format=output_format, line_width=line_width
        )
    return buffer.getvalue()


if __name__ == "__main__":
    sys.exit(main())
