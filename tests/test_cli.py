"""Tests for the genome-variant-toolkit command line interface."""

from __future__ import annotations

import sys

import pytest

from genome_variant.cli import main


def run(argv, stdin: str = ""):
    from io import StringIO

    old_in, old_out, old_err = sys.stdin, sys.stdout, sys.stderr
    sys.stdin = StringIO(stdin)
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    try:
        try:
            code = main(argv)
        except SystemExit as exc:
            code = exc.code
    finally:
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdin, sys.stdout, sys.stderr = old_in, old_out, old_err
    return code, out, err


def test_no_subcommand_prints_help_and_succeeds() -> None:
    code, out, err = run([])
    assert code == 0
    assert "usage:" in out
    assert err == ""


def test_version() -> None:
    code, out, err = run(["version"])
    assert code == 0
    assert out.strip() == "0.1.0"
    assert err == ""


def test_normalize_fasta_stdin_wrapping() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--line-width", "4"],
        ">a gene desc\nacgt\naCGT\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">a gene desc\nACGT\nACGT\n"


def test_normalize_fastq_fixed_layout_ignores_width() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--line-width", "2"],
        "@r1\nAC\nGT\n+r1\nIIII\n@r2\nA\n+\n!\n",
    )
    assert code == 0
    assert out == "@r1\nACGT\n+r1\nIIII\n@r2\nA\n+r2\n!\n"


def test_fasta_to_fastq_is_usage_error() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--output-format", "fastq"],
        ">a\nACGT\n",
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert "no quality" in err


def test_fastq_to_fasta_conversion() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--output-format", "fasta"],
        "@x\nacgt\n+x\nIIII\n",
    )
    assert code == 0
    assert out == ">x\nACGT\n"


def test_bad_line_width_exit_2_single_line() -> None:
    for value in ("0", "-3", "abc"):
        code, out, err = run(["normalize-sequences", "-", "--line-width", value], ">a\nA\n")
        assert code == 2
        assert err.count("\n") == 1
        assert out == ""


def test_format_errors_exit_2() -> None:
    code, _, err = run(["normalize-sequences", "-"], "")
    assert code == 2 and err.endswith("\n") and err.count("\n") == 1

    code, _, err = run(["normalize-sequences", "-"], ">a\nACGZ\n")
    assert code == 2 and "position 4" in err

    code, _, err = run(["normalize-sequences", "-"], "@a\nACGT\n+\nIII\n")
    assert code == 2 and "incomplete" in err


def test_input_file_error_exit_1(tmp_path) -> None:
    code, out, err = run(["normalize-sequences", str(tmp_path / "missing.fa")])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_output_file_error_exit_1(tmp_path) -> None:
    code, _, err = run(
        [
            "normalize-sequences",
            "-",
            "--output",
            str(tmp_path / "no-such-dir" / "out.fa"),
        ],
        ">a\nACGT\n",
    )
    assert code == 1


def test_file_roundtrip_is_byte_stable(tmp_path) -> None:
    source = tmp_path / "in.fq"
    source.write_text("@x\nacgt\n+x\nIIII\n")
    first = tmp_path / "out1.fa"
    second = tmp_path / "out2.fa"
    code1, _, err1 = run(
        ["normalize-sequences", str(source), "--output", str(first)]
    )
    code2, _, err2 = run(
        ["normalize-sequences", str(source), "--output", str(second)]
    )
    assert code1 == code2 == 0
    assert err1 == err2 == ""
    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes() == b"@x\nACGT\n+x\nIIII\n"


def test_unknown_argument_exit_2_single_line() -> None:
    code, out, err = run(["normalize-sequences", "-", "--bogus"], ">a\nA\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_filter_reads_trims_and_filters_stdin() -> None:
    fastq = (
        "@keep desc\n" + "A" * 30 + "\n+keep desc\n" + "I" * 30 + "\n"
        "@trim\n" + "TT" + "A" * 30 + "TT\n+\n" + "!!" + "I" * 30 + "!!\n"
        "@low\n" + "A" * 30 + "\n+\n" + "!" * 30 + "\n"
    )
    code, out, err = run(["filter-reads", "-"], fastq)
    assert code == 0
    assert err == ""
    assert out == (
        "@keep desc\n" + "A" * 30 + "\n+keep desc\n" + "I" * 30 + "\n"
        "@trim\n" + "A" * 30 + "\n+trim\n" + "I" * 30 + "\n"
    )


def test_filter_reads_threshold_options() -> None:
    fastq = "@r\n" + "A" * 10 + "\n+\n" + "5" * 10 + "\n"  # Phred 20
    code, out, err = run(
        ["filter-reads", "-", "--min-length", "5", "--min-end-quality", "21"],
        fastq,
    )
    assert code == 0
    assert out == ""
    assert err == ""

    code, out, err = run(
        ["filter-reads", "-", "--min-length", "5", "--min-mean-quality", "21"],
        fastq,
    )
    assert code == 0
    assert out == ""


def test_filter_reads_all_filtered_output_empty() -> None:
    code, out, err = run(
        ["filter-reads", "-", "--min-length", "100"],
        "@r\n" + "A" * 10 + "\n+\n" + "I" * 10 + "\n",
    )
    assert code == 0
    assert out == ""
    assert err == ""


def test_filter_reads_rejects_fasta_input() -> None:
    code, out, err = run(["filter-reads", "-"], ">a\nACGT\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_filter_reads_format_and_base_errors_exit_2() -> None:
    code, out, err = run(["filter-reads", "-"], "@a\nACGZ\n+\nIIII\n")
    assert code == 2
    assert err.count("\n") == 1

    code, out, err = run(["filter-reads", "-"], "@a\nACGT\n+\nIII\n")
    assert code == 2
    assert err.count("\n") == 1


def test_filter_reads_invalid_thresholds_exit_2() -> None:
    fastq = "@r\n" + "A" * 30 + "\n+\n" + "I" * 30 + "\n"
    for args in (
        ["--min-end-quality", "94"],
        ["--min-end-quality", "-1"],
        ["--min-mean-quality", "94"],
        ["--min-length", "0"],
        ["--min-length", "abc"],
    ):
        code, out, err = run(["filter-reads", "-", *args], fastq)
        assert code == 2, args
        assert out == ""
        assert err.count("\n") == 1


def test_filter_reads_input_file_error_exit_1(tmp_path) -> None:
    code, out, err = run(["filter-reads", str(tmp_path / "missing.fq")])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_filter_reads_output_file_byte_stable(tmp_path) -> None:
    source = tmp_path / "in.fq"
    source.write_text("@x some read\n" + "a" * 30 + "\n+x some read\n" + "I" * 30 + "\n")
    first = tmp_path / "out1.fq"
    second = tmp_path / "out2.fq"
    code1, _, err1 = run(["filter-reads", str(source), "--output", str(first)])
    code2, _, err2 = run(["filter-reads", str(source), "--output", str(second)])
    assert code1 == code2 == 0
    assert err1 == err2 == ""
    expected = b"@x some read\n" + b"A" * 30 + b"\n+x some read\n" + b"I" * 30 + b"\n"
    assert first.read_bytes() == second.read_bytes() == expected


def test_filter_reads_output_error_leaves_existing_file(tmp_path) -> None:
    target = tmp_path / "out.fq"
    target.write_bytes(b"previous content\n")
    # ReadQualityError surfaces mid-iteration: quality length mismatch is a
    # format error at read time; use a FASTA input to force a read error.
    source = tmp_path / "in.fq"
    source.write_text(">a\nACGT\n")
    code, _, err = run(["filter-reads", str(source), "--output", str(target)])
    assert code == 2
    assert err.count("\n") == 1
    assert target.read_bytes() == b"previous content\n"
    leftovers = [p for p in tmp_path.iterdir() if p.name.startswith(".filter-reads-")]
    assert leftovers == []


def test_filter_reads_output_write_error_exit_1_no_partial(tmp_path) -> None:
    code, _, err = run(
        ["filter-reads", "-", "--output", str(tmp_path / "no-dir" / "out.fq")],
        "@r\n" + "A" * 30 + "\n+\n" + "I" * 30 + "\n",
    )
    assert code == 1
    assert err.count("\n") == 1
