"""Tests for the deduplicate-reads subcommand."""

from __future__ import annotations

import sys

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


def test_fasta_deduplicates_and_keeps_earliest() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"],
        ">r1 first\nACGT\n>r2\nTTTT\n>r3 second\nACGT\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">r1 first\nACGT\n>r2\nTTTT\n"


def test_fasta_output_uses_default_line_width() -> None:
    sequence = "A" * 130
    code, out, err = run(["deduplicate-reads", "-"], f">r1\n{sequence}\n")
    assert code == 0
    assert err == ""
    lines = out.splitlines()
    assert lines[0] == ">r1"
    assert lines[1] == "A" * 60
    assert lines[2] == "A" * 60
    assert lines[3] == "A" * 10
    assert out.endswith("\n") and not out.endswith("\n\n")


def test_fastq_keeps_highest_quality_sum() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"],
        "@r1\nACGT\n+r1\n!!!!\n@r2 best\nACGT\n+r2 best\nIIII\n",
    )
    assert code == 0
    assert err == ""
    assert out == "@r2 best\nACGT\n+r2 best\nIIII\n"


def test_fastq_output_is_four_lines_per_record() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"],
        "@r1\nACGT\n+r1\nIIII\n@r2\nTT\n+r2\n##\n",
    )
    assert code == 0
    assert err == ""
    assert out == "@r1\nACGT\n+r1\nIIII\n@r2\nTT\n+r2\n##\n"


def test_canonical_groups_reverse_complements() -> None:
    code, out, err = run(
        ["deduplicate-reads", "--canonical", "-"],
        ">r1\nAACGTA\n>r2\nTACGTT\n>r3\nGG\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">r1\nAACGTA\n>r3\nGG\n"


def test_explicit_input_format() -> None:
    code, out, err = run(
        ["deduplicate-reads", "--input-format", "fasta", "-"],
        ">r1\nACGT\n>r2\nACGT\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">r1\nACGT\n"


def test_explicit_format_mismatch_is_an_error() -> None:
    code, out, err = run(
        ["deduplicate-reads", "--input-format", "fastq", "-"],
        ">r1\nACGT\n",
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_result_ordered_by_first_occurrence() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"],
        ">r1\nGGGG\n>r2\nAAAA\n>r3\nGGGG\n>r4\nCCCC\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">r1\nGGGG\n>r2\nAAAA\n>r4\nCCCC\n"


def test_same_identifier_different_sequences_not_merged() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"],
        ">r1\nACGT\n>r1\nTGCA\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">r1\nACGT\n>r1\nTGCA\n"


def test_empty_input_is_a_format_error() -> None:
    code, out, err = run(["deduplicate-reads", "-"], "")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_invalid_base_is_an_error() -> None:
    code, out, err = run(["deduplicate-reads", "-"], ">r1\nACGX\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_fastq_quality_length_mismatch_is_an_error() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"], "@r1\nACGT\n+r1\n!!!\n"
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_unknown_argument_is_an_error() -> None:
    code, out, err = run(["deduplicate-reads", "-", "--bogus"], ">r1\nACGT\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_missing_input_file_returns_1(tmp_path) -> None:
    code, out, err = run(
        ["deduplicate-reads", str(tmp_path / "missing.fasta")]
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_file_input_and_output_round_trip(tmp_path) -> None:
    source = tmp_path / "in.fastq"
    source.write_text(
        "@r1\nACGT\n+r1\n!!!!\n@r2\nACGT\n+r2\nIIII\n", encoding="utf-8"
    )
    target = tmp_path / "out.fastq"
    code, out, err = run(
        ["deduplicate-reads", str(source), "--output", str(target)]
    )
    assert code == 0
    assert out == ""
    assert err == ""
    assert target.read_text(encoding="utf-8") == "@r2\nACGT\n+r2\nIIII\n"


def test_failed_run_preserves_existing_output(tmp_path) -> None:
    source = tmp_path / "in.fasta"
    source.write_text(">r1\nACGX\n", encoding="utf-8")
    target = tmp_path / "out.fasta"
    target.write_text("previous\n", encoding="utf-8")
    code, out, err = run(
        ["deduplicate-reads", str(source), "--output", str(target)]
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert target.read_text(encoding="utf-8") == "previous\n"
    assert list(tmp_path.iterdir()) == [source, target]


def test_unwritable_output_returns_1(tmp_path) -> None:
    source = tmp_path / "in.fasta"
    source.write_text(">r1\nACGT\n", encoding="utf-8")
    code, out, err = run(
        [
            "deduplicate-reads",
            str(source),
            "--output",
            str(tmp_path / "missing-dir" / "out.fasta"),
        ]
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_repeated_runs_are_byte_identical(tmp_path) -> None:
    source = tmp_path / "in.fasta"
    source.write_text(
        ">r1\nGGGG\n>r2\nAAAA\n>r3\nGGGG\n>r4\nCCCC\n", encoding="utf-8"
    )
    first = tmp_path / "one.fasta"
    second = tmp_path / "two.fasta"
    for target in (first, second):
        code, out, err = run(
            ["deduplicate-reads", str(source), "--output", str(target)]
        )
        assert code == 0
        assert err == ""
    assert first.read_bytes() == second.read_bytes()
