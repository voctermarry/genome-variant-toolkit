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


def q(scores) -> str:
    return "".join(chr(s + 33) for s in scores)


def test_fasta_duplicates_collapse_to_earliest() -> None:
    text = ">a first\nACGT\n>b second\nACGT\n>c\nTTTT\n"
    code, out, err = run(["deduplicate-reads", "-"], text)
    assert code == 0
    assert err == ""
    assert out == ">a first\nACGT\n>c\nTTTT\n"


def test_fasta_description_and_wrapping_use_defaults() -> None:
    long_seq = "A" * 61
    text = f">a gene\n{long_seq}\n>b\n{long_seq}\n"
    code, out, err = run(["deduplicate-reads", "-"], text)
    assert code == 0
    assert out == f">a gene\n{'A' * 60}\nA\n"


def test_fastq_keeps_highest_quality_sum() -> None:
    low = f"@low\nACGT\n+low\n{q([1, 1, 1, 1])}\n"
    high = f"@high\nACGT\n+high\n{q([40, 40, 40, 40])}\n"
    other = f"@other\nTTTT\n+other\n{q([5, 5, 5, 5])}\n"
    code, out, err = run(["deduplicate-reads", "-"], low + high + other)
    assert code == 0
    assert err == ""
    assert out == high + other


def test_fastq_tie_keeps_earliest() -> None:
    first = f"@first\nACGT\n+first\n{q([10, 20, 30, 40])}\n"
    second = f"@second\nACGT\n+second\n{q([40, 30, 20, 10])}\n"
    code, out, err = run(["deduplicate-reads", "-"], first + second)
    assert code == 0
    assert out == first


def test_fastq_keeps_four_line_layout() -> None:
    text = "@a desc\nACGT\n+a desc\n" + q([30, 30, 30, 30]) + "\n"
    code, out, err = run(["deduplicate-reads", "-"], text + text)
    assert code == 0
    lines = out.splitlines()
    assert len(lines) == 4
    assert lines[0] == "@a desc"


def test_equal_identifier_different_sequence_is_not_merged() -> None:
    text = ">x\nACGT\n>x\nACGA\n"
    code, out, err = run(["deduplicate-reads", "-"], text)
    assert code == 0
    assert out == text


def test_output_order_tracks_first_group_appearance() -> None:
    text = ">a\nGG\n>b\nAA\n>c\nGG\n>d\nCC\n>e\nAA\n"
    code, out, err = run(["deduplicate-reads", "-"], text)
    assert code == 0
    assert [line for line in out.splitlines() if line.startswith(">")] == [
        ">a",
        ">b",
        ">d",
    ]


def test_non_canonical_keeps_reverse_complements() -> None:
    text = ">a\nAAAA\n>b\nTTTT\n"
    code, out, err = run(["deduplicate-reads", "-"], text)
    assert code == 0
    assert out == text


def test_canonical_groups_reverse_complements_fasta() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-", "--canonical"], ">a\nAAAA\n>b\nTTTT\n"
    )
    assert code == 0
    assert err == ""
    assert out == ">a\nAAAA\n"


def test_canonical_fastq_winner_across_orientation_keeps_original() -> None:
    fwd = f"@fwd\nAAAA\n+fwd\n{q([1, 1, 1, 1])}\n"
    rev = f"@rev\nTTTT\n+rev\n{q([9, 9, 9, 9])}\n"
    code, out, err = run(
        ["deduplicate-reads", "-", "--canonical"], fwd + rev
    )
    assert code == 0
    assert out == rev


def test_canonical_supports_full_iupac() -> None:
    # rc("R") = "Y", canonical key is "R"; the R record is the earlier one.
    code, out, err = run(
        ["deduplicate-reads", "-", "--canonical"], ">a\nR\n>b\nY\n"
    )
    assert code == 0
    assert out == ">a\nR\n"


def test_explicit_input_format_fasta() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-", "--input-format", "fasta"],
        ">a\nACGT\n>b\nACGT\n",
    )
    assert code == 0
    assert out == ">a\nACGT\n"


def test_explicit_wrong_format_is_a_usage_error() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-", "--input-format", "fastq"], ">a\nACGT\n"
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_empty_input_is_a_format_error() -> None:
    code, out, err = run(["deduplicate-reads", "-"], "")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_empty_input_with_explicit_format_is_a_format_error() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-", "--input-format", "fasta"], ""
    )
    assert code == 2
    assert out == ""


def test_invalid_sequence_is_a_single_line_error_without_data() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"], ">a\nACGT\n>b\nACGZ\n>c\nTTTT\n"
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert "'Z'" in err and "'b'" in err


def test_malformed_fastq_is_a_usage_error() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"], "@a\nACGT\n+\nIII\n"
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_quality_character_outside_ascii_range_is_a_usage_error() -> None:
    code, out, err = run(
        ["deduplicate-reads", "-"], "@a\nACGT\n+a\n" + chr(127) * 4 + "\n"
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_unknown_argument_is_a_usage_error() -> None:
    code, out, err = run(["deduplicate-reads", "-", "--bogus"], ">a\nA\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_missing_input_file_returns_one() -> None:
    code, out, err = run(["deduplicate-reads", "/nonexistent/input.fa"])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_stdout_uses_unix_newlines_and_single_trailing_newline() -> None:
    text = ">a\nACGT\n>b\nACGT\n"
    code, out, err = run(["deduplicate-reads", "-"], text)
    assert code == 0
    assert out.endswith("\n") and not out.endswith("\n\n")
    assert "\r" not in out


def test_result_is_byte_stable_across_runs(tmp_path) -> None:
    source = tmp_path / "in.fa"
    source.write_text(">a\nAAAA\n>b\nTTTT\n>c\nAAAA\n")
    first = tmp_path / "o1.fa"
    second = tmp_path / "o2.fa"
    code1, _, err1 = run(
        ["deduplicate-reads", str(source), "--canonical", "--output", str(first)]
    )
    code2, _, err2 = run(
        ["deduplicate-reads", str(source), "--canonical", "--output", str(second)]
    )
    assert code1 == code2 == 0
    assert err1 == err2 == ""
    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes() == b">a\nAAAA\n"


def test_output_file_is_written_atomically(tmp_path) -> None:
    target = tmp_path / "dedup.fa"
    code, out, err = run(
        ["deduplicate-reads", "-", "--output", str(target)],
        ">a\nACGT\n>b\nACGT\n",
    )
    assert code == 0
    assert out == ""
    assert err == ""
    assert target.read_text(encoding="utf-8") == ">a\nACGT\n"
    assert not list(tmp_path.glob("*.tmp"))


def test_failed_run_keeps_existing_output_and_leaves_no_partial_file(tmp_path) -> None:
    target = tmp_path / "dedup.fa"
    target.write_text("ORIGINAL\n", encoding="utf-8")
    code, out, err = run(
        ["deduplicate-reads", "-", "--output", str(target)],
        ">a\nACGT\n>b\nACGZ\n",
    )
    assert code == 2
    assert out == ""
    assert target.read_text(encoding="utf-8") == "ORIGINAL\n"
    assert [p.name for p in tmp_path.iterdir()] == ["dedup.fa"]


def test_error_before_output_keeps_existing_target_even_for_first_record(tmp_path) -> None:
    target = tmp_path / "dedup.fq"
    target.write_text("KEEP\n")
    code, out, _ = run(
        ["deduplicate-reads", "-", "--output", str(target)],
        "@a\nACGZ\n+a\n" + q([40, 40, 40, 40]) + "\n",
    )
    assert code == 2
    assert out == ""
    assert target.read_text() == "KEEP\n"
    assert [p.name for p in tmp_path.iterdir()] == ["dedup.fq"]


def test_unwritable_output_path_returns_one(tmp_path) -> None:
    code, out, err = run(
        ["deduplicate-reads", "-", "--output", str(tmp_path / "no" / "out.fa")],
        ">a\nACGT\n",
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_fastq_file_roundtrip(tmp_path) -> None:
    source = tmp_path / "in.fq"
    text = "@r desc\nACGT\n+r desc\n" + q([40, 40, 40, 40]) + "\n"
    source.write_text(text)
    target = tmp_path / "out.fq"
    code, _, err = run(
        ["deduplicate-reads", str(source), "--output", str(target)]
    )
    assert code == 0
    assert err == ""
    assert target.read_bytes() == source.read_bytes()
