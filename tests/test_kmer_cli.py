"""Tests for the kmer-index subcommand."""

from __future__ import annotations

import json
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


def parse_lines(out: str):
    return [json.loads(line) for line in out.splitlines()]


def test_basic_fasta_index() -> None:
    code, out, err = run(
        ["kmer-index", "-", "--k", "3"],
        ">r1 a description\nACGTAC\n>r2\nACG\n",
    )
    assert code == 0
    assert err == ""
    rows = parse_lines(out)
    assert [row["kmer"] for row in rows] == ["ACG", "GTA"]
    assert rows[0] == {
        "kmer": "ACG",
        "count": 3,
        "occurrences": [
            {"record": 0, "id": "r1", "position": 0, "strand": "+"},
            {"record": 0, "id": "r1", "position": 1, "strand": "-"},
            {"record": 1, "id": "r2", "position": 0, "strand": "+"},
        ],
    }
    assert rows[1]["count"] == 2


def test_output_is_compact_json_with_one_trailing_newline() -> None:
    code, out, err = run(["kmer-index", "-", "--k", "2"], ">r1\nACGT\n")
    assert code == 0
    assert err == ""
    assert out.endswith("\n") and not out.endswith("\n\n")
    for line in out.splitlines():
        assert '": "' not in line  # compact separators
        assert line.startswith('{"kmer":"')


def test_no_canonical_indexes_forward_windows_only() -> None:
    code, out, err = run(["kmer-index", "-", "--k", "3", "--no-canonical"], ">r1\nCGT\n")
    assert code == 0
    assert parse_lines(out) == [
        {
            "kmer": "CGT",
            "count": 1,
            "occurrences": [{"record": 0, "id": "r1", "position": 0, "strand": "+"}],
        }
    ]


def test_fastq_input_is_indexed_without_quality_filtering() -> None:
    fastq = "@r1\nACGT\n+\n!!!!\n@r2\nACGT\n+\nIIII\n"
    code, out, err = run(["kmer-index", "-", "--k", "4"], fastq)
    assert code == 0
    assert err == ""
    rows = parse_lines(out)
    assert rows[0]["count"] == 2  # low-quality read r1 is kept


def test_input_format_option() -> None:
    code, out, err = run(
        ["kmer-index", "-", "--k", "2", "--input-format", "fastq"],
        "@r1\nACGT\n+\nIIII\n",
    )
    assert code == 0
    assert err == ""
    code, out, err = run(
        ["kmer-index", "-", "--k", "2", "--input-format", "fasta"],
        "@r1\nACGT\n+\nIIII\n",
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_non_ascii_identifier_is_preserved() -> None:
    code, out, err = run(["kmer-index", "-", "--k", "4"], ">café\nACGT\n")
    assert code == 0
    assert "café" in out
    assert "\\u" not in out


def test_empty_index_produces_empty_output() -> None:
    code, out, err = run(["kmer-index", "-", "--k", "2"], ">r1\nNNN\n")
    assert code == 0
    assert out == ""
    assert err == ""


def test_missing_k_is_a_usage_error() -> None:
    code, out, err = run(["kmer-index", "-"], ">r1\nACGT\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_invalid_k_is_a_usage_error() -> None:
    for value in ("0", "-3", "abc"):
        code, out, err = run(["kmer-index", "-", "--k", value], ">r1\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1


def test_invalid_sequence_is_a_single_line_error_without_data() -> None:
    code, out, err = run(["kmer-index", "-", "--k", "2"], ">x\nACGZ\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert "'Z'" in err and "'x'" in err


def test_missing_input_file_returns_one() -> None:
    code, out, err = run(["kmer-index", "/nonexistent/input.fa", "--k", "3"])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_output_file_is_written_atomically(tmp_path) -> None:
    target = tmp_path / "index.jsonl"
    code, out, err = run(
        ["kmer-index", "-", "--k", "2", "--output", str(target)], ">r1\nACGT\n"
    )
    assert code == 0
    assert out == ""
    assert err == ""
    rows = parse_lines(target.read_text(encoding="utf-8"))
    assert [row["kmer"] for row in rows] == ["AC", "CG"]
    assert not list(tmp_path.glob("*.tmp"))


def test_failed_run_keeps_existing_output_and_leaves_no_partial_file(tmp_path) -> None:
    target = tmp_path / "index.jsonl"
    target.write_text("ORIGINAL\n", encoding="utf-8")
    code, out, err = run(
        ["kmer-index", "-", "--k", "2", "--output", str(target)], ">x\nACGZ\n"
    )
    assert code == 2
    assert out == ""
    assert target.read_text(encoding="utf-8") == "ORIGINAL\n"
    assert [p.name for p in tmp_path.iterdir()] == ["index.jsonl"]


def test_empty_index_writes_empty_file(tmp_path) -> None:
    target = tmp_path / "index.jsonl"
    code, out, err = run(
        ["kmer-index", "-", "--k", "2", "--output", str(target)], ">r1\nNNN\n"
    )
    assert code == 0
    assert target.read_text(encoding="utf-8") == ""


def test_unwritable_output_path_returns_one(tmp_path) -> None:
    code, out, err = run(
        ["kmer-index", "-", "--k", "2", "--output", str(tmp_path / "no" / "out.jsonl")],
        ">r1\nACGT\n",
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_existing_subcommands_still_work() -> None:
    code, out, err = run(["version"])
    assert code == 0
    code, out, err = run([])
    assert code == 0
    assert "kmer-index" in out
