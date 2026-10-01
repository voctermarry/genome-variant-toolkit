"""Tests for the align-pair subcommand."""

from __future__ import annotations

import json
import sys

from genome_variant.cli import main

KEYS = [
    "reference",
    "query",
    "mode",
    "score",
    "reference_start",
    "reference_end",
    "query_start",
    "query_end",
    "cigar",
    "aligned_reference",
    "aligned_query",
]


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


def test_fasta_reference_fastq_query_from_files(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "query.fq"
    reference.write_text(">ref gene\nACGT\n")
    query.write_text("@q desc\nACGT\n+q desc\nIIII\n")
    code, out, err = run(["align-pair", str(reference), str(query)])
    assert code == 0
    assert err == ""
    assert out.endswith("\n") and not out.endswith("\n\n")
    result = json.loads(out)
    assert list(result) == KEYS
    assert result == {
        "reference": "ref",
        "query": "q",
        "mode": "global",
        "score": 8,
        "reference_start": 0,
        "reference_end": 4,
        "query_start": 0,
        "query_end": 4,
        "cigar": "4=",
        "aligned_reference": "ACGT",
        "aligned_query": "ACGT",
    }


def test_reference_from_stdin_and_query_file(tmp_path) -> None:
    query = tmp_path / "query.fa"
    query.write_text(">q\nAGT\n")
    code, out, err = run(
        ["align-pair", "-", str(query)], ">ref\nACGT\n"
    )
    assert code == 0
    assert err == ""
    result = json.loads(out)
    assert result["reference"] == "ref"
    assert result["query"] == "q"
    assert result["cigar"] == "1=1D2="


def test_query_from_stdin_and_reference_file(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">ref\nACGT\n")
    code, out, err = run(
        ["align-pair", str(reference), "-"], "@q\nACGT\n+q\nIIII\n"
    )
    assert code == 0
    result = json.loads(out)
    assert result["cigar"] == "4="


def test_local_mode_from_files(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "q.fa"
    reference.write_text(">r\nNNACGTNN\n")
    query.write_text(">q\nRRACGTRR\n")
    code, out, err = run(
        [
            "align-pair", str(reference), str(query),
            "--mode", "local",
            "--match-score", "3",
            "--mismatch-penalty", "1",
            "--gap-open-penalty", "2",
            "--gap-extend-penalty", "0",
        ]
    )
    assert code == 0
    result = json.loads(out)
    assert result["mode"] == "local"
    assert result["score"] == 12
    assert result["cigar"] == "4="
    assert result["reference_start"] == 2
    assert result["query_end"] == 6


def test_output_is_single_line_compact_json_with_one_newline(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">r\nACGT\n")
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")
    code, out, err = run(["align-pair", str(reference), str(query)])
    assert code == 0
    assert out.count("\n") == 1
    assert '": ' not in out
    assert out.startswith('{"reference":"r"')


def test_empty_input_returns_2_without_data(tmp_path) -> None:
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")
    code, out, err = run(["align-pair", str(empty), str(query)])
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    code, out, err = run(["align-pair", str(query), str(empty)])
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_multiple_records_return_2(tmp_path) -> None:
    two = tmp_path / "two.fa"
    query = tmp_path / "q.fa"
    two.write_text(">a\nAC\n>b\nGT\n")
    query.write_text(">q\nACGT\n")
    code, out, err = run(["align-pair", str(two), str(query)])
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert "more than one" in err
    code, out, err = run(["align-pair", str(query), str(two)])
    assert code == 2
    assert out == ""
    assert "more than one" in err


def test_both_stdin_returns_2() -> None:
    code, out, err = run(["align-pair", "-", "-"], ">r\nACGT\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def test_format_and_sequence_errors_return_2(tmp_path) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">q\nACGT\n")

    malformed = tmp_path / "malformed.fa"
    malformed.write_text("ACGT\n")
    code, out, err = run(["align-pair", str(malformed), str(query)])
    assert code == 2 and out == "" and err.count("\n") == 1

    bad_base = tmp_path / "bad.fa"
    bad_base.write_text(">x\nACGZ\n")
    code, out, err = run(["align-pair", str(bad_base), str(query)])
    assert code == 2 and out == "" and err.count("\n") == 1
    assert "'Z'" in err and "'x'" in err

    code, out, err = run(
        ["align-pair", "-", str(query), "--reference-format", "fastq"],
        ">r\nACGT\n",
    )
    assert code == 2 and out == "" and err.count("\n") == 1


def test_missing_file_returns_1(tmp_path) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">q\nA\n")
    code, out, err = run(["align-pair", str(tmp_path / "missing.fa"), str(query)])
    assert code == 1 and out == "" and err.count("\n") == 1
    code, out, err = run(["align-pair", str(query), str(tmp_path / "missing.fa")])
    assert code == 1 and out == "" and err.count("\n") == 1


def test_invalid_scoring_arguments_return_2(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "q.fa"
    reference.write_text(">r\nACGT\n")
    query.write_text(">q\nACGT\n")
    base = ["align-pair", str(reference), str(query)]
    for argv in (
        base + ["--match-score", "0"],
        base + ["--match-score", "-2"],
        base + ["--mismatch-penalty", "-1"],
        base + ["--gap-open-penalty", "x"],
        base + ["--gap-extend-penalty", "-3"],
        base + ["--mode", "semi"],
    ):
        code, out, err = run(argv)
        assert code == 2, argv
        assert out == ""
        assert err.count("\n") == 1


def test_output_file_is_written_atomically(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "q.fa"
    reference.write_text(">r\nACGT\n")
    query.write_text(">q\nACGT\n")
    target = tmp_path / "out.json"
    code, out, err = run(
        ["align-pair", str(reference), str(query), "--output", str(target)]
    )
    assert code == 0
    assert out == "" and err == ""
    data = target.read_bytes()
    assert data.endswith(b"\n") and not data.endswith(b"\n\n")
    assert json.loads(data)["cigar"] == "4="
    assert not list(tmp_path.glob("*.tmp"))


def test_failed_run_keeps_existing_output(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "q.fa"
    reference.write_text(">r\nACGZ\n")
    query.write_text(">q\nACGT\n")
    target = tmp_path / "out.json"
    target.write_text("ORIGINAL\n")
    code, out, err = run(
        ["align-pair", str(reference), str(query), "--output", str(target)]
    )
    assert code == 2
    assert out == ""
    assert target.read_text() == "ORIGINAL\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "out.json",
        "q.fa",
        "ref.fa",
    ]
    assert not list(tmp_path.glob("*.tmp"))


def test_unwritable_output_returns_1(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "q.fa"
    reference.write_text(">r\nACGT\n")
    query.write_text(">q\nACGT\n")
    code, out, err = run(
        [
            "align-pair",
            str(reference),
            str(query),
            "--output",
            str(tmp_path / "no" / "out.json"),
        ]
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_output_is_byte_stable(tmp_path) -> None:
    reference = tmp_path / "ref.fa"
    query = tmp_path / "q.fa"
    reference.write_text(">r\nacgt\n")
    query.write_text(">q\nACGT\n")
    first = tmp_path / "o1.json"
    second = tmp_path / "o2.json"
    code1, _, err1 = run(
        ["align-pair", str(reference), str(query), "--output", str(first)]
    )
    code2, _, err2 = run(
        ["align-pair", str(reference), str(query), "--output", str(second)]
    )
    assert code1 == code2 == 0
    assert err1 == err2 == ""
    assert first.read_bytes() == second.read_bytes()
