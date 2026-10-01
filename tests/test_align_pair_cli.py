"""Tests for the align-pair subcommand."""

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


EXPECTED_KEYS = [
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


class TestAlignPairStdin:
    def test_basic_global(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nACGT\n")
        code, out, err = run(["align-pair", "-", str(query)], ">r\nACGT\n")
        assert code == 0
        assert err == ""
        assert out.endswith("\n") and not out.endswith("\n\n")
        payload = json.loads(out)
        assert list(payload) == EXPECTED_KEYS
        assert payload == {
            "reference": "r",
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

    def test_reference_can_be_file_and_query_stdin(self, tmp_path) -> None:
        reference = tmp_path / "r.fa"
        reference.write_text(">ref\nA\n")
        code, out, err = run(["align-pair", str(reference), "-"], ">qry\nA\n")
        assert code == 0
        assert json.loads(out)["cigar"] == "1="

    def test_local_mode(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nGGACGTGG\n")
        code, out, err = run(
            [
                "align-pair", "-", str(query),
                "--mode", "local",
                "--match", "2",
                "--mismatch", "3",
                "--gap-open", "5",
                "--gap-extend", "2",
            ],
            ">r\nCCACGTCC\n",
        )
        assert code == 0
        payload = json.loads(out)
        assert payload["mode"] == "local"
        assert payload["score"] == 8
        assert (
            payload["reference_start"],
            payload["reference_end"],
            payload["query_start"],
            payload["query_end"],
        ) == (2, 6, 2, 6)
        assert payload["cigar"] == "4="

    def test_local_no_positive_fragment(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nGT\n")
        code, out, err = run(
            ["align-pair", "-", str(query), "--mode", "local"],
            ">r\nAC\n",
        )
        assert code == 0
        payload = json.loads(out)
        assert payload["score"] == 0
        assert payload["cigar"] == ""
        assert payload["aligned_reference"] == ""
        assert payload["aligned_query"] == ""
        assert (
            payload["reference_start"],
            payload["reference_end"],
            payload["query_start"],
            payload["query_end"],
        ) == (0, 0, 0, 0)

    def test_fastq_quality_does_not_affect_scoring(self, tmp_path) -> None:
        query = tmp_path / "q.fq"
        query.write_text("@q\nACGT\n+q\n!!!!\n")
        code, out, err = run(["align-pair", "-", str(query)], ">r\nACGT\n")
        assert code == 0
        assert json.loads(out)["score"] == 8

    def test_explicit_input_formats(self, tmp_path) -> None:
        query = tmp_path / "q.dat"
        query.write_text("@q\nAC\n+q\nII\n")
        code, out, err = run(
            [
                "align-pair", "-", str(query),
                "--query-format", "fastq",
            ],
            ">r\nAC\n",
        )
        assert code == 0
        assert json.loads(out)["cigar"] == "2="

    def test_compact_separators_and_final_newline(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nA\n")
        code, out, err = run(["align-pair", "-", str(query)], ">r\nA\n")
        assert code == 0
        assert '": ' not in out
        assert ", " not in out
        assert out.count("\n") == 1


class TestAlignPairErrors:
    def test_both_stdin_rejected(self) -> None:
        code, out, err = run(["align-pair", "-", "-"], ">r\nA\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nA\n")
        code, out, err = run(["align-pair", "-", str(query)], "")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_multiple_records_rejected(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q1\nA\n>q2\nC\n")
        code, out, err = run(["align-pair", "-", str(query)], ">r\nA\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

        multi = tmp_path / "r.fa"
        multi.write_text(">r1\nA\n>r2\nC\n")
        one = tmp_path / "one.fa"
        one.write_text(">q\nA\n")
        code, out, err = run(["align-pair", str(multi), str(one)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nAZ\n")
        code, out, err = run(["align-pair", "-", str(query)], ">r\nAA\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "'Z'" in err

    def test_malformed_input_exit_2(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text("no header\nAC\n")
        code, out, err = run(["align-pair", "-", str(query)], ">r\nAC\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_file_exit_1(self, tmp_path) -> None:
        one = tmp_path / "one.fa"
        one.write_text(">r\nA\n")
        code, out, err = run(
            ["align-pair", str(one), str(tmp_path / "missing.fa")]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_scoring_values_exit_2(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nA\n")
        for option, value in [
            ("--match", "0"),
            ("--match", "-1"),
            ("--mismatch", "-2"),
            ("--gap-open", "-1"),
            ("--gap-extend", "x"),
            ("--mode", "weird"),
        ]:
            code, out, err = run(
                ["align-pair", "-", str(query), option, value], ">r\nA\n"
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1


class TestAlignPairOutput:
    def test_file_output_is_atomic_and_byte_stable(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nACGT\n")
        first = tmp_path / "a.json"
        second = tmp_path / "b.json"
        code1, _, err1 = run(
            ["align-pair", "-", str(query), "--output", str(first)],
            ">r\nACGT\n",
        )
        code2, _, err2 = run(
            ["align-pair", "-", str(query), "--output", str(second)],
            ">r\nACGT\n",
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        assert first.read_bytes().endswith(b"\n")
        assert not list(tmp_path.glob("*.tmp"))

    def test_failed_run_preserves_existing_output(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nAZ\n")
        target = tmp_path / "out.json"
        target.write_text("KEEP\n")
        code, out, err = run(
            ["align-pair", "-", str(query), "--output", str(target)],
            ">r\nAA\n",
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "KEEP\n"
        assert {p.name for p in tmp_path.iterdir()} == {"q.fa", "out.json"}

    def test_unwritable_output_exit_1(self, tmp_path) -> None:
        query = tmp_path / "q.fa"
        query.write_text(">q\nA\n")
        code, out, err = run(
            [
                "align-pair", "-", str(query),
                "--output", str(tmp_path / "no" / "out.json"),
            ],
            ">r\nA\n",
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_other_subcommands_unchanged(self) -> None:
        code, out, err = run(["version"])
        assert code == 0 and out.strip() == "0.1.0"
        code, out, err = run([])
        assert code == 0 and "align-pair" in out
