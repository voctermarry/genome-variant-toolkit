"""Tests for the map-reads subcommand."""

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
    "record",
    "id",
    "mapped",
    "reference_record",
    "reference",
    "reference_start",
    "reference_end",
    "query_start",
    "query_end",
    "strand",
    "score",
    "cigar",
]


class TestMapReadsBasic:
    def test_fasta_reads_from_stdin(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r1\nTTACGTAA\n>r2\nCCCC\n")
        code, out, err = run(
            ["map-reads", str(reference), "-"], ">q\nACGT\n"
        )
        assert code == 0
        assert err == ""
        assert out.endswith("\n") and not out.endswith("\n\n")
        lines = out.splitlines()
        assert len(lines) == 1
        payload = json.loads(lines[0])
        assert list(payload) == EXPECTED_KEYS
        assert payload == {
            "record": 0,
            "id": "q",
            "mapped": True,
            "reference_record": 0,
            "reference": "r1",
            "reference_start": 2,
            "reference_end": 6,
            "query_start": 0,
            "query_end": 4,
            "strand": "+",
            "score": 8,
            "cigar": "4=",
        }

    def test_fastq_reads_auto_detected(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nAAGTCCAA\n")
        reads = tmp_path / "reads.fq"
        reads.write_text("@q\nGGACAA\n+q\n!!!!!!\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        payload = json.loads(out)
        assert payload["mapped"] is True
        assert payload["strand"] == "-"
        assert (payload["query_start"], payload["query_end"]) == (0, 4)
        assert payload["score"] == 8

    def test_reference_from_stdin(self, tmp_path) -> None:
        reads = tmp_path / "reads.fa"
        reads.write_text(">q\nACGT\n")
        code, out, err = run(
            ["map-reads", "-", str(reads)], ">r\nTTACGTAA\n"
        )
        assert code == 0
        assert json.loads(out)["mapped"] is True

    def test_unmapped_read_still_emitted(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nCCCC\n")
        code, out, err = run(
            ["map-reads", str(reference), "-"], ">q\nAAAA\n"
        )
        assert code == 0
        assert err == ""
        payload = json.loads(out)
        assert list(payload) == EXPECTED_KEYS
        assert payload == {
            "record": 0,
            "id": "q",
            "mapped": False,
            "reference_record": None,
            "reference": None,
            "reference_start": None,
            "reference_end": None,
            "query_start": None,
            "query_end": None,
            "strand": None,
            "score": None,
            "cigar": None,
        }

    def test_one_line_per_read_in_input_order(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nCCGG\n")
        code, out, err = run(
            ["map-reads", str(reference), "-"],
            ">a\nAAAA\n>b\nCCGG\n>c\nTTTT\n",
        )
        assert code == 0
        lines = out.splitlines()
        assert len(lines) == 3
        payloads = [json.loads(line) for line in lines]
        assert [p["record"] for p in payloads] == [0, 1, 2]
        assert [p["id"] for p in payloads] == ["a", "b", "c"]
        assert [p["mapped"] for p in payloads] == [False, True, False]

    def test_min_score_option(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nAC\n")
        code, out, err = run(
            ["map-reads", str(reference), "-", "--min-score", "5"],
            ">q\nAC\n",
        )
        assert code == 0
        assert json.loads(out)["mapped"] is False
        code, out, err = run(
            ["map-reads", str(reference), "-", "--min-score", "4"],
            ">q\nAC\n",
        )
        assert code == 0
        assert json.loads(out)["mapped"] is True

    def test_scoring_options_match_align_pair(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(
            [
                "map-reads", str(reference), "-",
                "--match", "3",
                "--mismatch", "1",
                "--gap-open", "4",
                "--gap-extend", "1",
            ],
            ">q\nACGT\n",
        )
        assert code == 0
        assert json.loads(out)["score"] == 12

    def test_explicit_reads_format(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        reads = tmp_path / "reads.dat"
        reads.write_text("@q\nACGT\n+q\nIIII\n")
        code, out, err = run(
            ["map-reads", str(reference), str(reads), "--reads-format", "fastq"]
        )
        assert code == 0
        assert json.loads(out)["mapped"] is True

    def test_compact_separators_and_final_newline(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(["map-reads", str(reference), "-"], ">q\nACGT\n")
        assert code == 0
        assert '": ' not in out
        assert ", " not in out
        assert out.count("\n") == 1


class TestMapReadsErrors:
    def test_both_stdin_rejected(self) -> None:
        code, out, err = run(["map-reads", "-", "-"], ">r\nA\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text("")
        code, out, err = run(["map-reads", str(reference), "-"], ">q\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reads_exit_2(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(["map-reads", str(reference), "-"], "")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(["map-reads", str(reference), "-"], ">q\nAZ\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "'Z'" in err

    def test_malformed_reads_exit_2(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(
            ["map-reads", str(reference), "-"], "no header\nAC\n"
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_file_exit_1(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(
            ["map-reads", str(reference), str(tmp_path / "missing.fa")]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

        reads = tmp_path / "reads.fa"
        reads.write_text(">q\nACGT\n")
        code, out, err = run(
            ["map-reads", str(tmp_path / "missing.fa"), str(reads)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_option_values_exit_2(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        for option, value in [
            ("--min-score", "0"),
            ("--min-score", "-1"),
            ("--min-score", "x"),
            ("--match", "0"),
            ("--mismatch", "-2"),
            ("--gap-open", "-1"),
            ("--gap-extend", "x"),
        ]:
            code, out, err = run(
                ["map-reads", str(reference), "-", option, value], ">q\nAC\n"
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1


class TestMapReadsOutput:
    def test_file_output_is_atomic_and_byte_stable(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nAAGTCCAA\n")
        reads = tmp_path / "reads.fa"
        reads.write_text(">q1\nGGACAA\n>q2\nAAAA\n")
        first = tmp_path / "a.jsonl"
        second = tmp_path / "b.jsonl"
        code1, _, err1 = run(
            ["map-reads", str(reference), str(reads), "--output", str(first)]
        )
        code2, _, err2 = run(
            ["map-reads", str(reference), str(reads), "--output", str(second)]
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        assert first.read_bytes().endswith(b"\n")
        assert not list(tmp_path.glob("*.tmp"))

    def test_failed_run_preserves_existing_output(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        target = tmp_path / "out.jsonl"
        target.write_text("KEEP\n")
        code, out, err = run(
            ["map-reads", str(reference), "-", "--output", str(target)],
            ">q\nAZ\n",
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "KEEP\n"
        assert {p.name for p in tmp_path.iterdir()} == {"ref.fa", "out.jsonl"}

    def test_unwritable_output_exit_1(self, tmp_path) -> None:
        reference = tmp_path / "ref.fa"
        reference.write_text(">r\nACGT\n")
        code, out, err = run(
            [
                "map-reads", str(reference), "-",
                "--output", str(tmp_path / "no" / "out.jsonl"),
            ],
            ">q\nACGT\n",
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_other_subcommands_unchanged(self) -> None:
        code, out, err = run(["version"])
        assert code == 0 and out.strip() == "0.1.0"
        code, out, err = run([])
        assert code == 0 and "map-reads" in out
