"""Tests for the map-reads subcommand."""

from __future__ import annotations

import json
import sys


def run(argv, stdin: str = ""):
    from io import StringIO

    old_in, old_out, old_err = sys.stdin, sys.stdout, sys.stderr
    sys.stdin = StringIO(stdin)
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    try:
        try:
            from genome_variant.cli import main

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


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


class TestMapReads:
    def test_basic_multi_reference_jsonl(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGT\n>c2\nTTTTACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n>b\nTTTT\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert len(lines) == 2
        payloads = [json.loads(line) for line in lines]
        for payload in payloads:
            assert list(payload) == EXPECTED_KEYS
        assert payloads[0] == {
            "record": 0,
            "id": "a",
            "mapped": True,
            "reference_record": 0,
            "reference": "c1",
            "reference_start": 0,
            "reference_end": 4,
            "query_start": 0,
            "query_end": 4,
            "strand": "+",
            "score": 8,
            "cigar": "4=",
        }
        assert payloads[1]["reference_record"] == 1
        assert payloads[1]["id"] == "b"

    def test_fastq_reads_auto_detected_and_quality_ignored(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\n!!!!\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 0
        payload = json.loads(out)
        assert payload["mapped"] is True
        assert payload["score"] == 8
        assert payload["cigar"] == "4="

    def test_explicit_reads_format(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.dat", "@a\nACGT\n+a\nIIII\n")
        code, out, err = run(
            ["map-reads", str(reference), str(reads), "--reads-format", "fastq"]
        )
        assert code == 0
        assert json.loads(out)["cigar"] == "4="

    def test_reads_from_stdin_while_reference_is_file(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nTTTT\n")
        code, out, err = run(
            ["map-reads", str(reference), "-"],
            ">a\nAAAA\n",
        )
        assert code == 0
        payload = json.loads(out)
        assert payload["strand"] == "-"
        assert payload["score"] == 8

    def test_reference_from_stdin_while_reads_are_file(self, tmp_path) -> None:
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["map-reads", "-", str(reads)], ">c1\nACGT\n")
        assert code == 0
        assert json.loads(out)["reference"] == "c1"

    def test_unmapped_record_uses_only_record_id_mapped_false(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nTTTT\n")
        reads = write(tmp_path, "q.fa", ">lost\nGCGC\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 0
        payload = json.loads(out)
        assert payload == {
            "record": 0,
            "id": "lost",
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

    def test_min_score_filters(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            ["map-reads", str(reference), str(reads), "--min-score", "9"]
        )
        assert code == 0
        assert json.loads(out)["mapped"] is False

    def test_scoring_options_accepted(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            [
                "map-reads",
                str(reference),
                str(reads),
                "--match",
                "1",
                "--mismatch",
                "1",
                "--gap-open",
                "0",
                "--gap-extend",
                "0",
                "--min-score",
                "1",
            ]
        )
        assert code == 0
        assert json.loads(out)["score"] == 4

    def test_output_ends_with_exactly_one_newline(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n>b\nGCGC\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 0
        assert out.endswith("\n")
        assert not out.endswith("\n\n")
        assert out.count("\n") == 2


class TestMapReadsErrors:
    def test_both_inputs_stdin_rejected(self) -> None:
        code, out, err = run(["map-reads", "-", "-"], ">a\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, tmp_path) -> None:
        empty = write(tmp_path, "empty.fa", "")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["map-reads", str(empty), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reads_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        empty = write(tmp_path, "empty.fa", "")
        code, out, err = run(["map-reads", str(reference), str(empty)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_min_score_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        for value in ("0", "-1", "x"):
            code, out, err = run(
                ["map-reads", str(reference), str(reads), "--min-score", value]
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1

    def test_reference_must_be_fasta(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fq", "@c1\nACGT\n+c1\nIIII\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_malformed_reads_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", "no header\nACGT\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_read_base_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nAZ\n")
        code, out, err = run(["map-reads", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert "'Z'" in err

    def test_missing_reference_file_exit_1(self, tmp_path) -> None:
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            ["map-reads", str(tmp_path / "missing.fa"), str(reads)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1


class TestMapReadsOutput:
    def test_file_output_is_atomic_and_byte_stable(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGT\n>c2\nTTTT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n>b\nAAAA\n>c\nGCGC\n")
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
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nAZ\n")
        target = tmp_path / "out.jsonl"
        target.write_text("KEEP\n")
        code, out, err = run(
            ["map-reads", str(reference), str(reads), "--output", str(target)]
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "KEEP\n"
        assert {p.name for p in tmp_path.iterdir()} == {
            "r.fa",
            "q.fa",
            "out.jsonl",
        }

    def test_unwritable_output_exit_1(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            [
                "map-reads",
                str(reference),
                str(reads),
                "--output",
                str(tmp_path / "no" / "out.jsonl"),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_other_subcommands_still_listed(self) -> None:
        code, out, err = run([])
        assert code == 0
        assert "map-reads" in out
        assert "align-pair" in out
