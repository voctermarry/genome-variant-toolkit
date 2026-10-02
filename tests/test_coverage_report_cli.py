"""Tests for the coverage-report subcommand."""

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
    "reference_record",
    "reference",
    "length",
    "mapped_reads",
    "observed_bases",
    "covered_bases",
    "coverage_fraction",
    "mean_depth",
    "concordant_bases",
    "concordance_fraction",
]


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


class TestCoverageReport:
    def test_basic_multi_reference_jsonl(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGT\n>c2\nTTTTACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n>b\nTTTT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert len(lines) == 2
        payloads = [json.loads(line) for line in lines]
        for payload in payloads:
            assert list(payload) == EXPECTED_KEYS
        assert payloads[0] == {
            "reference_record": 0,
            "reference": "c1",
            "length": 8,
            "mapped_reads": 1,
            "observed_bases": 4,
            "covered_bases": 4,
            "coverage_fraction": "0.500000",
            "mean_depth": "0.500000",
            "concordant_bases": 4,
            "concordance_fraction": "1.000000",
        }
        assert payloads[1]["reference_record"] == 1
        assert payloads[1]["reference"] == "c2"
        assert payloads[1]["mapped_reads"] == 1

    def test_reference_without_hits_still_reported(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n>c2\nGGGG\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        payloads = [json.loads(line) for line in out.splitlines()]
        assert payloads[1] == {
            "reference_record": 1,
            "reference": "c2",
            "length": 4,
            "mapped_reads": 0,
            "observed_bases": 0,
            "covered_bases": 0,
            "coverage_fraction": "0.000000",
            "mean_depth": "0.000000",
            "concordant_bases": 0,
            "concordance_fraction": "0.000000",
        }

    def test_duplicate_identifiers_not_merged(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">dup\nACGT\n>dup\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        payloads = [json.loads(line) for line in out.splitlines()]
        assert [p["reference_record"] for p in payloads] == [0, 1]
        assert [p["reference"] for p in payloads] == ["dup", "dup"]
        assert [p["mapped_reads"] for p in payloads] == [1, 0]

    def test_fastq_reads_auto_detected_and_quality_ignored(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\n!!!!\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        payload = json.loads(out)
        assert payload["mapped_reads"] == 1
        assert payload["observed_bases"] == 4
        assert payload["concordance_fraction"] == "1.000000"

    def test_explicit_reads_format(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.dat", "@a\nACGT\n+a\nIIII\n")
        code, out, err = run(
            ["coverage-report", str(reference), str(reads), "--reads-format", "fastq"]
        )
        assert code == 0
        assert json.loads(out)["mapped_reads"] == 1

    def test_reads_from_stdin_while_reference_is_file(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nTTTT\n")
        code, out, err = run(
            ["coverage-report", str(reference), "-"],
            ">a\nAAAA\n",
        )
        assert code == 0
        payload = json.loads(out)
        assert payload["mapped_reads"] == 1
        assert payload["covered_bases"] == 4

    def test_reference_from_stdin_while_reads_are_file(self, tmp_path) -> None:
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["coverage-report", "-", str(reads)], ">c1\nACGT\n")
        assert code == 0
        assert json.loads(out)["reference"] == "c1"

    def test_mismatch_statistics(self, tmp_path) -> None:
        # 4=1X3=: the mismatched column adds depth but no concordance.
        reference = write(tmp_path, "r.fa", ">c1\nACGTTCGT\n")
        reads = write(tmp_path, "q.fa", ">x\nACGTACGT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        payload = json.loads(out)
        assert payload["observed_bases"] == 8
        assert payload["covered_bases"] == 8
        assert payload["concordant_bases"] == 7
        assert payload["concordance_fraction"] == "0.875000"

    def test_deletion_statistics(self, tmp_path) -> None:
        # 4=1D3=: the deleted reference position keeps zero depth.
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGT\n")
        reads = write(tmp_path, "q.fa", ">y\nACGTCGT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        payload = json.loads(out)
        assert payload["observed_bases"] == 7
        assert payload["covered_bases"] == 7
        assert payload["coverage_fraction"] == "0.875000"
        assert payload["mean_depth"] == "0.875000"

    def test_min_score_filters(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            ["coverage-report", str(reference), str(reads), "--min-score", "9"]
        )
        assert code == 0
        assert json.loads(out)["mapped_reads"] == 0

    def test_scoring_options_accepted(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            [
                "coverage-report",
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
        assert json.loads(out)["mapped_reads"] == 1

    def test_output_ends_with_exactly_one_newline(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n>c2\nTTTT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 0
        assert out.endswith("\n")
        assert not out.endswith("\n\n")
        assert out.count("\n") == 2

    def test_output_file_written_atomically(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        output = tmp_path / "out.jsonl"
        code, out, err = run(
            ["coverage-report", str(reference), str(reads), "--output", str(output)]
        )
        assert code == 0
        assert out == ""
        text = output.read_text()
        assert text.endswith("\n") and not text.endswith("\n\n")
        assert json.loads(text)["mapped_reads"] == 1
        # No temporary files remain.
        assert [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")] == []

    def test_failed_run_preserves_existing_output(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT-\n")
        output = write(tmp_path, "out.jsonl", "previous\n")
        code, out, err = run(
            ["coverage-report", str(reference), str(reads), "--output", str(output)]
        )
        assert code == 2
        assert out == ""
        assert output.read_text() == "previous\n"
        assert [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")] == []


class TestCoverageReportErrors:
    def test_both_inputs_stdin_rejected(self) -> None:
        code, out, err = run(["coverage-report", "-", "-"], ">a\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, tmp_path) -> None:
        empty = write(tmp_path, "empty.fa", "")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["coverage-report", str(empty), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reads_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        empty = write(tmp_path, "empty.fa", "")
        code, out, err = run(["coverage-report", str(reference), str(empty)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_min_score_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        for value in ("0", "-1", "x"):
            code, out, err = run(
                ["coverage-report", str(reference), str(reads), "--min-score", value]
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1

    def test_reference_must_be_fasta(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fq", "@c1\nACGT\n+c1\nIIII\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_read_symbol_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT-\n")
        code, out, err = run(["coverage-report", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_file_exit_1(self, tmp_path) -> None:
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        code, out, err = run(
            ["coverage-report", str(tmp_path / "missing.fa"), str(reads)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reads_file_exit_1(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        code, out, err = run(
            ["coverage-report", str(reference), str(tmp_path / "missing.fa")]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_unwritable_output_exit_1(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fa", ">a\nACGT\n")
        missing_dir = tmp_path / "nope" / "out.jsonl"
        code, out, err = run(
            ["coverage-report", str(reference), str(reads), "--output", str(missing_dir)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1
