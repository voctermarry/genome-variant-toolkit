"""Tests for the ``compare-variants`` command line entry point."""

from __future__ import annotations

import json
import os
import sys
from io import StringIO

import pytest

from genome_variant.cli import main

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"

REF_FASTA = ">chr1\nACGGGGGT\n>chr2\nTTTTAAAA\n"


def run(argv, stdin: str = ""):
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


def vcf(sample: str, rows: list[str]) -> str:
    text = META + HEADER.format(s=sample)
    text += "".join(row + "\n" for row in rows)
    return text


def row(chrom="chr1", pos=2, ref="C", alt="G", flt="PASS", gt="0/1") -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\tGT\t{gt}"


@pytest.fixture
def workspace(tmp_path):
    (tmp_path / "ref.fa").write_text(REF_FASTA)
    return tmp_path


def write_vcf(directory, name, sample, rows):
    path = directory / name
    path.write_text(vcf(sample, rows))
    return path


class TestCompareVariantsCli:
    def test_basic_stdout(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row()])
        candidate = write_vcf(
            workspace, "c.vcf", "s2", [row(gt="1/1"), row(chrom="chr2", pos=5, ref="A", alt="G")]
        )
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        assert err == ""
        assert out.endswith("\n")
        payload = json.loads(out)
        assert payload == {
            "baseline_sample": "s1",
            "candidate_sample": "s2",
            "total": 2,
            "matched": 0,
            "genotype_mismatch": 1,
            "baseline_only": 0,
            "candidate_only": 1,
            "concordance": "0.000000",
            "differences": [
                {
                    "chrom": "chr1", "pos": 2, "ref": "C", "alt": "G",
                    "status": "genotype_mismatch",
                    "baseline_gt": "0/1", "candidate_gt": "1/1",
                },
                {
                    "chrom": "chr2", "pos": 5, "ref": "A", "alt": "G",
                    "status": "candidate_only",
                    "baseline_gt": None, "candidate_gt": "0/1",
                },
            ],
        }

    def test_output_file(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row()])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        output = workspace / "out.json"
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa"), "--output", str(output)]
        )
        assert code == 0
        assert out == ""
        assert err == ""
        payload = json.loads(output.read_text())
        assert payload["matched"] == 1
        assert payload["concordance"] == "1.000000"

    def test_output_file_is_atomic_on_data_error(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row(flt="q10")])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        output = workspace / "out.json"
        output.write_text("previous\n")
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa"), "--output", str(output)]
        )
        assert code == 2
        assert output.read_text() == "previous\n"
        assert not list(workspace.glob(".compare-variants-*"))

    def test_stdin_baseline(self, workspace) -> None:
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", "-", str(candidate),
             "--reference", str(workspace / "ref.fa")],
            stdin=vcf("s1", [row()]),
        )
        assert code == 0
        assert json.loads(out)["matched"] == 1

    def test_stdin_reference(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row()])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate), "--reference", "-"],
            stdin=REF_FASTA,
        )
        assert code == 0
        assert json.loads(out)["matched"] == 1

    def test_two_stdin_inputs_rejected(self, workspace) -> None:
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", "-", str(candidate), "--reference", "-"],
            stdin="",
        )
        assert code == 2
        assert "standard input" in err

    def test_missing_baseline_file(self, workspace) -> None:
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(workspace / "missing.vcf"), str(candidate),
             "--reference", str(workspace / "ref.fa")]
        )
        assert code == 1
        assert err != ""

    def test_missing_reference_file(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row()])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "missing.fa")]
        )
        assert code == 1
        assert err != ""

    def test_malformed_vcf_is_usage_error(self, workspace) -> None:
        baseline = workspace / "b.vcf"
        baseline.write_text("##fileformat=VCFv4.2\nchr1\t2\n")
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert err != ""

    def test_constraint_violation_is_usage_error(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row(gt="0/0")])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert "GT" in err

    def test_reference_mismatch_is_usage_error(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row(ref="A")])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert err != ""

    def test_unwritable_output_is_exit_1(self, workspace) -> None:
        baseline = write_vcf(workspace, "b.vcf", "s1", [row()])
        candidate = write_vcf(workspace, "c.vcf", "s2", [row()])
        code, out, err = run(
            ["compare-variants", str(baseline), str(candidate),
             "--reference", str(workspace / "ref.fa"),
             "--output", str(workspace / "no-such-dir" / "out.json")]
        )
        assert code == 1
        assert err != ""

    def test_missing_arguments_is_exit_2(self) -> None:
        code, out, err = run(["compare-variants"])
        assert code == 2
        assert err != ""
