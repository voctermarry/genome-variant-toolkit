"""Tests for the ``normalize-vcf`` command line entry point."""

from __future__ import annotations

import sys
from io import StringIO

import pytest

from genome_variant.cli import main

HEADER = "##fileformat=VCFv4.2\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


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


@pytest.fixture
def reference_file(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(">chr1\nACGGT\n")
    return path


class TestNormalizeVcfCli:
    def test_normalize_stdin_vcf_to_stdout(self, reference_file) -> None:
        vcf = HEADER + COLUMNS + "chr1\t4\t.\tGT\tT\t.\t.\t.\n"
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(reference_file)], vcf
        )
        assert code == 0
        assert err == ""
        assert out == (
            HEADER
            + COLUMNS
            + "chr1\t2\t.\tCG\tC\t.\t.\t.\n"
        )

    def test_file_inputs_and_output_roundtrip(self, tmp_path, reference_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t2\t.\tC\tT\t.\t.\t.\n"
        )
        out_path = tmp_path / "out.vcf"
        code, out, err = run(
            [
                "normalize-vcf",
                str(vcf_path),
                "--reference",
                str(reference_file),
                "--output",
                str(out_path),
            ]
        )
        assert code == 0
        assert out == ""
        assert err == ""
        assert out_path.read_bytes() == (
            (HEADER + COLUMNS + "chr1\t2\t.\tC\tT\t.\t.\t.\n").encode()
        )

    def test_empty_record_vcf_writes_full_header(self, reference_file) -> None:
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(reference_file)],
            HEADER + COLUMNS,
        )
        assert code == 0
        assert err == ""
        assert out == HEADER + COLUMNS

    def test_vcf_and_reference_cannot_both_be_stdin(self) -> None:
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", "-"],
            HEADER + COLUMNS,
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_reference_can_be_stdin_when_vcf_is_file(self, tmp_path) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t2\t.\tC\tT\t.\t.\t.\n"
        )
        code, out, err = run(
            ["normalize-vcf", str(vcf_path), "--reference", "-"],
            ">chr1\nACGGT\n",
        )
        assert code == 0
        assert err == ""
        assert "chr1\t2\t.\tC\tT" in out

    def test_missing_vcf_file_exit_1(self, reference_file, tmp_path) -> None:
        code, out, err = run(
            [
                "normalize-vcf",
                str(tmp_path / "missing.vcf"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_file_exit_1(self, tmp_path) -> None:
        code, out, err = run(
            [
                "normalize-vcf",
                "-",
                "--reference",
                str(tmp_path / "missing.fa"),
            ],
            HEADER + COLUMNS,
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_malformed_vcf_exit_2_single_line_empty_stdout(self, reference_file) -> None:
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(reference_file)],
            HEADER + COLUMNS + "chr1\t0\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert err.endswith("\n")

    def test_bad_reference_fasta_exit_2(self, tmp_path) -> None:
        bad_reference = tmp_path / "ref.fa"
        bad_reference.write_text(">chr1\nACGZ\n")
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(bad_reference)],
            HEADER + COLUMNS + "chr1\t1\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_ref_mismatch_exit_2(self, reference_file) -> None:
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(reference_file)],
            HEADER + COLUMNS + "chr1\t3\t.\tC\tT\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "chr1" in err and "POS 3" in err

    def test_missing_chrom_exit_2(self, reference_file) -> None:
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(reference_file)],
            HEADER + COLUMNS + "chr9\t1\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert "chr9" in err

    def test_error_after_some_records_preserves_existing_output(
        self, tmp_path, reference_file
    ) -> None:
        target = tmp_path / "out.vcf"
        target.write_text("PREVIOUS\n")
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER
            + COLUMNS
            + "chr1\t2\t.\tC\tT\t.\t.\t.\n"
            + "chr1\t2\t.\tX\tY\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "normalize-vcf",
                str(vcf_path),
                "--reference",
                str(reference_file),
                "--output",
                str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"

    def test_output_is_byte_stable(self, tmp_path, reference_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t4\t.\tGT\tT\t.\t.\t.\n"
        )
        first, second = tmp_path / "o1.vcf", tmp_path / "o2.vcf"
        code1, _, err1 = run(
            [
                "normalize-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--output", str(first),
            ]
        )
        code2, _, err2 = run(
            [
                "normalize-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--output", str(second),
            ]
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        data = first.read_bytes()
        assert data.endswith(b"\n")
        assert not data.endswith(b"\n\n")

    def test_reference_required(self) -> None:
        code, out, err = run(["normalize-vcf", "-"], HEADER + COLUMNS)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_sample_columns_pass_through(self, reference_file) -> None:
        header = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        )
        row = "chr1\t4\t.\tGT\tT\t.\t.\t.\tGT:DP\t0/1:9\t1/1:2\n"
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(reference_file)],
            HEADER + header + row,
        )
        assert code == 0
        assert err == ""
        assert out.endswith(
            "chr1\t2\t.\tCG\tC\t.\t.\t.\tGT:DP\t0/1:9\t1/1:2\n"
        )
