"""Tests for the normalize-vcf subcommand."""

from __future__ import annotations

import sys

from genome_variant.cli import main

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
REFERENCE = ">chr1\nAAATTTAAA\n"


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


def make_reference(tmp_path) -> str:
    path = tmp_path / "ref.fa"
    path.write_text(REFERENCE)
    return str(path)


class TestNormalizeVcfCli:
    def test_stdin_to_stdout(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        vcf = "##fileformat=VCFv4.2\n" + HEADER + "\nchr1\t6\t.\tT\tTT\t.\t.\t.\n"
        code, out, err = run(["normalize-vcf", "-", "--reference", reference], vcf)
        assert code == 0
        assert err == ""
        assert out == (
            "##fileformat=VCFv4.2\n" + HEADER + "\nchr1\t3\t.\tA\tAT\t.\t.\t.\n"
        )

    def test_file_input_and_output_file(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(HEADER + "\nchr1\t2\t.\tAA\tAG\t.\t.\t.\n")
        out_path = tmp_path / "out.vcf"
        code, out, err = run(
            [
                "normalize-vcf",
                str(vcf_path),
                "--reference",
                reference,
                "--output",
                str(out_path),
            ]
        )
        assert code == 0
        assert out == ""
        assert err == ""
        assert out_path.read_text() == HEADER + "\nchr1\t3\t.\tA\tG\t.\t.\t.\n"

    def test_reference_may_come_from_stdin(self, tmp_path) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(HEADER + "\nchr1\t4\t.\tT\tA\t.\t.\t.\n")
        code, out, err = run(
            ["normalize-vcf", str(vcf_path), "--reference", "-"], REFERENCE
        )
        assert code == 0
        assert out == HEADER + "\nchr1\t4\t.\tT\tA\t.\t.\t.\n"

    def test_both_inputs_cannot_be_stdin(self) -> None:
        code, out, err = run(["normalize-vcf", "-", "--reference", "-"], "")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_reference_is_required(self, tmp_path) -> None:
        code, out, err = run(["normalize-vcf", "-"], "")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_vcf_format_error(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", reference],
            "chr1\t1\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":1:" in err and "column header" in err

    def test_reference_mismatch_error(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", reference],
            HEADER + "\nchr1\t1\t.\tC\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "chr1" in err and "'C'" in err and "'A'" in err

    def test_out_of_range_record_error(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", reference],
            HEADER + "\nchr1\t9\t.\tAA\tA\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":2:" in err and "chr1" in err

    def test_reference_sequence_format_error(self, tmp_path) -> None:
        bad_reference = tmp_path / "bad.fa"
        bad_reference.write_text(">chr1\nACGT@\n")
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(bad_reference)],
            HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_input_file(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        code, out, err = run(
            ["normalize-vcf", str(tmp_path / "nope.vcf"), "--reference", reference]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_file(self, tmp_path) -> None:
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", str(tmp_path / "nope.fa")],
            HEADER + "\n",
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_failed_run_preserves_existing_output(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        out_path = tmp_path / "out.vcf"
        out_path.write_text("previous contents\n")
        code, _, _ = run(
            [
                "normalize-vcf",
                "-",
                "--reference",
                reference,
                "--output",
                str(out_path),
            ],
            HEADER + "\nchr1\t1\t.\tC\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out_path.read_text() == "previous contents\n"
        assert list(tmp_path.glob(".normalize-vcf-*")) == []

    def test_empty_records_still_write_header(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        code, out, err = run(
            ["normalize-vcf", "-", "--reference", reference],
            "##fileformat=VCFv4.2\n" + HEADER + "\n",
        )
        assert code == 0
        assert err == ""
        assert out == "##fileformat=VCFv4.2\n" + HEADER + "\n"

    def test_output_is_byte_identical_across_runs(self, tmp_path) -> None:
        reference = make_reference(tmp_path)
        vcf = (
            "##fileformat=VCFv4.2\n"
            + HEADER
            + "\nchr1\t6\t.\tT\tTT\t.\t.\t.\nchr1\t2\t.\tAA\tAG,AC\t.\t.\t.\n"
        )
        first = run(["normalize-vcf", "-", "--reference", reference], vcf)
        second = run(["normalize-vcf", "-", "--reference", reference], vcf)
        assert first[0] == second[0] == 0
        assert first[1] == second[1]
        assert first[1].endswith("\n") and not first[1].endswith("\n\n")
