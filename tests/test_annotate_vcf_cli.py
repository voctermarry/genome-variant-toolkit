"""Tests for the ``annotate-vcf`` command line entry point."""

from __future__ import annotations

import sys
from io import StringIO

import pytest

from genome_variant.cli import main

HEADER = "##fileformat=VCFv4.2\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

# ATG AAA TTT GGG CCC TAA* TTA GGG
REFERENCE = ">chr1\nATGAAATTTGGGCCCTAATTAGGG\n"
PLUS_CDS = "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx1\n"
GVANN_META = '##INFO=<ID=GVANN,Number=.,Type=String,Description="Genomic Variant Annotation: ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|AA_CHANGE">\n'


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
    path.write_text(REFERENCE)
    return path


@pytest.fixture
def features_file(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text("##gff-version 3\n" + PLUS_CDS)
    return path


def vcf(row: str) -> str:
    return HEADER + COLUMNS + row


class TestAnnotateVcfCli:
    def test_stdin_vcf_to_stdout(self, reference_file, features_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            vcf("chr1\t1\t.\tA\tT\t.\t.\t.\n"),
        )
        assert code == 0
        assert err == ""
        assert out == (
            HEADER + GVANN_META + COLUMNS
            + "chr1\t1\t.\tA\tT\t.\t.\tGVANN=T|START_LOST|HIGH|tx1|1|ATG>TTG|M1L\n"
        )

    def test_file_inputs_and_output_roundtrip(self, tmp_path, reference_file, features_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(vcf("chr1\t7\t.\tT\tC\t.\t.\t.\n"))
        out_path = tmp_path / "out.vcf"
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", str(features_file),
                "--output", str(out_path),
            ]
        )
        assert code == 0
        assert out == ""
        assert err == ""
        assert out_path.read_bytes() == (
            (
                HEADER + GVANN_META + COLUMNS
                + "chr1\t7\t.\tT\tC\t.\t.\tGVANN=C|MISSENSE|MODERATE|tx1|7|TTT>CTT|F3L\n"
            ).encode()
        )

    def test_empty_record_vcf_writes_full_header(self, reference_file, features_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS,
        )
        assert code == 0
        assert err == ""
        assert out == HEADER + GVANN_META + COLUMNS

    @pytest.mark.parametrize(
        "argv_extra",
        [
            ["-", "--reference", "-", "--features", "genes.gff3"],
            ["-", "--reference", "ref.fa", "--features", "-"],
            ["in.vcf", "--reference", "-", "--features", "-"],
        ],
    )
    def test_at_most_one_stdin_input(self, argv_extra, tmp_path) -> None:
        argv = ["annotate-vcf"]
        for token in argv_extra:
            argv.append(str(tmp_path / token) if token.endswith((".gff3", ".fa", ".vcf")) else token)
        code, out, err = run(argv, HEADER + COLUMNS)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_reference_stdin_vcf_file(self, tmp_path, features_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(vcf("chr1\t1\t.\tA\tT\t.\t.\t.\n"))
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", "-",
                "--features", str(features_file),
            ],
            REFERENCE,
        )
        assert code == 0
        assert "START_LOST" in out

    def test_features_stdin_vcf_file(self, tmp_path, reference_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(vcf("chr1\t1\t.\tA\tT\t.\t.\t.\n"))
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", "-",
            ],
            PLUS_CDS,
        )
        assert code == 0
        assert "START_LOST" in out

    def test_missing_vcf_file_exit_1(self, reference_file, features_file, tmp_path) -> None:
        code, out, err = run(
            [
                "annotate-vcf", str(tmp_path / "missing.vcf"),
                "--reference", str(reference_file),
                "--features", str(features_file),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_file_exit_1(self, features_file, tmp_path) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(tmp_path / "missing.fa"),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS,
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_features_file_exit_1(self, reference_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", "/no/such/dir/genes.gff3",
            ],
            HEADER + COLUMNS,
        )
        assert code == 1
        assert out == ""

    def test_malformed_vcf_exit_2_empty_stdout(self, reference_file, features_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS + "chr1\t0\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_ref_mismatch_exit_2(self, reference_file, features_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            vcf("chr1\t1\t.\tC\tT\t.\t.\t.\n"),
        )
        assert code == 2
        assert out == ""
        assert "chr1" in err and "POS 1" in err

    def test_bad_gff_exit_2_with_line_number(self, tmp_path, reference_file) -> None:
        bad_features = tmp_path / "genes.gff3"
        bad_features.write_text(PLUS_CDS + "chr1\ttest\tCDS\t1\t18\t.\t+\t9\tID=c2;Parent=tx2\n")
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(bad_features),
            ],
            vcf("chr1\t1\t.\tA\tT\t.\t.\t.\n"),
        )
        assert code == 2
        assert out == ""
        assert ":2:" in err

    def test_gvann_already_declared_exit_2(self, reference_file, features_file) -> None:
        text = (
            HEADER
            + '##INFO=<ID=GVANN,Number=.,Type=String,Description="x">\n'
            + COLUMNS
            + "chr1\t1\t.\tA\tT\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            text,
        )
        assert code == 2
        assert out == ""
        assert "GVANN" in err

    def test_bad_reference_fasta_exit_2(self, tmp_path, features_file) -> None:
        bad_reference = tmp_path / "ref.fa"
        bad_reference.write_text(">chr1\nACGZ\n")
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(bad_reference),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS,
        )
        assert code == 2
        assert out == ""

    def test_error_preserves_existing_output(
        self, tmp_path, reference_file, features_file
    ) -> None:
        target = tmp_path / "out.vcf"
        target.write_text("PREVIOUS\n")
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(vcf("chr1\t99\t.\tA\tT\t.\t.\t.\n"))
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", str(features_file),
                "--output", str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "PREVIOUS\n"

    def test_output_is_byte_stable(self, tmp_path, reference_file, features_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(vcf("chr1\t1\t.\tA\tG,T\t.\t.\tDP=1\n"))
        first, second = tmp_path / "o1.vcf", tmp_path / "o2.vcf"
        code1, _, err1 = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", str(features_file),
                "--output", str(first),
            ]
        )
        code2, _, err2 = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", str(features_file),
                "--output", str(second),
            ]
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        data = first.read_bytes()
        assert data == second.read_bytes()
        assert data.endswith(b"\n")
        assert not data.endswith(b"\n\n")

    def test_sample_columns_pass_through(self, reference_file, features_file) -> None:
        columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        row = "chr1\t7\t.\tT\tC\t.\t.\tDP=3\tGT:DP\t0/1:9\t1/1:2\n"
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + columns + row,
        )
        assert code == 0
        assert out.endswith("\tGT:DP\t0/1:9\t1/1:2\n")

    def test_required_arguments(self) -> None:
        code, out, err = run(["annotate-vcf", "-"], HEADER + COLUMNS)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
