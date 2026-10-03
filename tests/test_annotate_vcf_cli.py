"""Tests for the ``annotate-vcf`` command line entry point."""

from __future__ import annotations

import sys
from io import StringIO

import pytest

from genome_variant.cli import main

HEADER = "##fileformat=VCFv4.2\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
# ATG AAA TTT GGG CCC : five codons on chr1.
REFERENCE = ">chr1\nATGAAATTTGGGCCC\n"
CDS = "chr1\ttest\tCDS\t1\t15\t.\t+\t0\tID=c1;Parent=t1\n"


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
    path = tmp_path / "features.gff3"
    path.write_text(CDS)
    return path


class TestAnnotateVcfCli:
    def test_stdin_vcf_to_stdout(self, reference_file, features_file) -> None:
        vcf = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            vcf,
        )
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert lines[0] == HEADER.rstrip("\n")
        assert lines[1].startswith("##INFO=<ID=GVANN,")
        assert lines[2] == COLUMNS.rstrip("\n")
        assert lines[3] == (
            "chr1\t4\t.\tA\tT\t.\t.\t"
            "GVANN=T|STOP_GAINED|HIGH|t1|4|AAA>TAA|K>*"
        )

    def test_file_inputs_and_output_roundtrip(
        self, tmp_path, reference_file, features_file
    ) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t2\t.\tT\tC\t.\t.\t.\n"
        )
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
        text = out_path.read_text()
        assert text.endswith(
            "chr1\t2\t.\tT\tC\t.\t.\t"
            "GVANN=C|START_LOST|HIGH|t1|2|ATG>ACG|M>T\n"
        )
        # Exactly one GVANN declaration and one trailing newline.
        assert text.count("##INFO=<ID=GVANN,") == 1
        assert text.endswith("\n") and not text.endswith("\n\n")

    def test_empty_record_vcf_writes_full_header(
        self, reference_file, features_file
    ) -> None:
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
        lines = out.splitlines()
        assert lines[0] == HEADER.rstrip("\n")
        assert lines[1].startswith("##INFO=<ID=GVANN,")
        assert lines[2] == COLUMNS.rstrip("\n")
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_only_one_input_may_be_stdin(self, tmp_path) -> None:
        ref = tmp_path / "ref.fa"
        ref.write_text(REFERENCE)
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", "-",
                "--features", str(tmp_path / "f.gff3"),
            ],
            HEADER + COLUMNS,
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_all_three_stdin_rejected(self) -> None:
        code, out, err = run(
            ["annotate-vcf", "-", "--reference", "-", "--features", "-"],
            HEADER + COLUMNS,
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_features_can_be_stdin(
        self, tmp_path, reference_file
    ) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", "-",
            ],
            CDS,
        )
        assert code == 0
        assert err == ""
        assert "STOP_GAINED" in out

    def test_reference_can_be_stdin(self, tmp_path, features_file) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", "-",
                "--features", str(features_file),
            ],
            REFERENCE,
        )
        assert code == 0
        assert err == ""
        assert "STOP_GAINED" in out

    def test_missing_vcf_file_exit_1(
        self, tmp_path, reference_file, features_file
    ) -> None:
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

    def test_missing_reference_file_exit_1(self, tmp_path, features_file) -> None:
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

    def test_missing_features_file_exit_1(self, tmp_path, reference_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(tmp_path / "missing.gff3"),
            ],
            HEADER + COLUMNS,
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_required_arguments(self) -> None:
        code, out, err = run(["annotate-vcf", "-"])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_malformed_vcf_exit_2_empty_stdout(
        self, reference_file, features_file
    ) -> None:
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
        assert err.endswith("\n")

    def test_malformed_gff_exit_2_empty_stdout(
        self, tmp_path, reference_file
    ) -> None:
        bad_features = tmp_path / "bad.gff3"
        bad_features.write_text("chr1\ttest\tCDS\t1\t9\t.\t+\t0\n")
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(bad_features),
            ],
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert str(bad_features) in err and ":1:" in err

    def test_contradictory_phase_exit_2(
        self, tmp_path, reference_file
    ) -> None:
        bad_features = tmp_path / "bad.gff3"
        bad_features.write_text(
            "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr1\ttest\tCDS\t10\t15\t.\t+\t1\tParent=t1\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(bad_features),
            ],
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert ":2:" in err

    def test_ref_mismatch_exit_2(self, reference_file, features_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS + "chr1\t3\t.\tC\tT\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert "chr1" in err and "POS 3" in err

    def test_missing_chrom_exit_2(self, reference_file, features_file) -> None:
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS + "chr9\t1\t.\tA\tG\t.\t.\t.\n",
        )
        assert code == 2
        assert out == ""
        assert "chr9" in err

    def test_existing_gvann_exit_2(self, reference_file, features_file) -> None:
        vcf = (
            HEADER
            + '##INFO=<ID=GVANN,Number=.,Type=String,Description="x">\n'
            + COLUMNS
            + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            vcf,
        )
        assert code == 2
        assert out == ""
        assert "GVANN" in err

    def test_error_preserves_existing_output(
        self, tmp_path, reference_file
    ) -> None:
        features = tmp_path / "f.gff3"
        features.write_text(CDS)
        target = tmp_path / "out.vcf"
        target.write_text("PREVIOUS\n")
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t99\t.\tA\tG\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", str(features),
                "--output", str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"

    def test_gff_error_preserves_existing_output(
        self, tmp_path, reference_file
    ) -> None:
        bad_features = tmp_path / "bad.gff3"
        bad_features.write_text("chr1\ttest\tCDS\t1\t9\t.\t+\tX\tParent=t1\n")
        target = tmp_path / "out.vcf"
        target.write_text("PREVIOUS\n")
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", str(vcf_path),
                "--reference", str(reference_file),
                "--features", str(bad_features),
                "--output", str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "PREVIOUS\n"

    def test_output_is_byte_stable(
        self, tmp_path, reference_file, features_file
    ) -> None:
        vcf_path = tmp_path / "in.vcf"
        vcf_path.write_text(
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT,G\t.\t.\tDP=1\n"
        )
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
        assert first.read_bytes() == second.read_bytes()
        data = first.read_bytes()
        assert data.endswith(b"\n") and not data.endswith(b"\n\n")

    def test_sample_columns_pass_through(
        self, reference_file, features_file
    ) -> None:
        header = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        )
        row = "chr1\t4\t.\tA\tT\t.\t.\tDP=9\tGT:DP\t0/1:9\t1/1:2\n"
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + header + row,
        )
        assert code == 0
        assert err == ""
        assert out.endswith(
            "chr1\t4\t.\tA\tT\t.\t.\tDP=9;"
            "GVANN=T|STOP_GAINED|HIGH|t1|4|AAA>TAA|K>*\tGT:DP\t0/1:9\t1/1:2\n"
        )

    def test_unsupported_and_non_coding_pass_through(
        self, reference_file, features_file
    ) -> None:
        rows = (
            # In-CDS SNV.
            "chr1\t4\t.\tA\tT\t.\t.\t.\n"
            # Symbolic ALT is unsupported but the record is retained.
            "chr1\t4\t.\tA\t<DEL>\t.\t.\t.\n"
        )
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            HEADER + COLUMNS + rows,
        )
        assert code == 0
        assert err == ""
        data_lines = [line for line in out.splitlines() if line.startswith("chr1")]
        assert data_lines[0].endswith(
            "GVANN=T|STOP_GAINED|HIGH|t1|4|AAA>TAA|K>*"
        )
        assert data_lines[1].endswith(
            "GVANN=<DEL>|UNSUPPORTED|MODIFIER|.|.|.|."
        )

    def test_normalized_frameshift_insertion(self, reference_file, features_file) -> None:
        # A normalized pure insertion: C inserted between coding
        # positions 3 and 4 of ATG|AAA|TTT|GGG|CCC.
        vcf = HEADER + COLUMNS + "chr1\t3\t.\tG\tGC\t.\t.\t.\n"
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            vcf,
        )
        assert code == 0
        assert err == ""
        assert out.splitlines()[3] == (
            "chr1\t3\t.\tG\tGC\t.\t.\t"
            "GVANN=GC|FRAMESHIFT|HIGH|t1|4|.|."
        )

    def test_inframe_deletion_is_moderate(self, reference_file, features_file) -> None:
        # In-frame deletion of codon AAA (coding positions 4-6).
        vcf = HEADER + COLUMNS + "chr1\t3\t.\tGAAA\tG\t.\t.\t.\n"
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features_file),
            ],
            vcf,
        )
        assert code == 0
        assert err == ""
        assert out.splitlines()[3] == (
            "chr1\t3\t.\tGAAA\tG\t.\t.\t"
            "GVANN=G|INFRAME_DELETION|MODERATE|t1|4|.|."
        )

    def test_cross_fragment_indel_exit_0_with_unsupported(
        self, tmp_path, reference_file
    ) -> None:
        # Two CDS fragments meeting at coding positions 9/10; a
        # boundary-spanning insertion is UNSUPPORTED, not an error.
        features = tmp_path / "split.gff3"
        features.write_text(
            "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr1\ttest\tCDS\t10\t15\t.\t+\t0\tParent=t1\n"
        )
        vcf = HEADER + COLUMNS + "chr1\t9\t.\tT\tTC\t.\t.\t.\n"
        code, out, err = run(
            [
                "annotate-vcf", "-",
                "--reference", str(reference_file),
                "--features", str(features),
            ],
            vcf,
        )
        assert code == 0
        assert err == ""
        assert out.splitlines()[3] == (
            "chr1\t9\t.\tT\tTC\t.\t.\t"
            "GVANN=TC|UNSUPPORTED|MODIFIER|t1|.|.|."
        )

