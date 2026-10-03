"""Tests for the ``summarize-variants`` command line entry point."""

from __future__ import annotations

import json
import sys
from io import StringIO

import pytest

from genome_variant.cli import main

HEADER = "##fileformat=VCFv4.2\n"
REFERENCE = ">chr1\nACGGTACGGT\n>chr2\nCCCAAAGGG\n"


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


def vcf_text(sample: str, rows: tuple[str, ...] = ()) -> str:
    columns = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    return HEADER + columns + "".join(rows)


def write_vcf(directory, name: str, sample: str, rows: tuple[str, ...] = ()):
    path = directory / name
    path.write_text(vcf_text(sample, rows))
    return path


@pytest.fixture
def reference_file(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(REFERENCE)
    return path


@pytest.fixture
def two_sample_case(tmp_path):
    """A manifest with two samples sharing one variant."""
    write_vcf(
        tmp_path,
        "s1.vcf",
        "s1",
        (
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:10:6,4\n",
            "chr1\t5\t.\tT\tA\t.\t.\t.\tGT:DP:AD\t1/1:8:0,8\n",
        ),
    )
    write_vcf(
        tmp_path,
        "s2.vcf",
        "s2",
        (
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t1/1:4:0,4\n",
            "chr2\t3\t.\tC\tA\t.\tPASS\t.\tGT:DP:AD\t0/1:3:2,1\n",
        ),
    )
    manifest = tmp_path / "manifest.jsonl"
    manifest.write_text(
        '{"sample":"s1","vcf":"s1.vcf"}\n{"sample":"s2","vcf":"s2.vcf"}\n'
    )
    return manifest


class TestSummarizeVariantsCli:
    def test_merge_two_samples_in_manifest_order(
        self, tmp_path, reference_file, two_sample_case
    ) -> None:
        code, out, err = run(
            [
                "summarize-variants",
                str(two_sample_case),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert len(lines) == 3
        first, second, third = (json.loads(line) for line in lines)
        assert first == {
            "chrom": "chr1",
            "pos": 2,
            "ref": "C",
            "alt": "T",
            "sample_count": 2,
            "allele_count": 3,
            "depth": 14,
            "samples": [
                {"sample": "s1", "gt": "0/1", "dp": 10, "ad": [6, 4]},
                {"sample": "s2", "gt": "1/1", "dp": 4, "ad": [0, 4]},
            ],
        }
        assert second == {
            "chrom": "chr1",
            "pos": 5,
            "ref": "T",
            "alt": "A",
            "sample_count": 1,
            "allele_count": 2,
            "depth": 8,
            "samples": [{"sample": "s1", "gt": "1/1", "dp": 8, "ad": [0, 8]}],
        }
        assert third == {
            "chrom": "chr2",
            "pos": 3,
            "ref": "C",
            "alt": "A",
            "sample_count": 1,
            "allele_count": 1,
            "depth": 3,
            "samples": [{"sample": "s2", "gt": "0/1", "dp": 3, "ad": [2, 1]}],
        }
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_records_are_normalized_before_merging(
        self, tmp_path, reference_file
    ) -> None:
        # Same variant expressed minimized (s1) and unminimized (s2).
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            ("chr1\t2\t.\tCG\tC\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        write_vcf(
            tmp_path, "s2.vcf", "s2",
            ("chr1\t4\t.\tGT\tT\t.\tPASS\t.\tGT:DP:AD\t1/1:7:0,7\n",),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text(
            '{"sample":"s1","vcf":"s1.vcf"}\n{"sample":"s2","vcf":"s2.vcf"}\n'
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert err == ""
        (only,) = (json.loads(line) for line in out.splitlines())
        assert only["chrom"] == "chr1"
        assert (only["pos"], only["ref"], only["alt"]) == (2, "CG", "C")
        assert only["sample_count"] == 2
        assert [s["sample"] for s in only["samples"]] == ["s1", "s2"]

    def test_results_follow_reference_record_order(
        self, tmp_path, reference_file
    ) -> None:
        # chr2 is listed before chr1 in the VCF but after it in the reference.
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            (
                "chr2\t3\t.\tC\tA\t.\tPASS\t.\tGT:DP:AD\t0/1:3:2,1\n",
                "chr1\t9\t.\tG\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:3:2,1\n",
                "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:3:2,1\n",
            ),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert err == ""
        keys = [
            (row["chrom"], row["pos"]) for row in map(json.loads, out.splitlines())
        ]
        assert keys == [("chr1", 2), ("chr1", 9), ("chr2", 3)]

    def test_empty_manifest_produces_empty_output(
        self, tmp_path, reference_file
    ) -> None:
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text("\n   \n")
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_manifest_from_stdin_resolves_against_cwd(
        self, tmp_path, reference_file, monkeypatch
    ) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            ("chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        monkeypatch.chdir(tmp_path)
        code, out, err = run(
            ["summarize-variants", "-", "--reference", str(reference_file)],
            '{"sample":"s1","vcf":"s1.vcf"}\n',
        )
        assert code == 0
        assert err == ""
        assert json.loads(out)["samples"][0]["sample"] == "s1"

    def test_manifest_and_reference_cannot_both_be_stdin(self) -> None:
        code, out, err = run(["summarize-variants", "-", "--reference", "-"], "")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_reference_can_be_stdin(self, tmp_path) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            ("chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", "-"], REFERENCE
        )
        assert code == 0
        assert err == ""
        assert json.loads(out)["chrom"] == "chr1"

    def test_non_ascii_sample_name_preserved(
        self, tmp_path, reference_file
    ) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "样本一",
            ("chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"样本一","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert err == ""
        assert "样本一" in out
        assert "\\u" not in out

    @pytest.mark.parametrize(
        "line",
        (
            "not json",
            '["s1", "s1.vcf"]',
            '{"sample":"s1"}',
            '{"sample":"s1","vcf":"s1.vcf","extra":1}',
            '{"sample":1,"vcf":"s1.vcf"}',
            '{"sample":"","vcf":"s1.vcf"}',
            '{"sample":"a\tb","vcf":"s1.vcf"}',
        ),
    )
    def test_bad_manifest_line_exit_2_with_line_number(
        self, tmp_path, reference_file, line
    ) -> None:
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text("\n" + line + "\n")
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":2:" in err

    def test_duplicate_sample_exit_2(self, tmp_path, reference_file) -> None:
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text(
            '{"sample":"s1","vcf":"a.vcf"}\n{"sample":"s1","vcf":"b.vcf"}\n'
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":2:" in err and "s1" in err

    def test_missing_manifest_exit_1(self, tmp_path, reference_file) -> None:
        code, out, err = run(
            [
                "summarize-variants",
                str(tmp_path / "missing.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_exit_1(self, tmp_path) -> None:
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text("")
        code, out, err = run(
            [
                "summarize-variants",
                str(manifest),
                "--reference",
                str(tmp_path / "missing.fa"),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_vcf_exit_1(self, tmp_path, reference_file) -> None:
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"missing.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_sample_name_mismatch_exit_2(self, tmp_path, reference_file) -> None:
        write_vcf(tmp_path, "s1.vcf", "other")
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "s1.vcf" in err and "s1" in err

    @pytest.mark.parametrize(
        "row",
        (
            # two ALT alleles
            "chr1\t2\t.\tC\tT,G\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",
            # symbolic ALT
            "chr1\t2\t.\tC\t<DEL>\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",
            # FILTER not PASS or .
            "chr1\t2\t.\tC\tT\t.\tq10\t.\tGT:DP:AD\t0/1:5:3,2\n",
            # FORMAT without AD
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP\t0/1:5\n",
            # GT not 0/1 or 1/1
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/0:5:5,0\n",
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t1/2:5:3,2\n",
            # DP not a non-negative integer
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:-5:3,2\n",
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:x:3,2\n",
            # AD without exactly two values
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3\n",
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2,0\n",
            # AD does not sum to DP
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,3\n",
        ),
    )
    def test_bad_record_exit_2_with_vcf_source_and_line(
        self, tmp_path, reference_file, row
    ) -> None:
        write_vcf(tmp_path, "s1.vcf", "s1", (row,))
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "s1.vcf:3:" in err

    def test_duplicate_normalized_key_within_sample_exit_2(
        self, tmp_path, reference_file
    ) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            (
                "chr1\t2\t.\tCG\tC\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",
                "chr1\t4\t.\tGT\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",
            ),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "s1.vcf:4:" in err

    def test_ref_mismatch_exit_2(self, tmp_path, reference_file) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            ("chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "chr1" in err and "POS 2" in err

    def test_missing_chrom_exit_2(self, tmp_path, reference_file) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            ("chr9\t1\t.\tA\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert "chr9" in err

    def test_malformed_vcf_exit_2(self, tmp_path, reference_file) -> None:
        (tmp_path / "s1.vcf").write_text(HEADER + "chr1\t2\t.\tC\tT\n")
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_output_file_roundtrip_and_byte_stability(
        self, tmp_path, reference_file, two_sample_case
    ) -> None:
        first, second = tmp_path / "o1.jsonl", tmp_path / "o2.jsonl"
        argv = ["summarize-variants", str(two_sample_case), "--reference",
                str(reference_file)]
        code1, out1, err1 = run(argv + ["--output", str(first)])
        code2, out2, err2 = run(argv + ["--output", str(second)])
        assert code1 == code2 == 0
        assert out1 == out2 == err1 == err2 == ""
        data = first.read_bytes()
        assert data == second.read_bytes()
        assert data.endswith(b"\n") and not data.endswith(b"\n\n")
        for line in data.decode("utf-8").splitlines():
            json.loads(line)

    def test_failure_preserves_existing_output(
        self, tmp_path, reference_file
    ) -> None:
        write_vcf(
            tmp_path, "s1.vcf", "s1",
            ("chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:5:3,2\n",),
        )
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s1.vcf"}\n')
        target = tmp_path / "out.jsonl"
        target.write_text("PREVIOUS\n")
        code, out, err = run(
            [
                "summarize-variants", str(manifest),
                "--reference", str(reference_file),
                "--output", str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "PREVIOUS\n"

    def test_reference_required(self, tmp_path) -> None:
        manifest = tmp_path / "manifest.jsonl"
        manifest.write_text("")
        code, out, err = run(["summarize-variants", str(manifest)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
