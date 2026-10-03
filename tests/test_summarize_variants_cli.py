"""Tests for the ``summarize-variants`` command line entry point."""

from __future__ import annotations

import json
import sys
from io import StringIO

import pytest

from genome_variant.cli import main

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"

REF_FASTA = ">chr1\nACGGT\n>chr2\nTTTTT\n"


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


def vcf(sample: str, rows: list[str], meta: str = META) -> str:
    lines = meta + HEADER.format(s=sample)
    lines += "".join(row + "\n" for row in rows)
    return lines


def row(chrom="chr1", pos=2, ref="C", alt="G", flt="PASS", gt="0/1",
        dp=10, ad=(6, 4), sample_format="GT:DP:AD") -> str:
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\t{sample_format}\t"
        f"{gt}:{dp}:{ad[0]},{ad[1]}"
    )


@pytest.fixture
def reference_file(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(REF_FASTA)
    return path


@pytest.fixture
def workspace(tmp_path, reference_file):
    vcfs = tmp_path / "vcfs"
    vcfs.mkdir()
    return tmp_path


def write_vcf(path, sample, rows, meta=META):
    path.write_text(vcf(sample, rows, meta))


class TestSummarizeVariantsCli:
    def test_single_sample_basic_merge(self, workspace, reference_file) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row(gt="0/1", dp=10, ad=(6, 4))])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"vcfs/s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert err == ""
        data = [json.loads(line) for line in out.splitlines()]
        assert data == [
            {
                "chrom": "chr1",
                "pos": 2,
                "ref": "C",
                "alt": "G",
                "sample_count": 1,
                "allele_count": 1,
                "depth": 10,
                "samples": [
                    {"sample": "s1", "gt": "0/1", "dp": 10, "ad": [6, 4]}
                ],
            }
        ]
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_merges_equivalent_calls_across_samples(
        self, workspace, reference_file
    ) -> None:
        # s1: deletion written at pos4 normalizes to chr1:2 CG>C.
        write_vcf(
            workspace / "vcfs" / "s1.vcf",
            "s1",
            [row(pos=4, ref="GT", alt="T", gt="1/1", dp=8, ad=(0, 8))],
        )
        # s2: the same haplotype already minimal at pos2.
        write_vcf(
            workspace / "vcfs" / "s2.vcf",
            "s2",
            [row(pos=2, ref="CG", alt="C", gt="0/1", dp=10, ad=(6, 4))],
        )
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
            '{"sample":"s2","vcf":"vcfs/s2.vcf"}\n'
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert err == ""
        lines = [json.loads(line) for line in out.splitlines()]
        assert len(lines) == 1
        merged = lines[0]
        assert (merged["chrom"], merged["pos"], merged["ref"], merged["alt"]) == (
            "chr1", 2, "CG", "C"
        )
        assert merged["sample_count"] == 2
        assert merged["allele_count"] == 3  # 1/1 (2) + 0/1 (1)
        assert merged["depth"] == 18
        assert [s["sample"] for s in merged["samples"]] == ["s1", "s2"]

    def test_samples_keep_manifest_order_not_vcf_order(
        self, workspace, reference_file
    ) -> None:
        # Manifest order b, a; output samples must be b then a.
        write_vcf(workspace / "vcfs" / "a.vcf", "a", [row()])
        write_vcf(workspace / "vcfs" / "b.vcf", "b", [row(gt="1/1", dp=5, ad=(0, 5))])
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            '{"sample":"b","vcf":"vcfs/b.vcf"}\n'
            '{"sample":"a","vcf":"vcfs/a.vcf"}\n'
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        merged = json.loads(out)
        assert [s["sample"] for s in merged["samples"]] == ["b", "a"]
        assert merged["allele_count"] == 3
        assert merged["depth"] == 15

    def test_output_sorted_by_reference_then_pos_ref_alt(
        self, workspace, reference_file
    ) -> None:
        write_vcf(
            workspace / "vcfs" / "s1.vcf",
            "s1",
            [
                row(chrom="chr2", pos=1, ref="T", alt="A"),
                row(chrom="chr1", pos=3, ref="G", alt="A"),
                row(chrom="chr1", pos=1, ref="A", alt="T"),
            ],
        )
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"vcfs/s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        keys = [(d["chrom"], d["pos"], d["ref"], d["alt"]) for d in
                (json.loads(line) for line in out.splitlines())]
        assert keys == [
            ("chr1", 1, "A", "T"),
            ("chr1", 3, "G", "A"),
            ("chr2", 1, "T", "A"),
        ]

    def test_only_actual_callers_listed(self, workspace, reference_file) -> None:
        # s1 carries two variants, s2 carries only one: each summary lists
        # only the samples calling it.
        write_vcf(
            workspace / "vcfs" / "s1.vcf",
            "s1",
            [row(pos=1, ref="A", alt="T"), row(pos=2, ref="C", alt="G")],
        )
        write_vcf(workspace / "vcfs" / "s2.vcf", "s2", [row(pos=2, ref="C", alt="G")])
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
            '{"sample":"s2","vcf":"vcfs/s2.vcf"}\n'
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        lines = [json.loads(line) for line in out.splitlines()]
        by_pos = {d["pos"]: d for d in lines}
        assert [s["sample"] for s in by_pos[1]["samples"]] == ["s1"]
        assert by_pos[1]["sample_count"] == 1
        assert [s["sample"] for s in by_pos[2]["samples"]] == ["s1", "s2"]
        assert by_pos[2]["sample_count"] == 2

    def test_blank_manifest_lines_ignored(self, workspace, reference_file) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row()])
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            "\n"
            '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
            "\n"
            "   \n"
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert out.count("\n") == 1

    def test_empty_manifest_produces_empty_output(
        self, workspace, reference_file
    ) -> None:
        manifest = workspace / "m.jsonl"
        manifest.write_text("\n\n")
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_stdin_manifest_resolves_vs_cwd(self, workspace, reference_file, monkeypatch) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row()])
        monkeypatch.chdir(workspace)
        code, out, err = run(
            ["summarize-variants", "-", "--reference", str(reference_file)],
            '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n',
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "chr1"

    def test_file_manifest_resolves_vs_its_directory(
        self, workspace, reference_file, monkeypatch
    ) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row()])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"vcfs/s1.vcf"}\n')
        # Running from an unrelated cwd must still find the VCF.
        monkeypatch.chdir("/")
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "chr1"

    def test_reference_can_be_stdin_with_file_manifest(
        self, workspace, reference_file
    ) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row()])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"vcfs/s1.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", "-"],
            REF_FASTA,
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "chr1"

    def test_manifest_and_reference_cannot_both_be_stdin(self) -> None:
        code, out, err = run(
            ["summarize-variants", "-", "--reference", "-"],
            '{"sample":"s1","vcf":"x"}\n',
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_non_ascii_sample_preserved_verbatim(
        self, workspace, reference_file
    ) -> None:
        write_vcf(workspace / "vcfs" / "u.vcf", "样", [row()])
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            '{"sample":"样","vcf":"vcfs/u.vcf"}\n', encoding="utf-8"
        )
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert "样" in out
        assert "\\u" not in out

    def test_output_file_roundtrip_and_single_newline(
        self, workspace, reference_file
    ) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row()])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"vcfs/s1.vcf"}\n')
        target = workspace / "out.jsonl"
        code, out, err = run(
            [
                "summarize-variants", str(manifest),
                "--reference", str(reference_file),
                "--output", str(target),
            ]
        )
        assert code == 0
        assert out == "" and err == ""
        data = target.read_bytes()
        assert data.endswith(b"\n") and not data.endswith(b"\n\n")

    def test_empty_result_writes_empty_file(self, workspace, reference_file) -> None:
        manifest = workspace / "m.jsonl"
        manifest.write_text("")
        target = workspace / "out.jsonl"
        target.write_text("OLD\n")
        code, _, err = run(
            [
                "summarize-variants", str(manifest),
                "--reference", str(reference_file),
                "--output", str(target),
            ]
        )
        assert code == 0
        assert err == ""
        assert target.read_bytes() == b""

    def test_byte_stable_across_runs(self, workspace, reference_file) -> None:
        write_vcf(workspace / "vcfs" / "s1.vcf", "s1", [row()])
        write_vcf(
            workspace / "vcfs" / "s2.vcf",
            "s2",
            [row(pos=4, ref="GT", alt="T", gt="1/1", dp=8, ad=(0, 8))],
        )
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
            '{"sample":"s2","vcf":"vcfs/s2.vcf"}\n'
        )
        argv = ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        _, out1, err1 = run(argv)
        _, out2, err2 = run(argv)
        assert err1 == err2 == ""
        assert out1 == out2

    # -- manifest validation ------------------------------------------------

    def test_bad_json_reports_line_number(self, reference_file) -> None:
        code, out, err = run(
            ["summarize-variants", "-", "--reference", str(reference_file)],
            '{"sample":"a","vcf":"x"}\n{not json}\n',
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":2:" in err

    def test_line_that_is_not_object(self, reference_file) -> None:
        code, _, err = run(
            ["summarize-variants", "-", "--reference", str(reference_file)],
            '["sample","vcf"]\n',
        )
        assert code == 2
        assert ":1:" in err

    @pytest.mark.parametrize(
        "text",
        [
            '{"sample":"a"}\n',                       # missing vcf
            '{"vcf":"x"}\n',                          # missing sample
            '{"sample":"a","vcf":"x","extra":1}\n',   # extra field
            '{"sample":"","vcf":"x"}\n',             # empty sample
            '{"sample":"a\tb","vcf":"x"}\n',         # tab in sample
            '{"sample":3,"vcf":"x"}\n',              # non-string sample
            '{"sample":"a","vcf":""}\n',             # empty vcf
            '{"sample":"a","vcf":4}\n',              # non-string vcf
        ],
    )
    def test_bad_manifest_fields_exit_2(self, text, reference_file) -> None:
        code, out, err = run(
            ["summarize-variants", "-", "--reference", str(reference_file)],
            text,
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":1:" in err

    def test_duplicate_sample_reports_second_line(
        self, workspace, reference_file
    ) -> None:
        code, _, err = run(
            ["summarize-variants", "-", "--reference", str(reference_file)],
            '{"sample":"a","vcf":"x"}\n{"sample":"a","vcf":"y"}\n',
        )
        assert code == 2
        assert ":2:" in err
        assert "duplicate" in err

    # -- VCF structural / genotype validation --------------------------------

    def test_two_samples_in_vcf_rejected(self, workspace, reference_file) -> None:
        text = (
            META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:3:1,2\t0/1:3:1,2\n"
        )
        target = workspace / "s.vcf"
        target.write_text(text)
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert ":1:" in err and "one sample" in err

    def test_vcf_sample_name_must_match_manifest(
        self, workspace, reference_file
    ) -> None:
        write_vcf(workspace / "s.vcf", "other", [row()])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert "does not match" in err

    @pytest.mark.parametrize(
        "alt",
        ["G,T", "<DEL>", ".", "*"],
    )
    def test_single_plain_alt_required(
        self, alt, workspace, reference_file
    ) -> None:
        write_vcf(workspace / "s.vcf", "s1", [row(alt=alt)])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""

    def test_filter_must_be_pass_or_dot(self, workspace, reference_file) -> None:
        write_vcf(workspace / "s.vcf", "s1", [row(flt="LowQual")])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert "FILTER" in err

    @pytest.mark.parametrize(
        "gt,dp,ad",
        [
            ("0/0", "10", "6,4"),
            ("1/2", "10", "0,10"),
            ("0|1", "10", "6,4"),
            ("0/1", "10", "6"),       # one AD value
            ("0/1", "10", "6,2,2"),   # three AD values
            ("0/1", "10", "6,5"),     # AD sum != DP
            ("0/1", "-1", "0,-1"),    # negative
            ("0/1", "10x", "6,4"),    # not an integer
            ("0/1", "10", "a,4"),     # non-numeric AD
        ],
    )
    def test_bad_genotype_fields_exit_2(
        self, gt, dp, ad, workspace, reference_file
    ) -> None:
        text = (
            META
            + HEADER.format(s="s1")
            + f"chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT:DP:AD\t{gt}:{dp}:{ad}\n"
        )
        target = workspace / "s.vcf"
        target.write_text(text)
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_format_missing_required_field(self, workspace, reference_file) -> None:
        text = (
            META
            + HEADER.format(s="s1")
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT:DP\t0/1:10\n"
        )
        target = workspace / "s.vcf"
        target.write_text(text)
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert "AD" in err

    def test_format_field_order_and_extra_fields_allowed(
        self, workspace, reference_file
    ) -> None:
        text = (
            META
            + HEADER.format(s="s1")
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tDP:AD:GT:GQ\t10:6,4:0/1:99\n"
        )
        target = workspace / "s.vcf"
        target.write_text(text)
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 0
        assert err == ""
        assert json.loads(out)["samples"][0]["gt"] == "0/1"

    def test_vcf_parse_error_reports_source_and_line(
        self, workspace, reference_file
    ) -> None:
        target = workspace / "s.vcf"
        target.write_text(META + "chr1\t2\t.\tC\tG\n")  # no #CHROM header
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert str(target) in err

    def test_duplicate_normalized_key_in_one_vcf(
        self, workspace, reference_file
    ) -> None:
        write_vcf(
            workspace / "s.vcf",
            "s1",
            [
                row(pos=4, ref="GT", alt="T", gt="1/1", dp=8, ad=(0, 8)),
                row(pos=2, ref="CG", alt="C", gt="1/1", dp=4, ad=(0, 4)),
            ],
        )
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert "duplicate" in err

    # -- reference semantics -------------------------------------------------

    def test_ref_mismatch_exit_2(self, workspace, reference_file) -> None:
        write_vcf(workspace / "s.vcf", "s1", [row(pos=3, ref="C", alt="T")])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert "POS 3" in err

    def test_missing_chrom_exit_2(self, workspace, reference_file) -> None:
        write_vcf(workspace / "s.vcf", "s1", [row(chrom="chr9", ref="A", alt="G")])
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"s.vcf"}\n')
        code, _, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 2
        assert "chr9" in err

    # -- file handling -------------------------------------------------------

    def test_missing_manifest_file_exit_1(self, workspace, reference_file) -> None:
        code, out, err = run(
            [
                "summarize-variants", str(workspace / "nope.jsonl"),
                "--reference", str(reference_file),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_vcf_file_exit_1(self, workspace, reference_file) -> None:
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"vcfs/gone.vcf"}\n')
        code, out, err = run(
            ["summarize-variants", str(manifest), "--reference", str(reference_file)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_file_exit_1(self, workspace) -> None:
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"s1","vcf":"x.vcf"}\n')
        code, out, err = run(
            [
                "summarize-variants", str(manifest),
                "--reference", str(workspace / "gone.fa"),
            ]
        )
        assert code == 1
        assert out == ""

    def test_failure_preserves_existing_output(
        self, workspace, reference_file
    ) -> None:
        write_vcf(workspace / "good.vcf", "s1", [row()])
        # bad.vcf has an invalid genotype; it appears after the good sample.
        bad = workspace / "bad.vcf"
        bad.write_text(
            META
            + HEADER.format(s="s2")
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT:DP:AD\t0/2:3:1,2\n"
        )
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            '{"sample":"s1","vcf":"good.vcf"}\n'
            '{"sample":"s2","vcf":"bad.vcf"}\n'
        )
        target = workspace / "out.jsonl"
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
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"

    def test_reference_required(self) -> None:
        code, _, err = run(["summarize-variants", "-"], "")
        assert code == 2
        assert err.count("\n") == 1
