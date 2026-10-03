"""Tests for the ``compare-variants`` command line entry point."""

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
    return meta + HEADER.format(s=sample) + "".join(row + "\n" for row in rows)


def row(chrom="chr1", pos=2, ref="C", alt="G", flt="PASS", gt="0/1",
        sample_format="GT") -> str:
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\t{sample_format}\t{gt}"
    )


@pytest.fixture
def reference_file(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(REF_FASTA)
    return path


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


class TestCompareVariantsCli:
    def test_basic_comparison_json(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row(gt="0/1")]))
        cand = write(
            tmp_path,
            "c.vcf",
            vcf("S2", [row(gt="0/1"), row(pos=3, ref="G", alt="A", gt="1/1")]),
        )
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 0
        assert err == ""
        data = json.loads(out)
        assert list(data) == [
            "baseline_sample", "candidate_sample", "total", "matched",
            "genotype_mismatch", "baseline_only", "candidate_only",
            "concordance", "differences",
        ]
        assert data["baseline_sample"] == "S1"
        assert data["candidate_sample"] == "S2"
        assert data["total"] == 2
        assert data["matched"] == 1
        assert data["candidate_only"] == 1
        assert data["concordance"] == "0.500000"
        assert data["differences"] == [
            {
                "chrom": "chr1", "pos": 3, "ref": "G", "alt": "A",
                "status": "candidate_only",
                "baseline_gt": None, "candidate_gt": "1/1",
            }
        ]
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_genotype_mismatch_payload(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row(gt="0/1")]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row(gt="1/1")]))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 0
        data = json.loads(out)
        assert data["genotype_mismatch"] == 1
        assert data["matched"] == 0
        difference = data["differences"][0]
        assert list(difference) == [
            "chrom", "pos", "ref", "alt", "status",
            "baseline_gt", "candidate_gt",
        ]
        assert difference["status"] == "genotype_mismatch"
        assert difference["baseline_gt"] == "0/1"
        assert difference["candidate_gt"] == "1/1"

    def test_empty_inputs_concordance_one(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", []))
        cand = write(tmp_path, "c.vcf", vcf("S2", []))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 0
        assert err == ""
        assert out == (
            '{"baseline_sample":"S1","candidate_sample":"S2","total":0,'
            '"matched":0,"genotype_mismatch":0,"baseline_only":0,'
            '"candidate_only":0,"concordance":"1.000000","differences":[]}\n'
        )

    def test_equivalent_normalized_calls_match(self, tmp_path, reference_file) -> None:
        base = write(
            tmp_path, "b.vcf",
            vcf("S1", [row(pos=4, ref="GT", alt="T", gt="1/1")]),
        )
        cand = write(
            tmp_path, "c.vcf",
            vcf("S2", [row(pos=2, ref="CG", alt="C", gt="1/1")]),
        )
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 0
        assert json.loads(out)["matched"] == 1

    def test_differences_sorted_by_reference_then_pos_ref_alt(
        self, tmp_path, reference_file
    ) -> None:
        base = write(
            tmp_path,
            "b.vcf",
            vcf(
                "S1",
                [
                    row(chrom="chr2", pos=1, ref="T", alt="A"),
                    row(chrom="chr1", pos=3, ref="G", alt="T"),
                    row(chrom="chr1", pos=3, ref="G", alt="A"),
                    row(chrom="chr1", pos=1, ref="A", alt="T"),
                ],
            ),
        )
        cand = write(tmp_path, "c.vcf", vcf("S2", []))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 0
        keys = [
            (d["chrom"], d["pos"], d["ref"], d["alt"])
            for d in json.loads(out)["differences"]
        ]
        assert keys == [
            ("chr1", 1, "A", "T"),
            ("chr1", 3, "G", "A"),
            ("chr1", 3, "G", "T"),
            ("chr2", 1, "T", "A"),
        ]

    def test_byte_stable_across_record_order(
        self, tmp_path, reference_file
    ) -> None:
        base_rows = [
            row(chrom="chr2", pos=1, ref="T", alt="A"),
            row(chrom="chr1", pos=3, ref="G", alt="A"),
            row(chrom="chr1", pos=1, ref="A", alt="T"),
        ]
        base_a = write(tmp_path, "a1.vcf", vcf("S1", base_rows))
        base_b = write(tmp_path, "a2.vcf", vcf("S1", list(reversed(base_rows))))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        argv_tail = [str(cand), "--reference", str(reference_file)]
        _, out1, err1 = run(["compare-variants", str(base_a), *argv_tail])
        _, out2, err2 = run(["compare-variants", str(base_b), *argv_tail])
        assert err1 == err2 == ""
        assert out1 == out2

    def test_non_ascii_sample_preserved(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", vcf("基线", []))
        cand = write(tmp_path, "c.vcf", vcf("候选", []))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 0
        assert "基线" in out and "候选" in out
        assert "\\u" not in out

    # -- standard input ------------------------------------------------------

    def test_baseline_from_stdin(self, tmp_path, reference_file) -> None:
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            [
                "compare-variants", "-", str(cand),
                "--reference", str(reference_file),
            ],
            vcf("S1", [row()]),
        )
        assert code == 0
        assert json.loads(out)["matched"] == 1

    def test_reference_from_stdin(self, tmp_path) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            ["compare-variants", str(base), str(cand), "--reference", "-"],
            REF_FASTA,
        )
        assert code == 0
        assert json.loads(out)["matched"] == 1

    @pytest.mark.parametrize(
        "baseline,candidate,reference",
        [("-", "-", "ref.fa"), ("-", "c.vcf", "-"), ("b.vcf", "-", "-")],
    )
    def test_at_most_one_stdin(
        self, baseline, candidate, reference, tmp_path, reference_file
    ) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        ref = "-" if reference == "-" else str(reference_file)
        argv = [
            "compare-variants",
            "-" if baseline == "-" else str(base),
            "-" if candidate == "-" else str(cand),
            "--reference", ref,
        ]
        code, out, err = run(argv, REF_FASTA)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    # -- output handling ------------------------------------------------------

    def test_output_file_roundtrip(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        target = tmp_path / "out.json"
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
                "--output", str(target),
            ]
        )
        assert code == 0
        assert out == "" and err == ""
        data = target.read_bytes()
        assert data.endswith(b"\n") and not data.endswith(b"\n\n")
        json.loads(data)

    def test_failure_preserves_existing_output(
        self, tmp_path, reference_file
    ) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row(flt="LowQual")]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        target = tmp_path / "out.json"
        target.write_text("PREVIOUS\n")
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
                "--output", str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"

    def test_reference_required(self, tmp_path) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, _, err = run(["compare-variants", str(base), str(cand)])
        assert code == 2
        assert err.count("\n") == 1

    # -- data errors ---------------------------------------------------------

    def test_wrong_sample_count_exit_2(self, tmp_path, reference_file) -> None:
        text = (
            META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t0/1\n"
        )
        base = write(tmp_path, "b.vcf", text)
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert "one sample" in err

    @pytest.mark.parametrize(
        "record",
        [
            "chr1\t2\t.\tC\tG\t.\tLowQual\t.\tGT\t0/1",
            "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/0",
            "chr1\t2\t.\tC\tG,T\t.\tPASS\t.\tGT\t0/1",
            "chr1\t2\t.\tC\t<DEL>\t.\tPASS\t.\tGT\t0/1",
            "chr1\t2\t.\tC\t.\t.\tPASS\t.\tGT\t0/1",
            "chr1\t9\t.\tC\tG\t.\tPASS\t.\tGT\t0/1",
        ],
    )
    def test_bad_records_exit_2(
        self, record, tmp_path, reference_file
    ) -> None:
        base = write(tmp_path, "b.vcf", META + HEADER.format(s="S1") + record + "\n")
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_duplicate_normalized_key_exit_2(
        self, tmp_path, reference_file
    ) -> None:
        base = write(
            tmp_path,
            "b.vcf",
            vcf(
                "S1",
                [
                    row(pos=4, ref="GT", alt="T", gt="1/1"),
                    row(pos=2, ref="CG", alt="C", gt="1/1"),
                ],
            ),
        )
        cand = write(tmp_path, "c.vcf", vcf("S2", []))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 2
        assert "duplicate" in err

    def test_malformed_vcf_exit_2(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", "garbage\n")
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert str(base) in err

    def test_missing_chrom_exit_2(self, tmp_path, reference_file) -> None:
        base = write(
            tmp_path, "b.vcf",
            vcf("S1", [row(chrom="chr9", ref="A", alt="G")]),
        )
        cand = write(tmp_path, "c.vcf", vcf("S2", []))
        code, _, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
            ]
        )
        assert code == 2
        assert "chr9" in err

    def test_bad_reference_fasta_exit_2(
        self, tmp_path, reference_file
    ) -> None:
        bad_ref = write(tmp_path, "bad.fa", "not a fasta\n")
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(bad_ref),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    # -- file access ---------------------------------------------------------

    @pytest.mark.parametrize("which", ["baseline", "candidate", "reference"])
    def test_missing_input_file_exit_1(
        self, which, tmp_path, reference_file
    ) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        argv = [
            "compare-variants",
            str(base if which != "baseline" else tmp_path / "gone.vcf"),
            str(cand if which != "candidate" else tmp_path / "gone.vcf"),
            "--reference",
            str(reference_file if which != "reference" else tmp_path / "gone.fa"),
        ]
        code, out, err = run(argv)
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_unwritable_output_exit_1(self, tmp_path, reference_file) -> None:
        base = write(tmp_path, "b.vcf", vcf("S1", [row()]))
        cand = write(tmp_path, "c.vcf", vcf("S2", [row()]))
        code, out, err = run(
            [
                "compare-variants", str(base), str(cand),
                "--reference", str(reference_file),
                "--output", str(tmp_path / "missing-dir" / "out.json"),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1
