"""Tests for the ``genome_variant.compare`` public interface."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.compare import (
    VariantComparison,
    VariantComparisonError,
    VariantDifference,
    compare_variants,
    render_comparison,
)
from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import ReferenceMismatchError, VcfFormatError, read_vcf

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"

REF_RECORDS = (
    SequenceRecord("chr1", "ACGGT"),
    SequenceRecord("chr2", "TTTTT"),
)


def vcf(sample: str, rows: list[str], meta: str = META) -> str:
    return meta + HEADER.format(s=sample) + "".join(row + "\n" for row in rows)


def row(chrom="chr1", pos=2, ref="C", alt="G", flt="PASS", gt="0/1",
        sample_format="GT", sample_value=None) -> str:
    value = gt if sample_value is None else sample_value
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\t{sample_format}\t{value}"
    )


def parse(text: str):
    return read_vcf(StringIO(text))


def compare(base_text, cand_text, reference=REF_RECORDS):
    return compare_variants(parse(base_text), parse(cand_text), reference)


class TestComparisonCategories:
    def test_exact_match_counts_as_matched(self) -> None:
        result = compare(
            vcf("S1", [row(gt="0/1")]),
            vcf("S2", [row(gt="0/1")]),
        )
        assert result == VariantComparison(
            baseline_sample="S1",
            candidate_sample="S2",
            total=1,
            matched=1,
            genotype_mismatch=0,
            baseline_only=0,
            candidate_only=0,
            concordance="1.000000",
            differences=(),
        )

    def test_same_key_different_gt_is_genotype_mismatch(self) -> None:
        result = compare(
            vcf("S1", [row(gt="0/1")]),
            vcf("S2", [row(gt="1/1")]),
        )
        assert result.matched == 0
        assert result.genotype_mismatch == 1
        assert result.total == 1
        assert result.concordance == "0.000000"
        assert result.differences == (
            VariantDifference(
                "chr1", 2, "C", "G",
                "genotype_mismatch", "0/1", "1/1",
            ),
        )

    def test_only_categories_keep_one_sided_gt(self) -> None:
        result = compare(
            vcf("S1", [row(pos=2, gt="0/1")]),
            vcf("S2", [row(pos=3, ref="G", alt="A", gt="1/1")]),
        )
        assert result.baseline_only == 1
        assert result.candidate_only == 1
        assert result.total == 2
        assert result.concordance == "0.000000"
        statuses = {d.status for d in result.differences}
        assert statuses == {"baseline_only", "candidate_only"}
        by_status = {d.status: d for d in result.differences}
        assert by_status["baseline_only"].baseline_gt == "0/1"
        assert by_status["baseline_only"].candidate_gt is None
        assert by_status["candidate_only"].baseline_gt is None
        assert by_status["candidate_only"].candidate_gt == "1/1"

    def test_mixed_categories_and_six_decimal_concordance(self) -> None:
        # 4 matched, 2 mismatched, 1 baseline-only, 1 candidate-only.
        base_rows = [
            row(chrom="chr1", pos=1, ref="A", alt="T", gt="0/1"),
            row(chrom="chr1", pos=2, ref="C", alt="G", gt="0/1"),
            row(chrom="chr1", pos=3, ref="G", alt="A", gt="0/1"),
            row(chrom="chr1", pos=4, ref="G", alt="T", gt="0/1"),
            row(chrom="chr1", pos=5, ref="T", alt="A", gt="0/1"),
            row(chrom="chr2", pos=1, ref="T", alt="A", gt="0/1"),
            row(chrom="chr2", pos=2, ref="T", alt="C", gt="0/1"),
        ]
        cand_rows = [
            row(chrom="chr1", pos=1, ref="A", alt="T", gt="0/1"),
            row(chrom="chr1", pos=2, ref="C", alt="G", gt="0/1"),
            row(chrom="chr1", pos=3, ref="G", alt="A", gt="0/1"),
            row(chrom="chr1", pos=4, ref="G", alt="T", gt="0/1"),
            row(chrom="chr1", pos=5, ref="T", alt="A", gt="1/1"),
            row(chrom="chr2", pos=1, ref="T", alt="A", gt="1/1"),
            row(chrom="chr2", pos=3, ref="T", alt="G", gt="1/1"),
        ]
        result = compare(vcf("S1", base_rows), vcf("S2", cand_rows))
        assert (result.matched, result.genotype_mismatch,
                result.baseline_only, result.candidate_only) == (4, 2, 1, 1)
        assert result.total == 8
        assert result.concordance == f"{4 / 8:.6f}"

    def test_both_empty_gives_perfect_concordance(self) -> None:
        result = compare(vcf("S1", []), vcf("S2", []))
        assert result.total == 0
        assert result.matched == 0
        assert result.concordance == "1.000000"
        assert result.differences == ()


class TestNormalizationSemantics:
    def test_equivalent_representations_match(self) -> None:
        # Baseline: deletion anchored at pos 4 (GT>T); candidate: already
        # minimal at pos 2 (CG>C).  Both normalize to chr1:2 CG>C.
        result = compare(
            vcf("S1", [row(pos=4, ref="GT", alt="T", gt="1/1")]),
            vcf("S2", [row(pos=2, ref="CG", alt="C", gt="1/1")]),
        )
        assert result.matched == 1
        assert result.differences == ()

    def test_input_record_order_does_not_change_result(self) -> None:
        base_rows = [
            row(chrom="chr2", pos=1, ref="T", alt="A"),
            row(chrom="chr1", pos=5, ref="T", alt="A"),
            row(chrom="chr1", pos=2, ref="C", alt="G"),
        ]
        cand_rows = list(reversed(base_rows))
        first = compare(vcf("S1", base_rows), vcf("S2", cand_rows))
        second = compare(vcf("S1", cand_rows), vcf("S2", base_rows))
        assert render_comparison(first) == render_comparison(second)

    def test_differences_sorted_by_reference_order_pos_ref_alt(self) -> None:
        base_rows = [
            row(chrom="chr2", pos=1, ref="T", alt="A"),
            row(chrom="chr1", pos=3, ref="G", alt="T"),
            row(chrom="chr1", pos=3, ref="G", alt="A"),
            row(chrom="chr1", pos=1, ref="A", alt="T"),
        ]
        result = compare(vcf("S1", base_rows), vcf("S2", []))
        keys = [(d.chrom, d.pos, d.ref, d.alt) for d in result.differences]
        assert keys == [
            ("chr1", 1, "A", "T"),
            ("chr1", 3, "G", "A"),
            ("chr1", 3, "G", "T"),
            ("chr2", 1, "T", "A"),
        ]

    def test_duplicate_normalized_key_raises(self) -> None:
        rows = [
            row(pos=4, ref="GT", alt="T", gt="1/1"),
            row(pos=2, ref="CG", alt="C", gt="1/1"),
        ]
        with pytest.raises(VariantComparisonError, match="duplicate normalized"):
            compare(vcf("S1", rows), vcf("S2", []))


class TestRecordConstraints:
    def test_wrong_sample_count_rejected(self) -> None:
        text = (
            META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t0/1\n"
        )
        with pytest.raises(VariantComparisonError, match="exactly one sample"):
            compare(text, vcf("S2", [row()]))

    def test_sample_less_vcf_rejected(self) -> None:
        text = (
            META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t2\t.\tC\tG\t.\tPASS\t.\n"
        )
        with pytest.raises(VariantComparisonError, match="exactly one sample"):
            compare(text, vcf("S2", [row()]))

    @pytest.mark.parametrize("flt", ["LowQual", "PASS;x", ""])
    def test_filter_must_be_pass_or_dot(self, flt) -> None:
        with pytest.raises(VariantComparisonError, match="FILTER"):
            compare(vcf("S1", [row(flt=flt)]), vcf("S2", [row()]))

    @pytest.mark.parametrize("alt", ["G,T", "<DEL>", ".", "*", "G,<DEL>"])
    def test_single_plain_alt_required(self, alt) -> None:
        with pytest.raises(VariantComparisonError):
            compare(vcf("S1", [row(alt=alt)]), vcf("S2", [row()]))

    @pytest.mark.parametrize("gt", ["0/0", "1/2", "0|1", "1|1", "./1", "1/0", ".", ""])
    def test_genotype_subset(self, gt) -> None:
        with pytest.raises(VariantComparisonError, match="GT"):
            compare(vcf("S1", [row(gt=gt)]), vcf("S2", [row()]))

    def test_extra_format_fields_allowed(self) -> None:
        text = vcf(
            "S1",
            [row(sample_format="DP:GT:GQ", sample_value="10:0/1:99")],
        )
        result = compare(text, vcf("S2", [row()]))
        assert result.matched == 1

    def test_missing_gt_field_rejected(self) -> None:
        text = vcf("S1", [row(sample_format="DP", sample_value="10")])
        with pytest.raises(VariantComparisonError, match="GT"):
            compare(text, vcf("S2", [row()]))

    def test_duplicate_format_ids_rejected(self) -> None:
        text = vcf("S1", [row(sample_format="GT:GT", sample_value="0/1:0/1")])
        with pytest.raises(VariantComparisonError, match="duplicate field"):
            compare(text, vcf("S2", [row()]))

    def test_sample_field_count_must_match_format(self) -> None:
        text = vcf("S1", [row(sample_format="GT:DP", sample_value="0/1")])
        with pytest.raises(VariantComparisonError, match="fields but FORMAT"):
            compare(text, vcf("S2", [row()]))

    def test_dot_filter_is_accepted(self) -> None:
        result = compare(
            vcf("S1", [row(flt=".")]), vcf("S2", [row(flt=".")])
        )
        assert result.matched == 1


class TestErrors:
    def test_malformed_vcf_raises_vcf_format_error(self) -> None:
        with pytest.raises(VcfFormatError):
            compare("not a vcf\n", vcf("S2", [row()]))

    def test_record_out_of_bounds_raises_vcf_format_error(self) -> None:
        with pytest.raises(VcfFormatError, match="out of bounds"):
            compare(
                vcf("S1", [row(pos=9, ref="C", alt="G")]),
                vcf("S2", [row()]),
            )

    def test_missing_chrom_raises_reference_mismatch(self) -> None:
        with pytest.raises(ReferenceMismatchError, match="no matching record"):
            compare(
                vcf("S1", [row(chrom="chr9", ref="A", alt="G")]),
                vcf("S2", [row()]),
            )

    def test_duplicate_chrom_raises_reference_mismatch(self) -> None:
        reference = (
            SequenceRecord("chr1", "ACGGT"),
            SequenceRecord("chr1", "ACGGT"),
        )
        with pytest.raises(ReferenceMismatchError, match="more than once"):
            compare(
                vcf("S1", [row()]), vcf("S2", [row()]), reference=reference
            )

    def test_ref_disagreement_raises_reference_mismatch(self) -> None:
        with pytest.raises(ReferenceMismatchError, match="does not match"):
            compare(
                vcf("S1", [row(pos=3, ref="C", alt="T")]),
                vcf("S2", [row()]),
            )

    def test_missing_file_preserves_oserror(self, tmp_path) -> None:
        missing = tmp_path / "gone.vcf"
        with pytest.raises(OSError):
            compare_variants(str(missing), vcf("S2", [row()]), REF_RECORDS)


class TestSources:
    def test_paths_and_streams_accepted(self, tmp_path) -> None:
        base = tmp_path / "b.vcf"
        cand = tmp_path / "c.vcf"
        ref = tmp_path / "r.fa"
        base.write_text(vcf("S1", [row()]))
        cand.write_text(vcf("S2", [row()]))
        ref.write_text(">chr1\nACGGT\n>chr2\nTTTTT\n")
        result = compare_variants(str(base), StringIO(cand.read_text()), str(ref))
        assert result.matched == 1
        assert result.baseline_sample == "S1"
        assert result.candidate_sample == "S2"

    def test_reference_record_iterable_accepted(self) -> None:
        result = compare_variants(
            parse(vcf("S1", [row()])),
            parse(vcf("S2", [row()])),
            iter(REF_RECORDS),
        )
        assert result.matched == 1

    def test_invalid_source_type(self) -> None:
        with pytest.raises(TypeError):
            compare_variants(42, vcf("S2", [row()]), REF_RECORDS)


class TestRendering:
    def test_field_order_and_trailing_newline(self) -> None:
        text = render_comparison(compare(vcf("S1", []), vcf("S2", [])))
        assert text.endswith("\n") and not text.endswith("\n\n")
        assert list(__import__("json").loads(text)) == [
            "baseline_sample",
            "candidate_sample",
            "total",
            "matched",
            "genotype_mismatch",
            "baseline_only",
            "candidate_only",
            "concordance",
            "differences",
        ]

    def test_difference_item_field_order(self) -> None:
        result = compare(vcf("S1", [row(gt="0/1")]), vcf("S2", [row(gt="1/1")]))
        text = render_comparison(result)
        item = __import__("json").loads(text)["differences"][0]
        assert list(item) == [
            "chrom",
            "pos",
            "ref",
            "alt",
            "status",
            "baseline_gt",
            "candidate_gt",
        ]
        assert text == (
            '{"baseline_sample":"S1","candidate_sample":"S2","total":1,'
            '"matched":0,"genotype_mismatch":1,"baseline_only":0,'
            '"candidate_only":0,"concordance":"0.000000","differences":['
            '{"chrom":"chr1","pos":2,"ref":"C","alt":"G",'
            '"status":"genotype_mismatch","baseline_gt":"0/1",'
            '"candidate_gt":"1/1"}]}\n'
        )

    def test_non_ascii_sample_emitted_verbatim(self) -> None:
        text = render_comparison(
            compare(vcf("基线", []), vcf("候选", []))
        )
        assert "基线" in text and "候选" in text
        assert "\\u" not in text

    def test_byte_stable_for_equivalent_inputs(self) -> None:
        a = compare(
            vcf("S1", [row(pos=4, ref="GT", alt="T")]),
            vcf("S2", [row(pos=2, ref="CG", alt="C")]),
        )
        b = compare(
            vcf("S1", [row(pos=2, ref="CG", alt="C")]),
            vcf("S2", [row(pos=4, ref="GT", alt="T")]),
        )
        # Sample names differ, so compare the call-driven parts; the
        # normalized key comparison is symmetric.
        payload_a = __import__("json").loads(render_comparison(a))
        payload_b = __import__("json").loads(render_comparison(b))
        for payload in (payload_a, payload_b):
            payload.pop("baseline_sample")
            payload.pop("candidate_sample")
        assert payload_a == payload_b
        assert payload_a["matched"] == 1
