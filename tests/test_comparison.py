"""Tests for the pairwise variant comparison API."""

from __future__ import annotations

import json
from io import StringIO

import pytest

from genome_variant.comparison import (
    VariantComparisonError,
    compare_variants,
    render_comparison,
)
from genome_variant.vcf import (
    ReferenceMismatchError,
    VcfFormatError,
    read_vcf,
)

REF_FASTA = ">chr1\nACGGGGGT\n>chr2\nTTTTAAAA\n"

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"


def vcf(sample: str, rows: list[str], header: str = HEADER) -> StringIO:
    text = META + header.format(s=sample)
    text += "".join(row + "\n" for row in rows)
    return StringIO(text)


def row(chrom="chr1", pos=2, ref="C", alt="G", flt="PASS", gt="0/1",
        fmt="GT") -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\t{fmt}\t{gt}"


def reference() -> StringIO:
    return StringIO(REF_FASTA)


class TestCompareVariants:
    def test_all_matched(self) -> None:
        result = compare_variants(
            vcf("s1", [row(), row(chrom="chr2", pos=5, ref="A", alt="G", gt="1/1")]),
            vcf("s2", [row(), row(chrom="chr2", pos=5, ref="A", alt="G", gt="1/1")]),
            reference(),
        )
        assert result.baseline_sample == "s1"
        assert result.candidate_sample == "s2"
        assert result.matched == 2
        assert result.genotype_mismatch == 0
        assert result.baseline_only == 0
        assert result.candidate_only == 0
        assert result.total == 2
        assert result.concordance == "1.000000"
        assert result.differences == ()

    def test_both_empty(self) -> None:
        result = compare_variants(vcf("s1", []), vcf("s2", []), reference())
        assert result.total == 0
        assert result.concordance == "1.000000"
        assert render_comparison(result) == (
            '{"baseline_sample":"s1","candidate_sample":"s2","total":0,'
            '"matched":0,"genotype_mismatch":0,"baseline_only":0,'
            '"candidate_only":0,"concordance":"1.000000","differences":[]}\n'
        )

    def test_all_categories(self) -> None:
        result = compare_variants(
            vcf("s1", [
                row(),
                row(pos=3, ref="G", alt="A", gt="0/1"),
                row(chrom="chr2", pos=1, ref="T", alt="C", gt="1/1"),
            ]),
            vcf("s2", [
                row(),
                row(pos=3, ref="G", alt="A", gt="1/1"),
                row(chrom="chr2", pos=5, ref="A", alt="G", gt="0/1"),
            ]),
            reference(),
        )
        assert result.matched == 1
        assert result.genotype_mismatch == 1
        assert result.baseline_only == 1
        assert result.candidate_only == 1
        assert result.total == 4
        assert result.concordance == "0.250000"
        assert [(d.status, d.baseline_gt, d.candidate_gt) for d in result.differences] == [
            ("genotype_mismatch", "0/1", "1/1"),
            ("baseline_only", "1/1", None),
            ("candidate_only", None, "0/1"),
        ]

    def test_equivalent_representations_match(self) -> None:
        # Insertion of a G in the chr1 G-run (pos 3-7): written at
        # different positions, with shared-context padding, both sides
        # normalize to the same left-aligned minimal key.
        result = compare_variants(
            vcf("s1", [row(pos=4, ref="G", alt="GG")]),
            vcf("s2", [row(pos=6, ref="GG", alt="GGG")]),
            reference(),
        )
        assert result.matched == 1
        assert result.differences == ()

    def test_input_order_does_not_matter(self) -> None:
        rows_a = [
            row(chrom="chr2", pos=5, ref="A", alt="G", gt="1/1"),
            row(),
            row(pos=3, ref="G", alt="A", gt="0/1"),
        ]
        rows_b = [
            row(pos=3, ref="G", alt="A", gt="1/1"),
            row(chrom="chr2", pos=5, ref="A", alt="G", gt="1/1"),
        ]
        first = compare_variants(vcf("s1", rows_a), vcf("s2", rows_b), reference())
        second = compare_variants(
            vcf("s1", list(reversed(rows_a))),
            vcf("s2", list(reversed(rows_b))),
            reference(),
        )
        assert render_comparison(first) == render_comparison(second)

    def test_differences_sorted_by_reference_order(self) -> None:
        result = compare_variants(
            vcf("s1", [
                row(chrom="chr2", pos=5, ref="A", alt="G"),
                row(pos=4, ref="G", alt="A"),
                row(pos=3, ref="G", alt="A"),
            ]),
            vcf("s2", []),
            reference(),
        )
        assert [(d.chrom, d.pos) for d in result.differences] == [
            ("chr1", 3),
            ("chr1", 4),
            ("chr2", 5),
        ]

    def test_parsed_documents_accepted(self) -> None:
        baseline = read_vcf(vcf("s1", [row()]))
        candidate = read_vcf(vcf("s2", [row()]))
        result = compare_variants(baseline, candidate, reference())
        assert result.matched == 1

    def test_reference_as_records(self) -> None:
        from genome_variant.sequence_io import read_sequences

        records = list(read_sequences(StringIO(REF_FASTA), format="fasta"))
        result = compare_variants(vcf("s1", [row()]), vcf("s2", [row()]), records)
        assert result.matched == 1

    def test_duplicate_normalized_keys_rejected(self) -> None:
        # Two representations of the same G insertion in one VCF.
        with pytest.raises(VariantComparisonError, match="duplicate normalized"):
            compare_variants(
                vcf("s1", [row(pos=4, ref="G", alt="GG"), row(pos=5, ref="G", alt="GG")]),
                vcf("s2", []),
                reference(),
            )

    def test_zero_samples_rejected(self) -> None:
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        text = META + header + "chr1\t2\t.\tC\tG\t.\tPASS\t.\n"
        with pytest.raises(VariantComparisonError, match="exactly one sample"):
            compare_variants(StringIO(text), vcf("s2", []), reference())

    def test_two_samples_rejected(self) -> None:
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        text = META + header + "chr1\t2\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n"
        with pytest.raises(VariantComparisonError, match="exactly one sample"):
            compare_variants(StringIO(text), vcf("s2", []), reference())

    def test_non_pass_filter_rejected(self) -> None:
        with pytest.raises(VariantComparisonError, match="FILTER"):
            compare_variants(
                vcf("s1", [row(flt="q10")]), vcf("s2", []), reference()
            )

    def test_dot_filter_accepted(self) -> None:
        result = compare_variants(
            vcf("s1", [row(flt=".")]), vcf("s2", [row()]), reference()
        )
        assert result.matched == 1

    def test_multiple_alts_rejected(self) -> None:
        with pytest.raises(VariantComparisonError, match="exactly one ALT"):
            compare_variants(
                vcf("s1", [row(alt="G,A")]), vcf("s2", []), reference()
            )

    def test_symbolic_alt_rejected(self) -> None:
        with pytest.raises(VariantComparisonError, match="ordinary sequence allele"):
            compare_variants(
                vcf("s1", [row(alt="<DEL>")]), vcf("s2", []), reference()
            )

    def test_bad_genotype_rejected(self) -> None:
        with pytest.raises(VariantComparisonError, match="GT must be 0/1 or 1/1"):
            compare_variants(
                vcf("s1", [row(gt="0/0")]), vcf("s2", []), reference()
            )

    def test_missing_gt_field_rejected(self) -> None:
        with pytest.raises(VariantComparisonError, match="FORMAT must contain GT"):
            compare_variants(
                vcf("s1", [row(fmt="DP", gt="10")]), vcf("s2", []), reference()
            )

    def test_ref_mismatch_raises(self) -> None:
        with pytest.raises(ReferenceMismatchError, match="REF does not match"):
            compare_variants(
                vcf("s1", [row(ref="A")]), vcf("s2", []), reference()
            )

    def test_unknown_chrom_raises(self) -> None:
        with pytest.raises(ReferenceMismatchError, match="no matching record"):
            compare_variants(
                vcf("s1", [row(chrom="chr9")]), vcf("s2", []), reference()
            )

    def test_duplicate_reference_chrom_raises(self) -> None:
        ref = StringIO(">chr1\nACGT\n>chr1\nACGT\n")
        with pytest.raises(ReferenceMismatchError, match="more than once"):
            compare_variants(vcf("s1", [row()]), vcf("s2", []), ref)

    def test_out_of_bounds_record_raises(self) -> None:
        with pytest.raises(VcfFormatError, match="out of bounds"):
            compare_variants(
                vcf("s1", [row(chrom="chr2", pos=8, ref="AA", alt="A")]),
                vcf("s2", []),
                reference(),
            )

    def test_malformed_vcf_raises(self) -> None:
        with pytest.raises(VcfFormatError):
            compare_variants(
                StringIO("##fileformat=VCFv4.2\nchr1\t2\t.\tC\tG\n"),
                vcf("s2", []),
                reference(),
            )

    def test_missing_path_raises_oserror(self, tmp_path) -> None:
        with pytest.raises(OSError):
            compare_variants(
                str(tmp_path / "missing.vcf"), vcf("s2", []), reference()
            )


class TestRenderComparison:
    def test_field_order_and_compactness(self) -> None:
        result = compare_variants(
            vcf("s1", [row(gt="0/1")]),
            vcf("s2", [row(gt="1/1")]),
            reference(),
        )
        text = render_comparison(result)
        assert text == (
            '{"baseline_sample":"s1","candidate_sample":"s2","total":1,'
            '"matched":0,"genotype_mismatch":1,"baseline_only":0,'
            '"candidate_only":0,"concordance":"0.000000","differences":['
            '{"chrom":"chr1","pos":2,"ref":"C","alt":"G",'
            '"status":"genotype_mismatch","baseline_gt":"0/1",'
            '"candidate_gt":"1/1"}]}\n'
        )
        # Round-trips through the JSON parser with the documented values.
        payload = json.loads(text)
        assert payload["concordance"] == "0.000000"
        assert payload["differences"][0]["baseline_gt"] == "0/1"

    def test_byte_identical_for_equivalent_inputs(self) -> None:
        first = compare_variants(
            vcf("s1", [row(pos=4, ref="G", alt="GG"), row()]),
            vcf("s2", []),
            reference(),
        )
        second = compare_variants(
            vcf("s1", [row(), row(pos=6, ref="GG", alt="GGG")]),
            vcf("s2", []),
            reference(),
        )
        assert render_comparison(first) == render_comparison(second)
