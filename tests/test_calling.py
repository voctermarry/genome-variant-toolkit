"""Tests for genome_variant.calling.call_variants."""

from __future__ import annotations

import math

import pytest

from genome_variant.calling import VariantCallingError, call_variants
from genome_variant.quality import ReadQualityError
from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import read_vcf

from io import StringIO


def rec(identifier, sequence, quality=None):
    return SequenceRecord(identifier, sequence, quality=quality)


def q(scores):
    return tuple(scores)


def records_of(document):
    return document.records


def render(document):
    from genome_variant.vcf import render_vcf

    return render_vcf(document)


class TestBasicCalls:
    def test_homozygous_variant(self) -> None:
        refs = [rec("c1", "ACGTACGTAAAA")]
        reads = [
            rec("a", "ACGTTCGTAAAA", q([40] * 12)),
            rec("b", "ACGTTCGTAAAA", q([40] * 12)),
        ]
        doc = call_variants(refs, reads)
        assert len(doc.records) == 1
        record = doc.records[0]
        assert (record.chrom, record.pos, record.ref, record.alt) == (
            "c1",
            5,
            "A",
            ("T",),
        )
        assert record.qual == "."
        assert record.filter == "PASS"
        assert record.info == "DP=2;AC=2;AF=1.000000"
        assert record.format_text == "GT:DP:AD"
        assert record.sample_text == ("1/1:2:0,2",)

    def test_heterozygous_genotype(self) -> None:
        refs = [rec("c1", "ACGTACGTAAAA")]
        reads = [
            rec("a", "ACGTTCGTAAAA", q([40] * 12)),
            rec("b", "ACGTTCGTAAAA", q([40] * 12)),
            rec("c", "ACGTACGTAAAA", q([40] * 12)),
            rec("d", "ACGTACGTAAAA", q([40] * 12)),
        ]
        doc = call_variants(refs, reads)
        record = doc.records[0]
        assert record.info == "DP=4;AC=2;AF=0.500000"
        assert record.sample_text == ("0/1:4:2,2",)

    def test_homozygous_threshold_boundary_is_inclusive(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        # 4 alt / 1 ref = 0.8 -> 1/1; 3 alt / 1 ref = 0.75 -> 0/1.
        homo_reads = [
            rec("r1", "AAAATAAAAAAA", q([40] * 12)),
            rec("r2", "AAAATAAAAAAA", q([40] * 12)),
            rec("r3", "AAAATAAAAAAA", q([40] * 12)),
            rec("r4", "AAAATAAAAAAA", q([40] * 12)),
            rec("r5", "AAAAAAAAAAAA", q([40] * 12)),
        ]
        record = call_variants(refs, homo_reads).records[0]
        assert record.info == "DP=5;AC=4;AF=0.800000"
        assert record.sample_text[0].startswith("1/1:")

        hetero_reads = homo_reads[:3] + homo_reads[4:]
        record = call_variants(refs, hetero_reads).records[0]
        assert record.info == "DP=4;AC=3;AF=0.750000"
        assert record.sample_text[0].startswith("0/1:")

    def test_min_alt_count_filters(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        reads = [
            rec("r1", "AAAATAAAAAAA", q([40] * 12)),
            rec("r2", "AAAAAAAAAAAA", q([40] * 12)),
            rec("r3", "AAAAAAAAAAAA", q([40] * 12)),
            rec("r4", "AAAAAAAAAAAA", q([40] * 12)),
        ]
        # One ALT, fraction 0.25 >= 0.2 but count 1 < 2.
        assert call_variants(refs, reads).records == ()
        record = call_variants(refs, reads, min_alt_count=1).records[0]
        assert record.sample_text == ("0/1:4:3,1",)

    def test_min_alt_fraction_filters(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        reads = [rec("r1", "AAAATAAAAAAA", q([40] * 12))] + [
            rec(f"r{i}", "AAAAAAAAAAAA", q([40] * 12)) for i in range(2, 11)
        ]
        # 1 alt out of 10: count threshold relaxed, fraction 0.1 < 0.2.
        assert call_variants(refs, reads, min_alt_count=1).records == ()
        record = call_variants(
            refs, reads, min_alt_count=1, min_alt_fraction=0.1
        ).records[0]
        assert record.info == "DP=10;AC=1;AF=0.100000"

    def test_alt_tie_breaks_in_acgt_order(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        # Two C alts and two G alts, tied for highest non-ref count -> C.
        reads = [
            rec("c1", "AAAACAAAAAAA", q([40] * 12)),
            rec("c2", "AAAACAAAAAAA", q([40] * 12)),
            rec("g1", "AAAAGAAAAAAA", q([40] * 12)),
            rec("g2", "AAAAGAAAAAAA", q([40] * 12)),
        ]
        record = call_variants(refs, reads).records[0]
        assert record.alt == ("C",)
        assert record.sample_text == ("0/1:4:0,2",)

    def test_only_most_frequent_non_reference_is_alt(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        reads = [
            rec("t1", "AAAATAAAAAAA", q([40] * 12)),
            rec("t2", "AAAATAAAAAAA", q([40] * 12)),
            rec("g1", "AAAAGAAAAAAA", q([40] * 12)),
        ]
        record = call_variants(refs, reads).records[0]
        assert record.alt == ("T",)
        # G observations are neither REF nor selected ALT: DP includes
        # them, AD reports ref and selected ALT only.
        assert record.info == "DP=3;AC=2;AF=0.666667"
        assert record.sample_text == ("0/1:3:0,2",)


class TestEvidenceRules:
    def test_low_quality_bases_do_not_count(self) -> None:
        refs = [rec("c1", "ACGTACGTAAAA")]
        # Both alt reads carry low quality exactly at the variant site.
        reads = [
            rec("a", "ACGTTCGTAAAA", q([40] * 4 + [10] + [40] * 7)),
            rec("b", "ACGTTCGTAAAA", q([40] * 4 + [10] + [40] * 7)),
        ]
        assert call_variants(refs, reads).records == ()
        record = call_variants(refs, reads, min_base_quality=10).records[0]
        assert record.alt == ("T",)

    def test_quality_threshold_is_inclusive(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        reads = [
            rec("a", "AAAATAAAAAAA", q([40] * 4 + [20] + [40] * 7)),
            rec("b", "AAAATAAAAAAA", q([40] * 4 + [20] + [40] * 7)),
        ]
        record = call_variants(refs, reads).records[0]
        assert record.info == "DP=2;AC=2;AF=1.000000"

    def test_unmapped_reads_are_ignored(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        reads = [
            rec("a", "GCGCGCGCGCGC", q([40] * 12)),
            rec("b", "GCGCGCGCGCGC", q([40] * 12)),
        ]
        assert call_variants(refs, reads).records == ()

    def test_ambiguous_read_base_ignored(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        reads = [
            rec("a", "AAAANAAAAAAA", q([40] * 12)),
            rec("b", "AAAANAAAAAAA", q([40] * 12)),
        ]
        assert call_variants(refs, reads).records == ()

    def test_ambiguous_reference_column_ignored(self) -> None:
        refs = [rec("c1", "AAAANAAAAAAA")]
        reads = [
            rec("a", "AAAATAAAAAAA", q([40] * 12)),
            rec("b", "AAAATAAAAAAA", q([40] * 12)),
        ]
        assert call_variants(refs, reads).records == ()

    def test_deletion_columns_contribute_no_evidence(self) -> None:
        # The short read aligns 4=4D4=; only read 2 covers the deleted
        # region, so the C at POS 5 has DP 1, not 2.
        refs = [rec("c1", "ACGTTTTTACGT")]
        reads = [
            rec("gapped", "ACGTACGT", q([40] * 8)),
            rec("full", "ACGTCTTTACGT", q([40] * 12)),
        ]
        doc = call_variants(
            refs, reads, min_alt_count=1, min_alt_fraction=0.0
        )
        assert len(doc.records) == 1
        record = doc.records[0]
        assert (record.pos, record.ref, record.alt) == (5, "T", ("C",))
        assert record.info == "DP=1;AC=1;AF=1.000000"

    def test_insertion_columns_contribute_no_evidence(self) -> None:
        refs = [rec("c1", "ACGTACGT")]
        reads = [rec("long", "ACGTTTTTACGT", q([40] * 12))]
        assert call_variants(refs, reads).records == ()

    def test_unaligned_portions_are_ignored(self) -> None:
        # A short perfect fragment inside a longer read's local
        # alignment: only the aligned reference positions get evidence.
        refs = [rec("c1", "TTTACGTTTT")]
        reads = [
            rec("a", "TTTACGTTTT", q([40] * 10)),
            rec("b", "TTTACGTTTT", q([40] * 10)),
        ]
        doc = call_variants(refs, reads)
        assert doc.records == ()

    def test_one_read_contributes_once_per_position(self) -> None:
        refs = [rec("c1", "ACGTACGTAAAA")]
        # Two distinct reads, each covers POS 5 exactly once.
        reads = [
            rec("a", "ACGTTCGTAAAA", q([40] * 12)),
            rec("b", "ACGTTCGTAAAA", q([40] * 12)),
        ]
        record = call_variants(refs, reads).records[0]
        assert record.info == "DP=2;AC=2;AF=1.000000"

    def test_reverse_strand_quality_maps_to_original_read(self) -> None:
        def rc(sequence):
            return sequence.translate(
                str.maketrans("ACGTN", "TGCAN")
            )[::-1]

        refs = [rec("c1", "AAAAAGAAAAA")]
        # This read only matches through the reverse strand; the variant
        # column on the RC read is index 5, original index 5.
        original = rc("AAAAATAAAAA")
        good = [
            rec("r1", original, q([40] * 11)),
            rec("r2", original, q([40] * 11)),
        ]
        record = call_variants(refs, good).records[0]
        assert (record.pos, record.ref, record.alt) == (6, "G", ("T",))

        bad_quality = [
            rec("r1", original, q([40] * 5 + [5] + [40] * 5)),
            rec("r2", original, q([40] * 5 + [5] + [40] * 5)),
        ]
        assert call_variants(refs, bad_quality).records == ()

    def test_reverse_strand_asymmetric_quality_index(self) -> None:
        def rc(sequence):
            return sequence.translate(
                str.maketrans("ACGTN", "TGCAN")
            )[::-1]

        refs = [rec("c1", "AAAAAGAAAAA")]
        # Mismatch at RC-read column 2 -> original index 8, variant at
        # POS 3 (A>T); the other mismatch (column 5, G>A) needs high
        # quality too, so isolate POS 3 filtering by zeroing index 8.
        original = rc("AATAAAAAAA" + "A")
        assert original[8] == "A"
        reads_good = [
            rec("r1", original, q([40] * 11)),
            rec("r2", original, q([40] * 11)),
        ]
        positions = {r.pos for r in call_variants(refs, reads_good).records}
        assert positions == {3, 6}

        reads_bad = [
            rec("r1", original, q(tuple(5 if i == 8 else 40 for i in range(11)))),
            rec("r2", original, q(tuple(5 if i == 8 else 40 for i in range(11)))),
        ]
        positions = {r.pos for r in call_variants(refs, reads_bad).records}
        assert positions == {6}


class TestOrderingAndHeader:
    def test_records_ordered_by_reference_then_position(self) -> None:
        refs = [
            rec("c1", "AAAAAAAAAAAA"),
            rec("c2", "CCCCCCCCCCCC"),
        ]
        # A>T reads only place on c1, C>T reads only on c2, so each read
        # maps unambiguously regardless of input order.
        reads = [
            rec("c2-p8", "CCCCCCCTCCCC", q([40] * 12)),
            rec("c2-p4", "CCCTCCCCCCCC", q([40] * 12)),
            rec("c1-p8", "AAAAAAATAAAA", q([40] * 12)),
            rec("c1-p4", "AAATAAAAAAAA", q([40] * 12)),
        ]
        doc = call_variants(refs, reads, min_alt_count=1)
        assert [(r.chrom, r.pos) for r in doc.records] == [
            ("c1", 4),
            ("c1", 8),
            ("c2", 4),
            ("c2", 8),
        ]

    def test_default_sample_name_and_header(self) -> None:
        doc = call_variants([rec("c1", "ACGT")], [])
        assert doc.header.samples == ("SAMPLE",)
        assert doc.header.meta_lines[0] == "##fileformat=VCFv4.2"
        assert doc.records == ()
        text = render(doc)
        assert text.endswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

    def test_custom_sample_name(self) -> None:
        doc = call_variants(
            [rec("c1", "ACGT")], [], sample_name="NA12878"
        )
        assert doc.header.samples == ("NA12878",)

    def test_no_candidates_still_writes_full_header(self) -> None:
        doc = call_variants(
            [rec("c1", "AAAAAAAAAAAA")],
            [rec("a", "AAAAAAAAAAAA", q([40] * 12))],
        )
        text = render(doc)
        assert text.startswith("##fileformat=VCFv4.2\n")
        assert text.count("\n") == len(doc.header.meta_lines) + 1
        assert "\n#CHROM" in text

    def test_output_is_readable_by_read_vcf(self) -> None:
        refs = [rec("c1", "ACGTACGTAAAA")]
        reads = [
            rec("a", "ACGTTCGTAAAA", q([40] * 12)),
            rec("b", "ACGTTCGTAAAA", q([40] * 12)),
            rec("c", "ACGTACGTAAAA", q([40] * 12)),
        ]
        text = render(call_variants(refs, reads))
        parsed = read_vcf(StringIO(text))
        assert parsed.header.samples == ("SAMPLE",)
        assert len(parsed.records) == 1
        assert parsed.records[0].ref == "A"
        assert parsed.records[0].alt == ("T",)

    def test_batching_does_not_change_result(self) -> None:
        refs = [rec("c1", "ACGTACGTAAAA")]
        reads = [
            rec(f"r{i}", "ACGTTCGTAAAA", q([40] * 12))
            for i in range(5)
        ]

        def generator():
            for record in reads:
                yield record

        assert render(call_variants(refs, list(reads))) == render(
            call_variants(refs, generator())
        )


class TestValidation:
    def test_thresholds_validated_before_consuming_inputs(self) -> None:
        with pytest.raises(ValueError):
            call_variants(iter([rec("c1", "ACGT")],), [], min_base_quality=94)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_base_quality=-1)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_base_quality=True)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_alt_count=0)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_alt_count=1.5)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_alt_fraction=-0.1)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_alt_fraction=1.1)
        with pytest.raises(ValueError):
            call_variants(
                [rec("c1", "ACGT")], [], min_alt_fraction=math.nan
            )
        with pytest.raises(ValueError):
            call_variants(
                [rec("c1", "ACGT")], [], min_alt_fraction=math.inf
            )
        with pytest.raises(ValueError):
            call_variants(
                [rec("c1", "ACGT")],
                [],
                min_alt_fraction=0.8,
                homozygous_fraction=0.5,
            )
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], min_score=0)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], match_score=0)
        with pytest.raises(ValueError):
            call_variants([rec("c1", "ACGT")], [], gap_open=-1)

    def test_boundary_fractions_accepted(self) -> None:
        refs = [rec("c1", "AAAAAAAAAAAA")]
        homo = [
            rec("a", "AAAATAAAAAAA", q([40] * 12)),
            rec("b", "AAAATAAAAAAA", q([40] * 12)),
        ]
        # Fraction exactly 1.0 with homozygous_fraction=1 is still 1/1.
        record = call_variants(
            refs, homo, min_alt_fraction=0.0, homozygous_fraction=1.0
        ).records[0]
        assert record.sample_text[0].startswith("1/1:")

        # min_alt_fraction=0 admits a call the default 0.2 would reject.
        weak = [
            rec("a", "AAAATAAAAAAA", q([40] * 12)),
            rec("b", "AAAAAAAAAAAA", q([40] * 12)),
            rec("c", "AAAAAAAAAAAA", q([40] * 12)),
            rec("d", "AAAAAAAAAAAA", q([40] * 12)),
        ]
        assert call_variants(refs, weak).records == ()
        record = call_variants(
            refs, weak, min_alt_count=1, min_alt_fraction=0.0
        ).records[0]
        assert record.info == "DP=4;AC=1;AF=0.250000"
        assert record.sample_text[0].startswith("0/1:")

    def test_empty_reference_raises_before_reads_consumed(self) -> None:
        def reads_boom():
            yield None  # pragma: no cover
            raise AssertionError("reads should not be consumed")

        with pytest.raises(VariantCallingError):
            call_variants([], reads_boom())

    def test_duplicate_reference_identifier(self) -> None:
        with pytest.raises(VariantCallingError):
            call_variants([rec("c1", "ACGT"), rec("c1", "TTTT")], [])

    @pytest.mark.parametrize("name", ["", "   ", "\t", "a\tb", "\n"])
    def test_bad_sample_name(self, name) -> None:
        with pytest.raises(VariantCallingError):
            call_variants([rec("c1", "ACGT")], [], sample_name=name)

    def test_fasta_read_without_quality_raises(self) -> None:
        with pytest.raises(ReadQualityError):
            call_variants(
                [rec("c1", "ACGTACGTAAAA")], [rec("a", "ACGTACGTAAAA")]
            )

    def test_quality_length_mismatch_raises(self) -> None:
        with pytest.raises(ReadQualityError):
            call_variants(
                [rec("c1", "ACGTACGTAAAA")],
                [rec("a", "ACGTACGTAAAA", q([40] * 5))],
            )

    def test_out_of_range_quality_raises(self) -> None:
        with pytest.raises(ReadQualityError):
            call_variants(
                [rec("c1", "ACGT")],
                [rec("a", "ACGT", q([40, 40, 94, 40]))],
            )
