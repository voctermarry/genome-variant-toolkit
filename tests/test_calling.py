"""Tests for genome_variant.calling.call_variants."""

from __future__ import annotations

import math
from io import StringIO

import pytest

from genome_variant.calling import VariantCallingError, call_variants
from genome_variant.quality import ReadQualityError
from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import normalize_vcf, read_vcf, render_vcf


def rec(identifier, sequence, quality=None):
    return SequenceRecord(identifier, sequence, quality=quality)


def quals(length, score=40):
    return tuple([score] * length)


def fq(identifier, sequence, score=40):
    return rec(identifier, sequence, quals(len(sequence), score))


def fields(document):
    return [
        (
            record.chrom,
            record.pos,
            record.ref,
            record.alt[0],
            record.qual,
            record.filter,
            record.info,
            record.format_text,
            record.sample_text[0],
        )
        for record in document.records
    ]


REF = "ACGTACGTACAA"


def alt_read(base, position=4, identifier="r"):
    sequence = list(REF)
    sequence[position] = base
    return fq(identifier, "".join(sequence))


class TestBasicCalls:
    def test_heterozygous_call_fields(self) -> None:
        refs = [rec("c1", REF)]
        reads = [
            alt_read("T", identifier="a"),
            alt_read("T", identifier="b"),
            fq("c", REF),
        ]
        document = call_variants(refs, reads)
        assert fields(document) == [
            (
                "c1",
                5,
                "A",
                "T",
                ".",
                "PASS",
                "DP=3;AC=2;AF=0.666667",
                "GT:DP:AD",
                "0/1:3:1,2",
            )
        ]

    def test_homozygous_call(self) -> None:
        refs = [rec("c1", REF)]
        reads = [alt_read("T", identifier=str(i)) for i in range(4)]
        document = call_variants(refs, reads)
        record = document.records[0]
        assert record.info == "DP=4;AC=4;AF=1.000000"
        assert record.sample_text == ("1/1:4:0,4",)

    def test_boundary_fractions_are_inclusive(self) -> None:
        # 1 ALT in 5 is exactly 0.2: callable as 0/1 at the default
        # count threshold lowered to 1.
        refs = [rec("c1", REF)]
        reads = [alt_read("T", identifier="a")] + [
            fq(f"r{i}", REF) for i in range(4)
        ]
        document = call_variants(refs, reads, min_alt_count=1)
        record = document.records[0]
        assert record.info == "DP=5;AC=1;AF=0.200000"
        assert record.sample_text[0].startswith("0/1:")

        # 4 ALT in 5 is exactly 0.8: genotype is 1/1.
        reads = [alt_read("T", identifier=str(i)) for i in range(4)] + [
            fq("r", REF)
        ]
        document = call_variants(refs, reads, min_alt_fraction=0.8)
        assert document.records[0].sample_text[0].startswith("1/1:")

    def test_alt_count_below_threshold_not_called(self) -> None:
        refs = [rec("c1", REF)]
        reads = [alt_read("T"), fq("r", REF), fq("s", REF)]
        assert call_variants(refs, reads).records == ()

    def test_alt_fraction_below_threshold_not_called(self) -> None:
        refs = [rec("c1", REF)]
        reads = [alt_read("T"), alt_read("T")] + [
            fq(f"r{i}", REF) for i in range(9)
        ]
        # 2/11 is below 0.2.
        document = call_variants(refs, reads, min_alt_fraction=0.2)
        assert document.records == ()

    def test_alt_tie_breaks_in_acgt_order(self) -> None:
        refs = [rec("c1", REF)]
        reads = (
            [alt_read("C", identifier="c1a"), alt_read("C", identifier="c2a")]
            + [alt_read("G", identifier="g1"), alt_read("G", identifier="g2")]
            + [alt_read("T", identifier="t1"), alt_read("T", identifier="t2")]
        )
        document = call_variants(refs, reads)
        assert len(document.records) == 1
        assert document.records[0].alt == ("C",)
        assert document.records[0].info == "DP=6;AC=2;AF=0.333333"

    def test_records_follow_reference_order_then_position(self) -> None:
        refs = [rec("c1", REF), rec("c2", "GGGACGTGGG")]
        c2_variant = list("GGGACGTGGG")
        c2_variant[3] = "T"
        late = list(REF)
        late[7] = "A"  # interior T->A at reference index 7 (POS 8)
        reads = [
            fq("late", "".join(late)),
            alt_read("T", position=4, identifier="early"),
            fq("x", "".join(c2_variant)),
        ]
        document = call_variants(refs, reads, min_alt_count=1)
        assert [(r.chrom, r.pos) for r in document.records] == [
            ("c1", 5),
            ("c1", 8),
            ("c2", 4),
        ]


class TestEvidence:
    def test_quality_threshold_filters_base_evidence(self) -> None:
        low_quality = SequenceRecord(
            "low",
            alt_read("T").sequence,
            quality=tuple(19 if i == 4 else 40 for i in range(len(REF))),
        )
        refs = [rec("c1", REF)]
        reads = [low_quality, fq("high", alt_read("T").sequence)]
        # The low-quality variant base is absent from DP: only one ALT
        # observation remains, so nothing is called.
        assert call_variants(refs, reads).records == ()

        # At threshold 19 the same base is admitted (threshold inclusive).
        document = call_variants(refs, reads, min_base_quality=19)
        assert document.records[0].info == "DP=2;AC=2;AF=1.000000"

    def test_reverse_strand_evidence_uses_forward_strand_base(self) -> None:
        # From the mapping suite: read ACGATA maps to TAACGTTTCC on the
        # reverse strand (its RC is TATCGT) with CIGAR 2=1X3=.  The
        # mismatch column has reference A against reverse-complement base
        # T, i.e. an A->T call at POS 3 even though the original read
        # carries a T there too (complemented back on the forward strand).
        refs = [rec("c1", "TAACGTTTCC")]
        reads = [fq(f"r{i}", "ACGATA") for i in range(3)]
        document = call_variants(refs, reads)
        assert len(document.records) == 1
        record = document.records[0]
        assert (record.pos, record.ref, record.alt[0]) == (3, "A", "T")
        assert record.info == "DP=3;AC=3;AF=1.000000"

    def test_reverse_strand_quality_maps_to_original_read(self) -> None:
        # Read ACGATA maps to TAACGTTTCC on '-' with the only mismatch at
        # original-read index 3; a low quality there suppresses that
        # read's ALT evidence.
        low = SequenceRecord(
            "low",
            "ACGATA",
            quality=tuple(10 if i == 3 else 40 for i in range(6)),
        )
        refs = [rec("c1", "TAACGTTTCC")]
        reads = [low, fq("b", "ACGATA"), fq("c", "ACGATA")]
        document = call_variants(refs, reads)
        record = document.records[0]
        assert (record.pos, record.alt[0]) == (3, "T")
        assert record.info == "DP=2;AC=2;AF=1.000000"

    def test_ambiguous_bases_never_count(self) -> None:
        ambiguous = list(REF)
        ambiguous[4] = "N"
        refs = [rec("c1", REF)]
        reads = [fq(f"r{i}", "".join(ambiguous)) for i in range(3)]
        assert call_variants(refs, reads).records == ()

    def test_ambiguous_reference_positions_skipped(self) -> None:
        reference = list(REF)
        reference[4] = "N"
        read = list(REF)
        read[4] = "T"
        refs = [rec("c1", "".join(reference))]
        reads = [fq(f"r{i}", "".join(read)) for i in range(3)]
        assert call_variants(refs, reads).records == ()

    def test_insertion_columns_contribute_nothing(self) -> None:
        # Three matching bases, an extra inserted A, then seven matching
        # bases: paying the gap penalty beats restarting at the suffix
        # (10*2-5 = 15 > 7*2 = 14), so the winning CIGAR contains the
        # insertion column, which consumes no reference position.  No SNV
        # evidence may arise from it.
        refs = [rec("c1", "ACGTACGTAC")]
        reads = [
            fq("a", "ACGATACGTAC"),
            fq("b", "ACGATACGTAC"),
            fq("c", "ACGTACGTAC"),
        ]
        assert call_variants(refs, reads).records == ()

    def test_deletion_columns_contribute_nothing(self) -> None:
        # Three matches, one deleted reference base, then six matches:
        # the gap path (9*2-5 = 13) beats the best suffix restart (6*2
        # = 12).  The deletion column consumes only the reference and
        # must never contribute an observation.
        refs = [rec("c1", "ACGTACGTAC")]
        reads = [
            fq("a", "ACGACGTAC"),
            fq("b", "ACGACGTAC"),
            fq("c", "ACGTACGTAC"),
        ]
        assert call_variants(refs, reads).records == ()

    def test_unmapped_reads_ignored(self) -> None:
        refs = [rec("c1", "ATATATATATAT")]
        reads = [fq("a", "GCGCGC"), fq("b", "CGCGCG")]
        document = call_variants(refs, reads)
        assert document.records == ()

    def test_one_observation_per_read_per_position(self) -> None:
        # Two reads give exactly two observations even though the variant
        # base also appears elsewhere in the read.
        refs = [rec("c1", REF)]
        sequence = alt_read("T").sequence
        document = call_variants(refs, [fq("a", sequence), fq("b", sequence)])
        assert document.records[0].info.startswith("DP=2;AC=2")


class TestQualityErrors:
    def test_read_without_quality_raises(self) -> None:
        refs = [rec("c1", REF)]
        with pytest.raises(ReadQualityError):
            call_variants(refs, [rec("a", REF)])

    def test_quality_length_mismatch_raises(self) -> None:
        refs = [rec("c1", REF)]
        bad = SequenceRecord("a", REF, quality=(40, 40))
        with pytest.raises(ReadQualityError):
            call_variants(refs, [bad])

    def test_quality_error_names_bad_record(self) -> None:
        refs = [rec("c1", REF)]
        good = fq("good", REF)
        bad = SequenceRecord("bad", REF, quality=(1, 2))
        with pytest.raises(ReadQualityError) as excinfo:
            call_variants(refs, [good, bad])
        assert "bad" in str(excinfo.value)


class TestValidation:
    def test_empty_reference_raises(self) -> None:
        with pytest.raises(VariantCallingError):
            call_variants([], [fq("a", "ACGT")])

    def test_empty_reference_before_reads_consumed(self) -> None:
        def boom():
            raise AssertionError("reads must not be consumed")
            yield  # pragma: no cover

        with pytest.raises(VariantCallingError):
            call_variants([], boom())

    def test_duplicate_reference_identifier_raises(self) -> None:
        with pytest.raises(VariantCallingError):
            call_variants([rec("x", "ACGT"), rec("x", "TTTT")], [])

    @pytest.mark.parametrize("name", ["", "   ", "\t", "a\tb"])
    def test_bad_sample_name_raises(self, name) -> None:
        with pytest.raises(VariantCallingError):
            call_variants([rec("c1", REF)], [], sample_name=name)

    def test_non_string_sample_name_raises(self) -> None:
        with pytest.raises(VariantCallingError):
            call_variants([rec("c1", REF)], [], sample_name=7)  # type: ignore[arg-type]

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"min_base_quality": -1},
            {"min_base_quality": 94},
            {"min_base_quality": True},
            {"min_base_quality": 1.5},
            {"min_base_quality": "20"},
            {"min_alt_count": 0},
            {"min_alt_count": -2},
            {"min_alt_count": True},
            {"min_alt_count": 1.0},
            {"min_alt_fraction": -0.1},
            {"min_alt_fraction": 1.1},
            {"min_alt_fraction": math.nan},
            {"min_alt_fraction": math.inf},
            {"min_alt_fraction": "0.5"},
            {"min_alt_fraction": True},
            {"homozygous_fraction": -0.1},
            {"homozygous_fraction": math.nan},
            {"max_indel_length": 0},
            {"max_indel_length": -3},
            {"max_indel_length": True},
            {"max_indel_length": 2.5},
            {"max_indel_length": "50"},
        ],
    )
    def test_invalid_thresholds_raise_before_consumption(self, kwargs) -> None:
        def boom():
            raise AssertionError("inputs must not be consumed")
            yield  # pragma: no cover

        with pytest.raises(ValueError):
            call_variants(boom(), boom(), **kwargs)

    def test_homozygous_below_alt_fraction_raises(self) -> None:
        with pytest.raises(ValueError):
            call_variants(
                [rec("c1", REF)],
                [],
                min_alt_fraction=0.7,
                homozygous_fraction=0.5,
            )

    def test_threshold_endpoints_accepted(self) -> None:
        document = call_variants(
            [rec("c1", REF)],
            [],
            min_base_quality=0,
            min_alt_count=1,
            min_alt_fraction=0.0,
            homozygous_fraction=0.0,
        )
        assert document.records == ()
        call_variants(
            [rec("c1", REF)],
            [],
            min_base_quality=93,
            min_alt_fraction=1.0,
            homozygous_fraction=1.0,
        )


class TestOutput:
    def test_default_sample_name_and_header(self) -> None:
        document = call_variants([rec("c1", REF)], [])
        assert document.header.samples == ("SAMPLE",)
        assert document.header.meta_lines[0] == "##fileformat=VCFv4.2"
        text = render_vcf(document)
        assert text.endswith("\n")
        assert not text.endswith("\n\n")
        assert text.splitlines()[-1] == (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        )

    def test_custom_sample_name(self) -> None:
        document = call_variants(
            [rec("c1", REF)], [], sample_name="tumor-sample"
        )
        assert document.header.samples == ("tumor-sample",)

    def test_output_round_trips_through_read_vcf(self) -> None:
        refs = [rec("c1", REF)]
        reads = [alt_read("T"), alt_read("T"), fq("c", REF)]
        text = render_vcf(call_variants(refs, reads))
        document = read_vcf(StringIO(text))
        assert document.header.samples == ("SAMPLE",)
        assert len(document.records) == 1
        record = document.records[0]
        assert (record.chrom, record.pos, record.ref, record.alt) == (
            "c1",
            5,
            "A",
            ("T",),
        )
        assert record.format_text == "GT:DP:AD"
        assert record.sample_text == ("0/1:3:1,2",)

    def test_batching_is_byte_stable(self) -> None:
        refs = [rec("c1", REF), rec("c2", "GGGACGTGGG")]
        c2_variant = list("GGGACGTGGG")
        c2_variant[3] = "T"
        read_list = [
            alt_read("T", identifier="a"),
            fq("b", REF),
            fq("c", "".join(c2_variant)),
        ]
        first = render_vcf(call_variants(list(refs), list(read_list)))
        second = render_vcf(
            call_variants(iter(list(refs)), iter(list(read_list)))
        )
        assert first == second


# A non-repetitive reference; the insertion of "TT" and the deletion of
# the "CG" at zero-based positions 15-16 both left-anchor on the "C" at
# position 14 (POS 15).
INDEL_REF = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"


def ins_read(identifier, inserted="TT", quality=40):
    sequence = INDEL_REF[:15] + inserted + INDEL_REF[15:]
    if isinstance(quality, int):
        quality = quals(len(sequence), quality)
    return rec(identifier, sequence, quality)


def del_read(identifier, quality=40):
    sequence = INDEL_REF[:15] + INDEL_REF[17:]
    if isinstance(quality, int):
        quality = quals(len(sequence), quality)
    return rec(identifier, sequence, quality)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class TestIndelCalls:
    def test_disabled_by_default(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [ins_read("a"), ins_read("b"), fq("w", INDEL_REF)]
        assert call_variants(refs, reads).records == ()
        assert call_variants(refs, reads, call_indels=False).records == ()

    def test_insertion_called(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [ins_read("a"), ins_read("b"), ins_read("c"), fq("w", INDEL_REF)]
        document = call_variants(refs, reads, call_indels=True)
        assert fields(document) == [
            (
                "c1",
                15,
                "C",
                "CTT",
                ".",
                "PASS",
                "DP=4;AC=3;AF=0.750000",
                "GT:DP:AD",
                "0/1:4:1,3",
            )
        ]

    def test_deletion_called(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [del_read("a"), del_read("b"), fq("w", INDEL_REF)]
        document = call_variants(refs, reads, call_indels=True)
        assert fields(document) == [
            (
                "c1",
                15,
                "CCG",
                "C",
                ".",
                "PASS",
                "DP=3;AC=2;AF=0.666667",
                "GT:DP:AD",
                "0/1:3:1,2",
            )
        ]

    def test_homozygous_insertion(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [ins_read(str(i)) for i in range(3)]
        document = call_variants(refs, reads, call_indels=True)
        record = document.records[0]
        assert record.info == "DP=3;AC=3;AF=1.000000"
        assert record.sample_text == ("1/1:3:0,3",)

    def test_reverse_strand_indel_evidence(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        sequence = ins_read("x").sequence
        reads = [
            fq("a", sequence),
            fq("b", reverse_complement(sequence)),
            fq("c", reverse_complement(sequence)),
            fq("w", INDEL_REF),
        ]
        document = call_variants(refs, reads, call_indels=True)
        assert fields(document)[0][1:4] == (15, "C", "CTT")
        assert document.records[0].info == "DP=4;AC=3;AF=0.750000"

    def test_reverse_strand_quality_maps_to_original_read(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        sequence = ins_read("x").sequence
        # The inserted bases sit at query indices 15 and 16 of the aligned
        # reverse complement; convert to original-read coordinates.
        quality = [40] * len(sequence)
        quality[len(sequence) - 1 - 15] = 5
        quality[len(sequence) - 1 - 16] = 5
        reads = [
            ins_read("a"),
            ins_read("b"),
            rec("c", reverse_complement(sequence), tuple(quality)),
            fq("w", INDEL_REF),
        ]
        document = call_variants(refs, reads, call_indels=True)
        # Read c still spans the boundary with good flanks (DP) but its
        # inserted bases fail the quality floor (no ALT evidence).
        assert document.records[0].info == "DP=4;AC=2;AF=0.500000"
        assert document.records[0].sample_text == ("0/1:4:2,2",)

    def test_inserted_base_quality_required(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        quality = [40] * len(ins_read("x").sequence)
        quality[15] = 5  # first inserted base
        reads = [
            ins_read("a"),
            ins_read("b"),
            ins_read("c", quality=tuple(quality)),
            fq("w", INDEL_REF),
        ]
        document = call_variants(refs, reads, call_indels=True)
        assert document.records[0].info == "DP=4;AC=2;AF=0.500000"

    def test_flanking_base_quality_required(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        quality = [40] * len(ins_read("x").sequence)
        quality[14] = 5  # left flanking read base
        reads = [
            ins_read("a"),
            ins_read("b"),
            ins_read("c", quality=tuple(quality)),
            fq("w", INDEL_REF),
        ]
        document = call_variants(refs, reads, call_indels=True)
        # The low-quality flank excludes the read from DP as well.
        assert document.records[0].info == "DP=3;AC=2;AF=0.666667"

    def test_conflicting_indel_excluded_from_depth(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [
            ins_read("a", "TT"),
            ins_read("b", "TT"),
            ins_read("c", "G"),
            ins_read("d", "G"),
            fq("w", INDEL_REF),
        ]
        document = call_variants(refs, reads, call_indels=True)
        # Each event's DP counts its own reads plus the indel-free read;
        # reads carrying the other insertion at the same boundary are
        # excluded.  Records sort by ALT within the shared POS.
        assert fields(document) == [
            (
                "c1",
                15,
                "C",
                "CG",
                ".",
                "PASS",
                "DP=3;AC=2;AF=0.666667",
                "GT:DP:AD",
                "0/1:3:1,2",
            ),
            (
                "c1",
                15,
                "C",
                "CTT",
                ".",
                "PASS",
                "DP=3;AC=2;AF=0.666667",
                "GT:DP:AD",
                "0/1:3:1,2",
            ),
        ]

    def test_left_alignment_merges_equivalent_events(self) -> None:
        # The same insertion into the A-run, produced at different raw
        # positions by reads of different lengths, normalizes to one
        # left-aligned event.
        refs = [rec("c2", "ACGTACAAAATCGAC")]
        reference = "ACGTACAAAATCGAC"
        reads = [
            fq("a", reference[0:6] + "A" + reference[6:15]),
            fq("b", reference[2:10] + "A" + reference[10:15]),
            fq("c", reference[1:8] + "A" + reference[8:14]),
        ]
        document = call_variants(refs, reads, call_indels=True)
        assert fields(document) == [
            (
                "c2",
                6,
                "C",
                "CA",
                ".",
                "PASS",
                "DP=3;AC=3;AF=1.000000",
                "GT:DP:AD",
                "1/1:3:0,3",
            )
        ]

    def test_max_indel_length_filters_long_events(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [
            ins_read("a", "TTTTTT"),
            ins_read("b", "TTTTTT"),
            fq("w", INDEL_REF),
        ]
        assert (
            call_variants(refs, reads, call_indels=True, max_indel_length=5).records
            == ()
        )
        document = call_variants(refs, reads, call_indels=True, max_indel_length=6)
        assert document.records[0].alt == ("CTTTTTT",)

    def test_non_acgt_insertion_ignored(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [ins_read("a", "NN"), ins_read("b", "NN")]
        assert call_variants(refs, reads, call_indels=True).records == ()

    def test_deletion_over_ambiguous_reference_ignored(self) -> None:
        # The deleted fragment contains an N.
        refs = [rec("c1", INDEL_REF[:15] + "N" + INDEL_REF[16:])]
        reads = [del_read("a"), del_read("b")]
        assert call_variants(refs, reads, call_indels=True).records == ()

    def test_deletion_with_ambiguous_anchor_ignored(self) -> None:
        # The left anchor base is an N.
        ambiguous = INDEL_REF[:14] + "N" + INDEL_REF[15:]
        refs = [rec("c1", ambiguous)]
        sequence = ambiguous[:15] + ambiguous[17:]
        reads = [fq("a", sequence), fq("b", sequence)]
        assert call_variants(refs, reads, call_indels=True).records == ()

    def test_thresholds_apply_to_indels(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [ins_read("a"), ins_read("b"), ins_read("c"), fq("w", INDEL_REF)]
        assert (
            call_variants(refs, reads, call_indels=True, min_alt_count=4).records
            == ()
        )
        assert (
            call_variants(
                refs,
                reads,
                call_indels=True,
                min_alt_fraction=0.9,
                homozygous_fraction=0.9,
            ).records
            == ()
        )

    def test_snv_and_indels_share_position_and_sort(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        snv = INDEL_REF[:14] + "A" + INDEL_REF[15:]
        reads = [
            ins_read("a", "TT"),
            ins_read("b", "TT"),
            ins_read("c", "G"),
            ins_read("d", "G"),
            fq("e", snv),
            fq("f", snv),
            fq("w", INDEL_REF),
        ]
        document = call_variants(refs, reads, call_indels=True)
        # SNV first (ALT "A"), then the two insertions by ALT; each
        # indel's DP excludes the reads carrying the other insertion.
        assert [(r.pos, r.ref, r.alt[0]) for r in document.records] == [
            (15, "C", "A"),
            (15, "C", "CG"),
            (15, "C", "CTT"),
        ]
        assert document.records[0].info == "DP=7;AC=2;AF=0.285714"
        assert document.records[1].info == "DP=5;AC=2;AF=0.400000"
        assert document.records[2].info == "DP=5;AC=2;AF=0.400000"

    def test_max_indel_length_endpoint_accepted(self) -> None:
        document = call_variants(
            [rec("c1", INDEL_REF)], [], call_indels=True, max_indel_length=1
        )
        assert document.records == ()

    def test_indel_output_round_trips_and_normalizes(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [
            ins_read("a"),
            ins_read("b"),
            del_read("c"),
            del_read("d"),
            fq("w", INDEL_REF),
        ]
        text = render_vcf(call_variants(refs, reads, call_indels=True))
        document = read_vcf(StringIO(text))
        assert len(document.records) == 2
        normalized = normalize_vcf(document, refs)
        assert render_vcf(normalized) == text

    def test_indel_batching_is_byte_stable(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [ins_read("a"), ins_read("b"), del_read("c"), fq("w", INDEL_REF)]
        first = render_vcf(
            call_variants(list(refs), list(reads), call_indels=True)
        )
        second = render_vcf(
            call_variants(iter(list(refs)), iter(list(reads)), call_indels=True)
        )
        assert first == second
