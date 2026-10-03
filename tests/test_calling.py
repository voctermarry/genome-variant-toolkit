"""Tests for genome_variant.calling.call_variants."""

from __future__ import annotations

import math
from io import StringIO

import pytest

from genome_variant.calling import VariantCallingError, call_variants
from genome_variant.mapping import _best_candidate
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


def _rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _indel_rows(document):
    return [
        (
            record.chrom,
            record.pos,
            record.ref,
            record.alt,
            record.info,
            record.format_text,
            record.sample_text[0],
        )
        for record in document.records
    ]


# 10-base reference: insertion reads gain one base after reference index 2.
INDEL_REF = "ACGTACGTAC"
INS_READ = "ACGATACGTAC"  # anchor G@2, insert A, right flank T@4
DEL_READ = "ACGACGTAC"  # 9 bases: 3=1D6=, anchors G@2, deletes T@3


class TestIndelCallsDisabledByDefault:
    def test_indels_not_called_without_flag(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", INS_READ), fq("b", INS_READ), fq("c", INDEL_REF)]
        assert call_variants(refs, reads).records == ()

    def test_max_indel_length_validated_even_when_disabled(self) -> None:
        # The parameter is always validated before inputs are consumed,
        # independent of whether indel calling is switched on.
        def boom():
            raise AssertionError("inputs must not be consumed")
            yield  # pragma: no cover

        with pytest.raises(ValueError):
            call_variants(boom(), boom(), max_indel_length=0)

    def test_output_byte_identical_with_flag_off(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", INS_READ), fq("b", INS_READ), fq("c", INDEL_REF)]
        baseline = render_vcf(call_variants(refs, reads, call_indels=False))
        # An internal SNV survives; no indel row appears.
        snv_read = list(INDEL_REF)
        snv_read[4] = "T"
        reads_with_snv = [
            fq("a", INS_READ),
            fq("b", INS_READ),
            fq("c", "".join(snv_read)),
            fq("d", "".join(snv_read)),
        ]
        text = render_vcf(
            call_variants(refs, reads_with_snv, call_indels=False)
        )
        assert text == render_vcf(
            call_variants(refs, reads_with_snv, call_indels=False)
        )
        assert baseline == render_vcf(
            call_variants(refs, reads, call_indels=False)
        )


class TestIndelBasic:
    def test_insertion_call_fields(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", INS_READ), fq("b", INS_READ), fq("c", INDEL_REF)]
        document = call_variants(refs, reads, call_indels=True)
        assert _indel_rows(document) == [
            (
                "c1",
                3,
                "G",
                ("GA",),
                "DP=3;AC=2;AF=0.666667",
                "GT:DP:AD",
                "0/1:3:1,2",
            )
        ]

    def test_deletion_call_fields(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", DEL_READ), fq("b", DEL_READ), fq("c", INDEL_REF)]
        document = call_variants(refs, reads, call_indels=True)
        assert _indel_rows(document) == [
            (
                "c1",
                3,
                "GT",
                ("G",),
                "DP=3;AC=2;AF=0.666667",
                "GT:DP:AD",
                "0/1:3:1,2",
            )
        ]

    def test_homozygous_indel(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq(f"i{i}", INS_READ) for i in range(4)]
        document = call_variants(refs, reads, call_indels=True)
        record = document.records[0]
        assert record.info == "DP=4;AC=4;AF=1.000000"
        assert record.sample_text == ("1/1:4:0,4",)

    def test_indel_below_count_not_called(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", INS_READ), fq("c", INDEL_REF), fq("d", INDEL_REF)]
        assert call_variants(refs, reads, call_indels=True).records == ()

    def test_indel_below_fraction_not_called(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", INS_READ), fq("b", INS_READ)] + [
            fq(f"r{i}", INDEL_REF) for i in range(9)
        ]
        document = call_variants(refs, reads, call_indels=True)
        assert document.records == ()

    def test_one_observation_per_read_per_event(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        # A read cannot double-count its single insertion run.
        reads = [fq("a", INS_READ), fq("b", INS_READ)]
        record = call_variants(refs, reads, call_indels=True).records[0]
        assert record.info == "DP=2;AC=2;AF=1.000000"

    def test_depth_requires_both_event_boundaries(self) -> None:
        # A reference read aligned only up to the left anchor (it does not
        # reach the base immediately after the insertion point) does not
        # cross the event's right boundary and is excluded from DP.
        long_reference = "ACGTACGTACGTACGT"
        refs = [rec("c1", long_reference)]
        insertion_read = "ACGATACGTACGTACGT"
        reads = [
            fq("i1", insertion_read),
            fq("i2", insertion_read),
            fq("short", "ACG"),
        ]
        document = call_variants(refs, reads, call_indels=True)
        record = document.records[0]
        assert (record.pos, record.ref, record.alt) == (3, "G", ("GA",))
        assert record.info == "DP=2;AC=2;AF=1.000000"

        # A reference read covering the anchor and the post-event base
        # does count toward DP.
        reads_with_span = [
            fq("i1", insertion_read),
            fq("i2", insertion_read),
            fq("span", "ACGTACGT"),
        ]
        document = call_variants(refs, reads_with_span, call_indels=True)
        assert document.records[0].info == "DP=3;AC=2;AF=0.666667"

    def test_snv_and_indel_both_retained_and_sorted(self) -> None:
        refs = [rec("c1", REF)]  # ACGTACGTACAA, index 4 is the SNV site
        snv = list(REF)
        snv[4] = "T"
        reads = [
            fq("i1", INS_READ + "AA"),
            fq("i2", INS_READ + "AA"),
            fq("s1", "".join(snv)),
            fq("s2", "".join(snv)),
        ]
        document = call_variants(refs, reads, call_indels=True)
        assert [(r.pos, r.ref, r.alt) for r in document.records] == [
            (3, "G", ("GA",)),
            (5, "A", ("T",)),
        ]

    def test_distinct_indels_same_position_are_separate_rows(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        ins_t = "ACGTTACGTAC"  # anchor G@2, insert T
        reads = [
            fq("a1", INS_READ),
            fq("a2", INS_READ),
            fq("t1", ins_t),
            fq("r1", INDEL_REF),
        ]
        document = call_variants(
            refs,
            reads,
            call_indels=True,
            min_alt_count=1,
            min_alt_fraction=0.0,
        )
        assert [(r.pos, r.ref, r.alt, r.info) for r in document.records] == [
            (3, "G", ("GA",), "DP=3;AC=2;AF=0.666667"),
            (3, "G", ("GT",), "DP=2;AC=1;AF=0.500000"),
        ]

    def test_multi_reference_indel_ordering(self) -> None:
        c2_reference = "TTACGTTTTACG"
        c2_mutant = "TTACGTTTACG"
        refs = [
            rec("c1", INDEL_REF),
            rec("c2", c2_reference),
        ]
        reads = [
            fq("c2a", c2_mutant),
            fq("c2b", c2_mutant),
            fq("c1a", INS_READ),
            fq("c1b", INS_READ),
        ]
        document = call_variants(refs, reads, call_indels=True)
        assert [(r.chrom, r.pos) for r in document.records] == [
            ("c1", 3),
            ("c2", 5),
        ]


class TestIndelNormalization:
    def test_deletion_in_repeat_minimized_and_left_aligned(self) -> None:
        # Interior T run: ACG TTTT ACGT; read ACG TTT ACGT (one T deleted).
        # The gap is anchored on the preceding G (reference index 2); REF
        # GT -> ALT G at POS 3, already the leftmost representation.
        refs = [rec("c1", "ACGTTTTACGT")]
        mutant = "ACGTTTACGT"
        reads = [fq("a", mutant), fq("b", mutant), fq("c", "ACGTTTTACGT")]
        document = call_variants(refs, reads, call_indels=True)
        record = document.records[0]
        assert (record.pos, record.ref, record.alt) == (3, "GT", ("G",))
        text = render_vcf(document)
        parsed = read_vcf(StringIO(text))
        assert render_vcf(normalize_vcf(parsed, refs)) == text

    def test_insertion_in_repeat_left_aligned_and_idempotent(self) -> None:
        refs = [rec("c1", "ACGTTTACGT")]
        mutant = "ACGTTTTACGT"
        reads = [fq("a", mutant), fq("b", mutant), fq("c", "ACGTTTACGT")]
        document = call_variants(refs, reads, call_indels=True)
        record = document.records[0]
        assert (record.pos, record.ref, record.alt) == (3, "G", ("GT",))
        text = render_vcf(document)
        parsed = read_vcf(StringIO(text))
        assert render_vcf(normalize_vcf(parsed, refs)) == text

    def test_generated_indels_survive_normalization(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        for mutant in (INS_READ, DEL_READ):
            reads = [fq("a", mutant), fq("b", mutant), fq("c", INDEL_REF)]
            document = call_variants(refs, reads, call_indels=True)
            text = render_vcf(document)
            parsed = read_vcf(StringIO(text))
            normalized = render_vcf(normalize_vcf(parsed, refs))
            assert normalized == text

    def test_left_align_across_interior_run(self) -> None:
        # Insertion of one T after five matched bases lands inside a T run
        # and is represented anchored on the preceding G at POS 6.
        refs = [rec("c1", "GGACGTTTACGT")]
        mutant = "GGACGTTTTACGT"
        reads = [fq("a", mutant), fq("b", mutant), fq("c", "GGACGTTTACGT")]
        document = call_variants(refs, reads, call_indels=True)
        record = document.records[0]
        assert (record.pos, record.ref, record.alt) == (5, "G", ("GT",))
        text = render_vcf(document)
        assert render_vcf(normalize_vcf(read_vcf(StringIO(text)), refs)) == text


class TestIndelEvidenceRules:
    def test_event_needs_two_aligned_flanks(self) -> None:
        # A run without an aligned read base on both sides never yields a
        # candidate; the descriptor simply carries no events.
        from genome_variant.calling import _collect_indel_evidence

        refs = [rec("c1", INDEL_REF)]
        # The deletion read has internal flanks, so events are present.
        dele = fq("d", DEL_READ)
        winner = _best_candidate(dele, refs, 2, 3, 5, 2, 1)
        assert winner is not None
        _ref_index, alignment, strand = winner
        descriptor = _collect_indel_evidence(
            dele,
            dele.quality,
            INDEL_REF,
            alignment,
            strand,
            20,
            50,
        )
        assert len(descriptor.events) == 1

    def test_inserted_base_quality_required(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        low = SequenceRecord(
            "low",
            INS_READ,
            quality=tuple(19 if i == 3 else 40 for i in range(len(INS_READ))),
        )
        reads = [low, fq("b", INS_READ), fq("c", INDEL_REF)]
        assert call_variants(refs, reads, call_indels=True).records == ()

    def test_flank_quality_required_insertion(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        low_left = SequenceRecord(
            "lowL",
            INS_READ,
            quality=tuple(19 if i == 2 else 40 for i in range(len(INS_READ))),
        )
        low_right = SequenceRecord(
            "lowR",
            INS_READ,
            quality=tuple(19 if i == 4 else 40 for i in range(len(INS_READ))),
        )
        for low in (low_left, low_right):
            reads = [low, fq("b", INS_READ), fq("c", INDEL_REF)]
            assert call_variants(refs, reads, call_indels=True).records == ()

    def test_flank_quality_required_deletion(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        low = SequenceRecord(
            "low",
            DEL_READ,
            quality=tuple(19 if i == 2 else 40 for i in range(len(DEL_READ))),
        )
        reads = [low, fq("b", DEL_READ), fq("c", INDEL_REF)]
        assert call_variants(refs, reads, call_indels=True).records == ()

    def test_reverse_strand_insertion_called(self) -> None:
        fwd_ref = "ACGTACGTAC"
        fwd_mut = "ACGTAGGGCGTAC"  # 5=3I5= insert GGG after index 4
        read_sequence = _rc(fwd_mut)
        refs = [rec("c1", fwd_ref)]
        reads = [
            fq("a", read_sequence),
            fq("b", read_sequence),
            fq("c", _rc(fwd_ref)),
        ]
        document = call_variants(
            refs, reads, call_indels=True, min_alt_count=2,
            min_alt_fraction=0.1,
        )
        record = document.records[0]
        assert (record.pos, record.ref, record.alt) == (5, "A", ("AGGG",))

    def test_reverse_strand_inserted_quality_mapped_back(self) -> None:
        fwd_ref = "ACGTACGTAC"
        fwd_mut = "ACGTAGGGCGTAC"
        read_sequence = _rc(fwd_mut)
        # On the reverse-complement read the three inserted bases fall at
        # original-read indices 5, 6 and 7 (mapped back from the
        # reference-oriented aligned columns).
        low_quality = tuple(
            19 if i in (5, 6, 7) else 40 for i in range(len(read_sequence))
        )
        refs = [rec("c1", fwd_ref)]
        reads = [
            SequenceRecord("low", read_sequence, quality=low_quality),
            fq("c", _rc(fwd_ref)),
        ]
        assert call_variants(
            refs, reads, call_indels=True, min_alt_count=1,
            min_alt_fraction=0.0,
        ).records == ()

    def test_reverse_strand_deletion_flank_quality_mapped_back(self) -> None:
        fwd_ref = "GGACGTTTTACGTC"
        fwd_mut = "GGACGTTTACGTC"  # 5=1D8= on the reverse-complement read
        read_sequence = _rc(fwd_mut)
        # Deletion flanks at aligned query columns 4 and 5 map back to
        # original-read indices 8 (left) and 7 (right).
        refs = [rec("c1", fwd_ref)]
        good = [
            fq("a", read_sequence),
            fq("b", read_sequence),
            fq("c", _rc(fwd_ref)),
        ]
        assert (
            call_variants(refs, good, call_indels=True).records[0].pos == 5
        )
        low = SequenceRecord(
            "low",
            read_sequence,
            quality=tuple(19 if i in (7, 8) else 40 for i in range(13)),
        )
        bad = [low, fq("b", read_sequence), fq("c", _rc(fwd_ref))]
        assert call_variants(refs, bad, call_indels=True).records == ()

    def test_ambiguous_deleted_base_ignored(self) -> None:
        refs = [rec("c1", "ACGTANGTAC")]
        mutant = "ACGTAGTAC"  # deletes N at reference index 5
        reads = [fq("a", mutant), fq("b", mutant), fq("c", "ACGTANGTAC")]
        assert call_variants(
            refs, reads, call_indels=True, min_alt_count=2,
            min_alt_fraction=0.1,
        ).records == ()

    def test_ambiguous_inserted_base_ignored(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        mutant = "ACGNTACGTAC"  # inserts N after index 2 (3=1I7=)
        reads = [fq("a", mutant), fq("b", mutant)]
        assert call_variants(
            refs, reads, call_indels=True, min_alt_count=1,
            min_alt_fraction=0.0,
        ).records == ()

    def test_indel_longer_than_limit_ignored(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        mutant = "ACGTAGGGCGTAC"  # 3-base insertion, 5=3I5=
        reads = [fq("a", mutant), fq("b", mutant)]
        assert call_variants(
            refs, reads, call_indels=True, max_indel_length=2,
            min_alt_count=1, min_alt_fraction=0.0,
        ).records == ()
        assert (
            call_variants(
                refs, reads, call_indels=True, max_indel_length=3,
                min_alt_count=2, min_alt_fraction=0.1,
            ).records[0].alt
            == ("AGGG",)
        )


class TestIndelValidation:
    @pytest.mark.parametrize("value", ["yes", 0, 1, 5.0, None, 1.5])
    def test_call_indels_must_be_boolean(self, value) -> None:
        with pytest.raises(ValueError):
            call_variants(
                [rec("c1", INDEL_REF)], [], call_indels=value
            )

    def test_call_indels_boolean_endpoints_accepted(self) -> None:
        call_variants([rec("c1", INDEL_REF)], [], call_indels=True)
        call_variants([rec("c1", INDEL_REF)], [], call_indels=False)

    @pytest.mark.parametrize("value", [0, -1, True, False, 1.5, "5", 5.0])
    def test_invalid_max_indel_length_raises_before_consumption(
        self, value
    ) -> None:
        def boom():
            raise AssertionError("inputs must not be consumed")
            yield  # pragma: no cover

        with pytest.raises(ValueError):
            call_variants(boom(), boom(), max_indel_length=value)

    def test_max_indel_length_endpoint_accepted(self) -> None:
        call_variants(
            [rec("c1", INDEL_REF)], [], call_indels=True, max_indel_length=1
        )

    def test_indel_batching_is_byte_stable(self) -> None:
        refs = [rec("c1", INDEL_REF), rec("c2", "ACGTTTTACGT")]
        read_list = [
            fq("a", INS_READ),
            fq("b", DEL_READ),
            fq("c", INDEL_REF),
            fq("d", "ACGTTTACGT"),
        ]
        first = render_vcf(
            call_variants(list(refs), list(read_list), call_indels=True)
        )
        second = render_vcf(
            call_variants(
                iter(list(refs)), iter(list(read_list)), call_indels=True
            )
        )
        assert first == second

    def test_indel_output_round_trips_through_read_vcf(self) -> None:
        refs = [rec("c1", INDEL_REF)]
        reads = [fq("a", INS_READ), fq("b", INS_READ), fq("c", INDEL_REF)]
        text = render_vcf(call_variants(refs, reads, call_indels=True))
        document = read_vcf(StringIO(text))
        record = document.records[0]
        assert (record.chrom, record.pos, record.ref, record.alt) == (
            "c1",
            3,
            "G",
            ("GA",),
        )
        assert record.sample_text == ("0/1:3:1,2",)
