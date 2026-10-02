"""Tests for genome_variant.mapping.map_reads."""

from __future__ import annotations

import pytest

from genome_variant.mapping import (
    MappingReferenceError,
    ReadMapping,
    map_reads,
)
from genome_variant.sequence_io import (
    SequenceRecord,
    SequenceValidationError,
)


def rec(identifier, sequence, quality=None):
    return SequenceRecord(identifier, sequence, quality=quality)


def mapped_dict(record, id_, ref_index, ref_id, rs, re_, qs, qe, strand, score, cigar):
    return {
        "record": record,
        "id": id_,
        "mapped": True,
        "reference_record": ref_index,
        "reference": ref_id,
        "reference_start": rs,
        "reference_end": re_,
        "query_start": qs,
        "query_end": qe,
        "strand": strand,
        "score": score,
        "cigar": cigar,
    }


UNMAPPED_KEYS_NULL = {
    "reference_record": None,
    "reference": None,
    "reference_start": None,
    "reference_end": None,
    "query_start": None,
    "query_end": None,
    "strand": None,
    "score": None,
    "cigar": None,
}


class TestBasicMapping:
    def test_forward_hit_against_multi_reference(self) -> None:
        refs = [rec("c1", "ACGTACGT"), rec("c2", "TTTTACGT")]
        result = list(map_reads(refs, [rec("r1", "ACGT")]))
        assert len(result) == 1
        m = result[0]
        assert m.to_dict() == mapped_dict(
            0, "r1", 0, "c1", 0, 4, 0, 4, "+", 8, "4="
        )

    def test_record_indices_follow_read_order(self) -> None:
        refs = [rec("c1", "ACGTACGT"), rec("c2", "TTTTACGT")]
        reads = [rec("a", "ACGT"), rec("b", "TTTT"), rec("c", "ACGT")]
        result = list(map_reads(refs, reads))
        assert [m.record for m in result] == [0, 1, 2]
        assert [m.id for m in result] == ["a", "b", "c"]
        # TTTT occurs on c2 only
        assert result[1].reference_record == 1
        assert result[1].strand == "+"

    def test_unmapped_read_is_kept_with_null_fields(self) -> None:
        # Only G/C bases: neither the read nor its reverse complement
        # matches a TTTT reference above the default threshold.
        refs = [rec("c1", "TTTT")]
        result = list(map_reads(refs, [rec("lost", "GCGC")]))
        assert len(result) == 1
        m = result[0]
        assert m.mapped is False
        assert m.record == 0
        assert m.id == "lost"
        payload = m.to_dict()
        assert list(payload) == [
            "record",
            "id",
            "mapped",
            "reference_record",
            "reference",
            "reference_start",
            "reference_end",
            "query_start",
            "query_end",
            "strand",
            "score",
            "cigar",
        ]
        assert payload == {"record": 0, "id": "lost", "mapped": False, **UNMAPPED_KEYS_NULL}

    def test_unmapped_read_between_mapped_reads(self) -> None:
        refs = [rec("c", "ACGT")]
        reads = [rec("a", "ACGT"), rec("b", "RR"), rec("c", "ACGT")]
        result = list(map_reads(refs, reads))
        assert [m.mapped for m in result] == [True, False, True]
        assert result[1].record == 1

    def test_min_score_threshold_filters(self) -> None:
        refs = [rec("c", "ACGT")]
        assert list(map_reads(refs, [rec("q", "ACGT")], min_score=8))[0].mapped
        result = list(map_reads(refs, [rec("q", "ACGT")], min_score=9))
        assert result[0].mapped is False


class TestReverseStrand:
    def test_reverse_complement_hit(self) -> None:
        # reference GGGACGTGGG holds GGGACGT at the start; a read equal to
        # that fragment's reverse complement (ACGTCCC) maps on '-'.
        refs = [rec("c", "GGGACGTGGG")]
        result = list(map_reads(refs, [rec("q", "ACGTCCC")]))
        m = result[0]
        assert m.strand == "-"
        assert (m.reference_start, m.reference_end) == (0, 7)
        assert m.cigar == "7="
        assert m.score == 14

    def test_minus_query_coordinates_convert_to_original_read(self) -> None:
        # RC of AAAAR is YTTTT; the TTTT fragment aligns to ref[3:7], and
        # on the RC read it occupies positions 1:5, converting back to
        # original-read coordinates 0:4 (the R at original position 4 is
        # excluded).
        refs = [rec("c", "GGGTTTT")]
        m = list(map_reads(refs, [rec("q", "AAAAR")]))[0]
        assert m.strand == "-"
        assert (m.reference_start, m.reference_end) == (3, 7)
        assert (m.query_start, m.query_end) == (0, 4)

    def test_palindrome_read_prefers_plus_strand(self) -> None:
        # ACGT is its own reverse complement; both strands tie on every
        # coordinate and CIGAR, so the forward strand must win.
        refs = [rec("c", "ACGT")]
        m = list(map_reads(refs, [rec("q", "ACGT")]))[0]
        assert m.strand == "+"
        assert (m.query_start, m.query_end) == (0, 4)

    def test_cigar_stays_oriented_to_reverse_complement(self) -> None:
        # read ACGATA reverses to TATCGTA; the TACGTA fragment aligns
        # against TAACGTA with one internal mismatch (position 2), so the
        # CIGAR reads along the reverse-complement read: 2=1X3=.
        refs = [rec("c", "TAACGTTTCC")]
        m = list(map_reads(refs, [rec("q", "ACGATA")]))[0]
        assert m.strand == "-"
        assert (m.reference_start, m.reference_end) == (0, 6)
        assert (m.query_start, m.query_end) == (0, 6)
        assert m.cigar == "2=1X3="
        assert m.score == 7


class TestTieBreaking:
    def test_higher_score_beats_reference_order(self) -> None:
        # c0 shares only a short fragment; c1 carries the full read.
        refs = [rec("c0", "ACGGGGGG"), rec("c1", "TTACGTAC")]
        m = list(map_reads(refs, [rec("q", "ACGT")]))[0]
        assert m.reference_record == 1
        assert m.score == 8

    def test_equal_score_uses_reference_index(self) -> None:
        refs = [rec("c0", "ACGT"), rec("c1", "ACGT")]
        m = list(map_reads(refs, [rec("q", "ACGT")]))[0]
        assert m.reference_record == 0
        assert m.reference == "c0"

    def test_equal_score_uses_smaller_reference_start(self) -> None:
        # read ACG (RC CGT): forward matches ref[0:3], reverse matches
        # ref[1:4], both score 6; the smaller reference start wins.
        refs = [rec("c", "ACGT")]
        m = list(map_reads(refs, [rec("q", "ACG")]))[0]
        assert (m.reference_start, m.reference_end) == (0, 3)
        assert m.strand == "+"
        assert (m.query_start, m.query_end) == (0, 3)


class TestQualityAndValidation:
    def test_fastq_quality_ignored(self) -> None:
        refs = [rec("c", "ACGT")]
        without = list(map_reads(refs, [rec("q", "ACGT")]))[0]
        with_quality = list(
            map_reads(refs, [rec("q", "ACGT", quality=(0, 0, 0, 0))])
        )[0]
        assert with_quality.score == without.score == 8
        assert with_quality.cigar == without.cigar

    def test_invalid_reference_base_raises(self) -> None:
        refs = [rec("c", "AZ")]
        with pytest.raises(SequenceValidationError):
            list(map_reads(refs, [rec("q", "A")]))

    def test_invalid_read_base_raises(self) -> None:
        refs = [rec("c", "AA")]
        with pytest.raises(SequenceValidationError):
            list(map_reads(refs, [rec("q", "AZ")]))

    def test_invalid_read_base_reached_lazily(self) -> None:
        refs = [rec("c", "ACGT")]
        good_then_bad = iter([rec("ok", "ACGT"), rec("bad", "AZ")])
        iterator = map_reads(refs, good_then_bad)
        first = next(iterator)
        assert first.mapped is True
        with pytest.raises(SequenceValidationError):
            next(iterator)


class TestParameters:
    def test_empty_references_raise_mapping_reference_error(self) -> None:
        with pytest.raises(MappingReferenceError):
            list(map_reads([], [rec("q", "A")]))

    def test_empty_references_before_reads_consumed(self) -> None:
        def boom():
            raise AssertionError("reads must not be consumed")
            yield  # pragma: no cover

        with pytest.raises(MappingReferenceError):
            list(map_reads([], boom()))

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"match_score": 0},
            {"match_score": True},
            {"mismatch_penalty": -1},
            {"mismatch_penalty": False},
            {"gap_open": -1},
            {"gap_extend": -2},
            {"min_score": 0},
            {"min_score": -1},
            {"min_score": True},
            {"min_score": 1.5},
        ],
    )
    def test_invalid_parameters_raise_value_error_before_consumption(
        self, kwargs
    ) -> None:
        def boom():
            raise AssertionError("inputs must not be consumed")
            yield  # pragma: no cover

        # references iterable would raise on consumption; the parameter
        # check happens even earlier.
        with pytest.raises(ValueError):
            list(map_reads(boom(), [rec("q", "A")], **kwargs))

    def test_min_score_default_is_one(self) -> None:
        # A single matching base scores 2 >= 1, so the read maps.
        m = list(map_reads([rec("c", "A")], [rec("q", "TA")]))[0]
        assert m.mapped is True
        assert m.score == 2

    def test_batching_does_not_change_results(self) -> None:
        refs = [rec("c0", "ACGTACGT"), rec("c1", "GGGTTTT")]
        reads = [
            rec(f"r{i}", s)
            for i, s in enumerate(
                ["ACGT", "ACGTCCC", "AAAAR", "TTTT", "GG", "AC", "ACGTAC"]
            )
        ]
        first = [m.to_dict() for m in map_reads(list(refs), list(reads))]
        second = [
            m.to_dict() for m in map_reads(iter(list(refs)), iter(list(reads)))
        ]
        assert first == second

    def test_result_is_frozen_dataclass(self) -> None:
        m = list(map_reads([rec("c", "ACGT")], [rec("q", "ACGT")]))[0]
        with pytest.raises(Exception):
            m.score = 0  # type: ignore[misc]
