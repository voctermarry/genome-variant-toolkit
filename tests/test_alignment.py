"""Tests for pairwise alignment (genome_variant.alignment)."""

from __future__ import annotations

import pytest

from genome_variant.alignment import PairwiseAlignment, align_pair
from genome_variant.sequence_io import SequenceRecord, SequenceValidationError


def rec(sequence: str, identifier: str = "x", quality=None) -> SequenceRecord:
    return SequenceRecord(identifier, sequence, quality=quality)


class TestGlobal:
    def test_identical_sequences(self) -> None:
        result = align_pair(rec("ACGT", "r"), rec("ACGT", "q"))
        assert result == PairwiseAlignment(
            reference="r",
            query="q",
            mode="global",
            score=8,
            reference_start=0,
            reference_end=4,
            query_start=0,
            query_end=4,
            cigar="4=",
            aligned_reference="ACGT",
            aligned_query="ACGT",
        )

    def test_mismatch_is_x(self) -> None:
        result = align_pair(
            rec("ACGT", "r"), rec("ACGA", "q"),
            match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2,
        )
        # 3 matches (6) minus one mismatch (3) = 3
        assert result.score == 3
        assert result.cigar == "3=1X"
        assert result.aligned_reference == "ACGT"
        assert result.aligned_query == "ACGA"

    def test_iupac_ambiguity_is_always_a_mismatch(self) -> None:
        result = align_pair(rec("R", "r"), rec("A", "q"), match_score=1, mismatch_penalty=1, gap_open=5, gap_extend=2)
        assert result.cigar == "1X"
        assert result.score == -1
        result = align_pair(rec("N", "r"), rec("N", "q"), match_score=1, mismatch_penalty=1, gap_open=5, gap_extend=2)
        assert result.cigar == "1="
        assert result.score == 1

    def test_deletion_run_is_affine(self) -> None:
        result = align_pair(
            rec("ACGT", "r"), rec("AGT", "q"),
            match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2,
        )
        # 3 matches (6) - opening (5) = 1
        assert result.score == 1
        assert result.cigar == "1=1D2="
        assert result.aligned_reference == "ACGT"
        assert result.aligned_query == "A-GT"

    def test_two_base_gap_pays_open_plus_one_extend(self) -> None:
        result = align_pair(
            rec("ACGT", "r"), rec("AT", "q"),
            match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2,
        )
        assert result.cigar == "1=2D1="
        assert result.score == 4 - (5 + 2)
        assert result.aligned_reference == "ACGT"
        assert result.aligned_query == "A--T"

    def test_insertion_consumes_query_only(self) -> None:
        result = align_pair(rec("AT", "r"), rec("ACGT", "q"), match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2)
        assert result.cigar == "1=2I1="
        assert result.aligned_reference == "A--T"
        assert result.aligned_query == "ACGT"
        assert (result.reference_start, result.reference_end) == (0, 2)
        assert (result.query_start, result.query_end) == (0, 4)

    def test_covers_full_sequences_with_empty_side(self) -> None:
        result = align_pair(rec("", "r"), rec("A", "q"), match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2)
        assert result.score == -5
        assert result.cigar == "1I"
        assert result.aligned_reference == "-"
        assert result.aligned_query == "A"
        assert result.reference_end == 0
        assert result.query_end == 1

        result = align_pair(rec("A", "r"), rec("", "q"), match_score=2, mismatch_penalty=3, gap_open=5, gap_extend=2)
        assert result.cigar == "1D"
        assert result.aligned_reference == "A"
        assert result.aligned_query == "-"

    def test_both_empty(self) -> None:
        result = align_pair(rec("", "r"), rec("", "q"))
        assert result.score == 0
        assert result.cigar == ""
        assert result.aligned_reference == ""
        assert result.aligned_query == ""

    def test_zero_penalties_score_zero(self) -> None:
        result = align_pair(
            rec("AC", "r"), rec("GT", "q"),
            match_score=1, mismatch_penalty=0, gap_open=0, gap_extend=0,
        )
        # Every path scores zero; only the score and full coverage matter.
        assert result.score == 0
        assert result.aligned_reference.replace("-", "") == "AC"
        assert result.aligned_query.replace("-", "") == "GT"
        assert len(result.aligned_reference) == len(result.aligned_query)

    def test_quality_values_are_ignored(self) -> None:
        fastq = SequenceRecord("q", "ACGT", quality=(1, 2, 3, 4))
        result = align_pair(rec("ACGT", "r"), fastq)
        assert result.score == 8
        assert result.cigar == "4="


class TestLocal:
    def test_best_positive_fragment(self) -> None:
        result = align_pair(
            rec("CCACGTCC", "r"), rec("GGACGTGG", "q"),
            mode="local", match_score=2, mismatch_penalty=3,
            gap_open=5, gap_extend=2,
        )
        assert result.score == 8
        assert (result.reference_start, result.reference_end) == (2, 6)
        assert (result.query_start, result.query_end) == (2, 6)
        assert result.cigar == "4="
        assert result.aligned_reference == "ACGT"
        assert result.aligned_query == "ACGT"

    def test_no_positive_fragment_returns_empty_alignment(self) -> None:
        result = align_pair(
            rec("AC", "r"), rec("GT", "q"),
            mode="local", match_score=2, mismatch_penalty=3,
            gap_open=5, gap_extend=2,
        )
        assert result.mode == "local"
        assert result.score == 0
        assert (
            result.reference_start,
            result.reference_end,
            result.query_start,
            result.query_end,
        ) == (0, 0, 0, 0)
        assert result.cigar == ""
        assert result.aligned_reference == ""
        assert result.aligned_query == ""

    def test_zero_score_fragment_is_not_returned(self) -> None:
        # Mismatch penalty zero: single differing pair scores 0, never > 0.
        result = align_pair(
            rec("A", "r"), rec("C", "q"),
            mode="local", match_score=1, mismatch_penalty=0,
            gap_open=0, gap_extend=0,
        )
        assert result.score == 0
        assert result.cigar == ""

    def test_tie_prefers_smaller_coordinates(self) -> None:
        # Two equally good matches "AA" at positions 0 and 2.
        result = align_pair(
            rec("AACC", "r"), rec("AA", "q"),
            mode="local", match_score=1, mismatch_penalty=1,
            gap_open=5, gap_extend=2,
        )
        assert result.score == 2
        assert (result.reference_start, result.reference_end) == (0, 2)
        assert (result.query_start, result.query_end) == (0, 2)
        assert result.cigar == "2="

    def test_local_fragment_can_include_internal_gap(self) -> None:
        result = align_pair(
            rec("CCACGTCC", "r"), rec("GGAGTGG", "q"),
            mode="local", match_score=2, mismatch_penalty=3,
            gap_open=1, gap_extend=1,
        )
        # Gapped "ACGT"/"AGT" scores 3*2 - 1 = 5, beating the best
        # ungapped common fragment "GT" (4).
        assert result.score == 5
        assert result.cigar == "1=1D2="
        assert result.aligned_reference == "ACGT"
        assert result.aligned_query == "A-GT"
        assert (result.reference_start, result.reference_end) == (2, 6)
        assert (result.query_start, result.query_end) == (2, 5)


class TestParameters:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"match_score": 0},
            {"match_score": -1},
            {"match_score": True},
            {"match_score": 1.0},
            {"mismatch_penalty": -1},
            {"mismatch_penalty": False},
            {"gap_open": -5},
            {"gap_open": True},
            {"gap_extend": -2},
            {"gap_extend": 2.0},
        ],
    )
    def test_invalid_scoring_parameters_raise_value_error(self, kwargs) -> None:
        with pytest.raises(ValueError):
            align_pair(rec("A", "r"), rec("A", "q"), **kwargs)

    def test_invalid_mode_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            align_pair(rec("A", "r"), rec("A", "q"), mode="semi-global")

    def test_parameters_validated_before_records(self) -> None:
        with pytest.raises(ValueError):
            align_pair(rec("Z", "r"), rec("Z", "q"), match_score=0)

    def test_invalid_base_raises_sequence_validation_error(self) -> None:
        with pytest.raises(SequenceValidationError):
            align_pair(rec("AZ", "r"), rec("A", "q"))
        with pytest.raises(SequenceValidationError):
            align_pair(rec("A", "r"), rec("AZ", "q"))


class TestDeterminism:
    def test_same_inputs_produce_equal_results(self) -> None:
        a = align_pair(rec("ACGTTAC", "r"), rec("AGTAC", "q"), "global", match_score=2, mismatch_penalty=1, gap_open=2, gap_extend=1)
        b = align_pair(rec("ACGTTAC", "r"), rec("AGTAC", "q"), "global", match_score=2, mismatch_penalty=1, gap_open=2, gap_extend=1)
        assert a == b

    def test_result_field_order_is_dataclass_order(self) -> None:
        result = align_pair(rec("A", "r"), rec("A", "q"))
        assert list(result.__dict__.keys()) == [
            "reference",
            "query",
            "mode",
            "score",
            "reference_start",
            "reference_end",
            "query_start",
            "query_end",
            "cigar",
            "aligned_reference",
            "aligned_query",
        ]
