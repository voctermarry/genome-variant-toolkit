"""Tests for pairwise alignment (genome_variant.alignment)."""

from __future__ import annotations

import json

import pytest

from genome_variant.alignment import align_pair
from genome_variant.sequence_io import SequenceRecord, SequenceValidationError

KEYS = [
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


def record(sequence: str, identifier: str = "r", quality=None) -> SequenceRecord:
    return SequenceRecord(identifier, sequence, quality=quality)


def align(reference: str, query: str, **kwargs):
    return align_pair(record(reference, "ref"), record(query, "qry"), **kwargs)


class TestGlobal:
    def test_identical_sequences(self) -> None:
        result = align("ACGT", "ACGT")
        assert list(result) == KEYS
        assert result["reference"] == "ref"
        assert result["query"] == "qry"
        assert result["mode"] == "global"
        assert result["score"] == 8
        assert (
            result["reference_start"],
            result["reference_end"],
            result["query_start"],
            result["query_end"],
        ) == (0, 4, 0, 4)
        assert result["cigar"] == "4="
        assert result["aligned_reference"] == "ACGT"
        assert result["aligned_query"] == "ACGT"

    def test_deletion_and_insertion_operations(self) -> None:
        # D consumes the reference: the gap character appears in the query.
        result = align("ACGT", "AGT")
        assert result["score"] == 4
        assert result["cigar"] == "1=1D2="
        assert result["aligned_reference"] == "ACGT"
        assert result["aligned_query"] == "A-GT"
        assert (result["reference_start"], result["reference_end"]) == (0, 4)
        assert (result["query_start"], result["query_end"]) == (0, 3)

        # I consumes the query: the gap character appears in the reference.
        result = align("AGT", "ACGT")
        assert result["cigar"] == "1=1I2="
        assert result["aligned_reference"] == "A-GT"
        assert result["aligned_query"] == "ACGT"
        assert result["reference_end"] == 3
        assert result["query_end"] == 4

    def test_mismatch_operation(self) -> None:
        result = align("ACGT", "ACGA")
        assert result["cigar"] == "3=1X"
        assert result["score"] == 4

    def test_adjacent_operations_are_merged(self) -> None:
        result = align("ACGTACGT", "ACGT", gap_open_penalty=3, gap_extend_penalty=1)
        # One length-4 gap costs open + 3 extends = 6.
        assert result["score"] == 2
        assert result["cigar"] == "1=4D3="
        assert result["aligned_reference"] == "ACGTACGT"
        assert result["aligned_query"] == "A----CGT"

    def test_gap_extend_is_charged_per_extra_column(self) -> None:
        # A length-3 gap costs open + 2 extends = 5 against 3 matches.
        result = align("ACGTAC", "AAC", gap_open_penalty=3, gap_extend_penalty=1)
        assert result["score"] == 1
        assert result["cigar"] == "1=3D2="

        # A length-3 internal gap against 2 matches: 4 - 5 = -1, and the
        # gap is one merged D operation, not three.
        result = align("AAAAA", "AA", gap_open_penalty=3, gap_extend_penalty=1)
        assert result["score"] == -1
        assert result["cigar"] == "1=3D1="

    def test_global_covers_full_sequences_with_leading_and_trailing_gaps(self) -> None:
        result = align("ACGT", "TACGTA", gap_open_penalty=2, gap_extend_penalty=1)
        assert result["reference_start"] == 0
        assert result["query_start"] == 0
        assert result["reference_end"] == 4
        assert result["query_end"] == 6
        assert result["aligned_reference"].replace("-", "") == "ACGT"
        assert result["aligned_query"].replace("-", "") == "TACGTA"
        assert len(result["aligned_reference"]) == len(result["aligned_query"])

    def test_only_identical_characters_count_as_matches(self) -> None:
        # N == N is an exact character match; compatible-but-different
        # IUPAC symbols are mismatches.
        assert align("NN", "NN")["cigar"] == "2="
        result = align("AR", "AA")
        assert result["cigar"] == "1=1X"
        result = align("RS", "AG")
        assert result["cigar"] == "2X"

    def test_empty_sequences(self) -> None:
        result = align("", "")
        assert result["score"] == 0
        assert result["cigar"] == ""
        assert result["aligned_reference"] == ""
        assert result["aligned_query"] == ""
        assert (result["reference_end"], result["query_end"]) == (0, 0)

        result = align("AC", "")
        assert result["cigar"] == "2D"
        assert result["aligned_reference"] == "AC"
        assert result["aligned_query"] == "--"

        result = align("", "AC")
        assert result["cigar"] == "2I"
        assert result["aligned_reference"] == "--"
        assert result["aligned_query"] == "AC"

    def test_fastq_quality_is_ignored(self) -> None:
        reference = record("ACGT", "ref")
        good = SequenceRecord("q", "ACGT", quality=(40, 40, 40, 40))
        bad = SequenceRecord("q", "ACGT", quality=(0, 0, 0, 0))
        assert align_pair(reference, good) == align_pair(reference, bad)


class TestLocal:
    def test_best_segment_coordinates_and_strings(self) -> None:
        result = align("NNACGTNN", "RRACGTRR", mode="local")
        assert result["score"] == 8
        assert result["cigar"] == "4="
        assert result["aligned_reference"] == "ACGT"
        assert result["aligned_query"] == "ACGT"
        assert (
            result["reference_start"],
            result["reference_end"],
            result["query_start"],
            result["query_end"],
        ) == (2, 6, 2, 6)

    def test_no_positive_segment_returns_empty_alignment(self) -> None:
        result = align(
            "AAAA",
            "CCCC",
            mode="local",
            match_score=1,
            mismatch_penalty=1,
            gap_open_penalty=1,
        )
        assert result["score"] == 0
        assert result["cigar"] == ""
        assert result["aligned_reference"] == ""
        assert result["aligned_query"] == ""
        assert (
            result["reference_start"],
            result["reference_end"],
            result["query_start"],
            result["query_end"],
        ) == (0, 0, 0, 0)

    def test_mismatch_zero_penalty_still_has_no_positive_seed_for_mismatch(self) -> None:
        result = align(
            "AAA", "CCC", mode="local", match_score=1, mismatch_penalty=0
        )
        assert result["score"] == 0
        assert result["cigar"] == ""

    def test_earliest_segment_wins_a_score_tie(self) -> None:
        # Two identical 4-base match windows separated by five IUPAC
        # mismatch columns; bridging them scores 16 - 10 = 6, below a
        # single window's 8, so exactly one window is returned and the
        # earlier one wins the tie.
        result = align("NNACGTNNNNNACGT", "RRACGTRRRRRACGT", mode="local")
        assert result["score"] == 8
        assert result["cigar"] == "4="
        assert result["reference_start"] == 2
        assert result["reference_end"] == 6
        assert result["query_start"] == 2
        assert result["query_end"] == 6

    def test_tie_uses_coordinates_then_cigar(self) -> None:
        # Zero gap penalties: leading gaps are free, so a segment may be
        # reported with starts pulled toward 0; coordinates are minimized
        # first.
        result = align(
            "TACG",
            "NACG",
            mode="local",
            mismatch_penalty=0,
            gap_open_penalty=0,
            gap_extend_penalty=0,
        )
        assert result["score"] == 6
        assert result["reference_start"] == 0
        assert result["query_start"] == 0
        assert result["cigar"] == "1D1I3="

        # In global mode the full-length alignment may be expressed as a
        # mismatch column ("1X1=") or as two length-1 gaps flanking the
        # match ("1D1=1I"); the scores and all coordinates are equal, so
        # the lexicographically smaller CIGAR wins.
        result = align(
            "AC",
            "CA",
            mode="global",
            mismatch_penalty=0,
            gap_open_penalty=0,
            gap_extend_penalty=0,
        )
        assert (result["reference_end"], result["query_end"]) == (2, 2)
        assert result["score"] == 2
        assert result["cigar"] == "1D1=1I"

    def test_segment_can_open_with_a_gap(self) -> None:
        result = align(
            "ACGT",
            "AAACGT",
            mode="local",
            mismatch_penalty=0,
            gap_open_penalty=0,
            gap_extend_penalty=0,
        )
        # With free gaps the alignment reaches the earliest possible
        # starts; the match at reference 0 anchors it there and the three
        # extra query bases are an internal insertion.  ("1=2I3=" ties
        # "3I4=" on score and all four coordinates but sorts earlier.)
        assert result["score"] == 8
        assert result["cigar"] == "1=2I3="
        assert result["reference_start"] == 0
        assert result["query_start"] == 0


class TestValidation:
    @pytest.mark.parametrize("value", [0, -1, 1.0, True, False, "2", None])
    def test_match_score_must_be_positive_non_boolean_int(self, value) -> None:
        with pytest.raises(ValueError):
            align_pair(record("A"), record("A"), match_score=value)

    @pytest.mark.parametrize("parameter", ["mismatch_penalty", "gap_open_penalty", "gap_extend_penalty"])
    @pytest.mark.parametrize("value", [-1, 1.0, True, False, 0.0, "0"])
    def test_penalties_must_be_non_negative_non_boolean_ints(self, parameter, value) -> None:
        with pytest.raises(ValueError):
            align_pair(record("A"), record("A"), **{parameter: value})

    def test_zero_penalties_are_accepted(self) -> None:
        result = align_pair(
            record("AC"), record("A"), mismatch_penalty=0,
            gap_open_penalty=0, gap_extend_penalty=0,
        )
        assert result["score"] == 2

    def test_invalid_mode(self) -> None:
        with pytest.raises(ValueError):
            align_pair(record("A"), record("A"), mode="semi-global")

    def test_invalid_base_raises_sequence_validation_error(self) -> None:
        with pytest.raises(SequenceValidationError):
            align_pair(record("ACGZ"), record("ACGT"))
        with pytest.raises(SequenceValidationError):
            align_pair(record("ACGT"), record("ACGZ"))

    def test_parameter_error_is_raised_before_records_are_examined(self) -> None:
        # Even though a record carries an invalid base, the illegal score
        # parameter must surface first as a plain ValueError.
        with pytest.raises(ValueError) as exc_info:
            align_pair(record("Z"), record("Z"), match_score=True)
        assert not isinstance(exc_info.value, SequenceValidationError)


class TestDeterminism:
    def test_same_inputs_produce_byte_identical_json(self) -> None:
        kwargs = dict(gap_open_penalty=3, gap_extend_penalty=1)
        first = json.dumps(align("ACGTACGT", "AGT", **kwargs), separators=(",", ":"))
        second = json.dumps(align("ACGTACGT", "AGT", **kwargs), separators=(",", ":"))
        assert first == second
        assert list(json.loads(first)) == KEYS

    def test_json_is_compact_and_key_order_is_fixed(self) -> None:
        encoded = json.dumps(align("ACGT", "ACGN"), separators=(",", ":"))
        assert '": ' not in encoded
        assert encoded.startswith('{"reference":"ref","query":"qry","mode":"global"')
