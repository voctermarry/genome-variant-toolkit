"""Tests for genome_variant.mapping.map_reads."""

from __future__ import annotations

import pytest

from genome_variant.mapping import (
    MappingReferenceError,
    ReadMapping,
    map_reads,
)
from genome_variant.sequence_io import SequenceRecord, SequenceValidationError


def rec(identifier: str, sequence: str) -> SequenceRecord:
    return SequenceRecord(identifier, sequence)


class TestMapReadsBasic:
    def test_forward_mapping(self) -> None:
        (result,) = map_reads([rec("r", "TTACGTAA")], [rec("q", "ACGT")])
        assert result == ReadMapping(
            record=0,
            id="q",
            mapped=True,
            reference_record=0,
            reference="r",
            reference_start=2,
            reference_end=6,
            query_start=0,
            query_end=4,
            strand="+",
            score=8,
            cigar="4=",
        )

    def test_reverse_strand_coordinates_converted(self) -> None:
        # reverse_complement("GGACAA") == "TTGTCC"; its "GTCC" fragment
        # aligns to reference [2, 6) at reverse-read [2, 6), which is
        # [0, 4) on the original read.
        (result,) = map_reads([rec("r", "AAGTCCAA")], [rec("q", "GGACAA")])
        assert result.mapped
        assert result.strand == "-"
        assert (result.reference_start, result.reference_end) == (2, 6)
        assert (result.query_start, result.query_end) == (0, 4)
        assert result.score == 8
        assert result.cigar == "4="

    def test_forward_strand_preferred_on_tie(self) -> None:
        # "ACGT" is a reverse-complement palindrome fragment here: both
        # strands give the same score and coordinates, "+" wins.
        (result,) = map_reads([rec("r", "ACGTACGT")], [rec("q", "ACGT")])
        assert result.mapped
        assert result.strand == "+"
        assert (result.reference_start, result.reference_end) == (0, 4)

    def test_reference_index_breaks_tie(self) -> None:
        (result,) = map_reads(
            [rec("r2", "ACGT"), rec("r1", "ACGT")], [rec("q", "ACGT")]
        )
        assert result.reference_record == 0
        assert result.reference == "r2"

    def test_duplicate_identifiers_not_merged(self) -> None:
        results = list(
            map_reads([rec("r", "ACGT")], [rec("q", "ACGT"), rec("q", "ACGT")])
        )
        assert [r.record for r in results] == [0, 1]
        assert all(r.id == "q" for r in results)
        assert all(r.mapped for r in results)

    def test_unmapped_read_still_reported(self) -> None:
        (result,) = map_reads([rec("r", "CCCC")], [rec("q", "AAAA")])
        assert result == ReadMapping(
            record=0,
            id="q",
            mapped=False,
            reference_record=None,
            reference=None,
            reference_start=None,
            reference_end=None,
            query_start=None,
            query_end=None,
            strand=None,
            score=None,
            cigar=None,
        )

    def test_min_score_filters_candidates(self) -> None:
        references = [rec("r", "AC")]
        reads = [rec("q", "AC")]
        (mapped,) = map_reads(references, reads, min_score=4)
        assert mapped.mapped and mapped.score == 4
        (unmapped,) = map_reads(references, reads, min_score=5)
        assert not unmapped.mapped

    def test_iupac_reverse_complement(self) -> None:
        # reverse_complement("RYSW") == "WSRY"; it aligns to the reference.
        (result,) = map_reads([rec("r", "TTWSRYTT")], [rec("q", "RYSW")])
        assert result.mapped
        assert result.strand == "-"
        assert (result.reference_start, result.reference_end) == (2, 6)

    def test_quality_ignored(self) -> None:
        read = SequenceRecord("q", "ACGT", quality=(0, 0, 0, 0))
        (result,) = map_reads([rec("r", "ACGT")], [read])
        assert result.mapped and result.score == 8

    def test_record_order_and_indices(self) -> None:
        reads = [rec("a", "AAAA"), rec("b", "CCGG"), rec("c", "TTTT")]
        results = list(map_reads([rec("r", "CCGG")], reads))
        assert [r.record for r in results] == [0, 1, 2]
        assert [r.id for r in results] == ["a", "b", "c"]
        assert [r.mapped for r in results] == [False, True, False]

    def test_batching_does_not_change_results(self) -> None:
        references = [rec("r1", "AAGTCCAA"), rec("r2", "ACGTACGT")]
        reads = [rec("q1", "GGACAA"), rec("q2", "ACGT"), rec("q3", "AAAA")]
        whole = list(map_reads(references, reads))
        assert whole == list(map_reads(iter(references), iter(reads)))
        assert whole == list(
            map_reads((r for r in references), (q for q in reads))
        )


class TestMapReadsErrors:
    def test_empty_reference_set(self) -> None:
        with pytest.raises(MappingReferenceError):
            list(map_reads([], [rec("q", "ACGT")]))

    def test_empty_reference_set_is_value_error(self) -> None:
        with pytest.raises(ValueError):
            list(map_reads(iter([]), iter([rec("q", "ACGT")])))

    def test_invalid_parameters_raise_before_consuming_reads(self) -> None:
        def exploding_reads():
            raise AssertionError("reads must not be consumed")
            yield

        references = [rec("r", "ACGT")]
        for kwargs in [
            {"min_score": 0},
            {"min_score": -1},
            {"min_score": True},
            {"min_score": 1.5},
            {"min_score": "2"},
            {"match_score": 0},
            {"match_score": True},
            {"mismatch_penalty": -1},
            {"gap_open": -1},
            {"gap_extend": -1},
        ]:
            with pytest.raises(ValueError):
                map_reads(references, exploding_reads(), **kwargs)

    def test_invalid_reference_base(self) -> None:
        with pytest.raises(SequenceValidationError):
            list(map_reads([rec("r", "AZ")], [rec("q", "ACGT")]))

    def test_invalid_read_base(self) -> None:
        with pytest.raises(SequenceValidationError):
            list(map_reads([rec("r", "ACGT")], [rec("q", "AZ")]))
