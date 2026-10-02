"""Tests for genome_variant.coverage.coverage_report."""

from __future__ import annotations

import pytest

from genome_variant.coverage import ReferenceCoverage, coverage_report
from genome_variant.mapping import MappingReferenceError
from genome_variant.sequence_io import (
    SequenceRecord,
    SequenceValidationError,
)


def rec(identifier, sequence, quality=None):
    return SequenceRecord(identifier, sequence, quality=quality)


class TestBasicCoverage:
    def test_full_match_single_reference(self) -> None:
        result = list(coverage_report([rec("c1", "ACGTACGT")], [rec("r1", "ACGT")]))
        assert result == [
            ReferenceCoverage(
                reference_record=0,
                reference="c1",
                length=8,
                mapped_reads=1,
                observed_bases=4,
                covered_bases=4,
                coverage_fraction="0.500000",
                mean_depth="0.500000",
                concordant_bases=4,
                concordance_fraction="1.000000",
            )
        ]

    def test_to_dict_key_order(self) -> None:
        (entry,) = coverage_report([rec("c1", "ACGT")], [rec("r1", "ACGT")])
        assert list(entry.to_dict()) == [
            "reference_record",
            "reference",
            "length",
            "mapped_reads",
            "observed_bases",
            "covered_bases",
            "coverage_fraction",
            "mean_depth",
            "concordant_bases",
            "concordance_fraction",
        ]

    def test_no_mapped_reads_still_reports(self) -> None:
        (entry,) = coverage_report([rec("c1", "AAAAAAAA")], [rec("r1", "CCCC")])
        assert entry.mapped_reads == 0
        assert entry.observed_bases == 0
        assert entry.covered_bases == 0
        assert entry.coverage_fraction == "0.000000"
        assert entry.mean_depth == "0.000000"
        assert entry.concordant_bases == 0
        assert entry.concordance_fraction == "0.000000"

    def test_empty_reads_reports_every_reference(self) -> None:
        result = list(
            coverage_report([rec("c1", "ACGT"), rec("c2", "TTTT")], [])
        )
        assert [entry.reference for entry in result] == ["c1", "c2"]
        assert all(entry.mapped_reads == 0 for entry in result)

    def test_duplicate_identifiers_not_merged(self) -> None:
        result = list(
            coverage_report(
                [rec("dup", "ACGTACGT"), rec("dup", "ACGTACGT")],
                [rec("r1", "ACGT")],
            )
        )
        assert len(result) == 2
        assert [entry.reference_record for entry in result] == [0, 1]
        assert [entry.reference for entry in result] == ["dup", "dup"]
        # The tie is resolved to the first reference record.
        assert result[0].mapped_reads == 1
        assert result[1].mapped_reads == 0

    def test_depth_stacks_over_reads(self) -> None:
        result = list(
            coverage_report(
                [rec("c1", "ACGT")],
                [rec("r1", "ACGT"), rec("r2", "ACGT")],
            )
        )
        (entry,) = result
        assert entry.mapped_reads == 2
        assert entry.observed_bases == 8
        assert entry.covered_bases == 4
        assert entry.coverage_fraction == "1.000000"
        assert entry.mean_depth == "2.000000"
        assert entry.concordant_bases == 8

    def test_quality_values_ignored(self) -> None:
        with_quality = list(
            coverage_report(
                [rec("c1", "ACGT")], [rec("r1", "ACGT", quality=(0, 5, 40, 93))]
            )
        )
        without_quality = list(
            coverage_report([rec("c1", "ACGT")], [rec("r1", "ACGT")])
        )
        assert with_quality == without_quality


class TestCigarSemantics:
    def test_mismatch_counts_depth_but_not_concordance(self) -> None:
        # 4=1X3=: the mismatched column adds depth but no concordance.
        (entry,) = coverage_report(
            [rec("c1", "ACGTTCGT")], [rec("r1", "ACGTACGT")]
        )
        assert entry.mapped_reads == 1
        assert entry.observed_bases == 8
        assert entry.covered_bases == 8
        assert entry.concordant_bases == 7
        assert entry.concordance_fraction == "0.875000"

    def test_deletion_adds_no_depth(self) -> None:
        # 4=1D3=: the deleted reference position keeps zero depth.
        (entry,) = coverage_report(
            [rec("c1", "ACGTACGT")], [rec("r1", "ACGTCGT")]
        )
        assert entry.observed_bases == 7
        assert entry.covered_bases == 7
        assert entry.coverage_fraction == "0.875000"
        assert entry.mean_depth == "0.875000"
        assert entry.concordant_bases == 7

    def test_insertion_has_no_reference_position(self) -> None:
        # 4=1I4=: the inserted base touches no reference position.
        (entry,) = coverage_report(
            [rec("c1", "ACGTACGT")], [rec("r1", "ACGTTACGT")]
        )
        assert entry.observed_bases == 8
        assert entry.covered_bases == 8
        assert entry.coverage_fraction == "1.000000"
        assert entry.concordant_bases == 8

    def test_reverse_strand_hit_counts(self) -> None:
        # The read's reverse complement matches the reference fully.
        (entry,) = coverage_report(
            [rec("c1", "AAAAGGGG")], [rec("r1", "CCCCTTTT")]
        )
        assert entry.mapped_reads == 1
        assert entry.observed_bases == 8
        assert entry.covered_bases == 8
        assert entry.concordant_bases == 8

    def test_iupac_identity_is_literal(self) -> None:
        # N against N is concordant; R against A is a literal mismatch.
        (entry,) = coverage_report(
            [rec("c1", "ACRTAC")], [rec("r1", "ACNTAC")]
        )
        assert entry.observed_bases == 6
        assert entry.concordant_bases == 5
        assert entry.concordance_fraction == "0.833333"


class TestAdjudication:
    def test_only_the_winning_mapping_counts(self) -> None:
        # The read matches both records fully; the tie goes to the first.
        result = list(
            coverage_report(
                [rec("c1", "ACGT"), rec("c2", "ACGT")], [rec("r1", "ACGT")]
            )
        )
        assert result[0].mapped_reads == 1
        assert result[1].mapped_reads == 1 - 1

    def test_min_score_leaves_read_unmapped(self) -> None:
        (entry,) = coverage_report(
            [rec("c1", "ACGTACGT")],
            [rec("r1", "ACGT")],
            min_score=9,
        )
        assert entry.mapped_reads == 0
        assert entry.observed_bases == 0

    def test_scoring_parameters_change_outcome(self) -> None:
        # With a mismatch penalty of zero the 4=1X3= alignment still wins,
        # but a huge gap open penalty makes the deletion alignment lose to
        # a pure-match fragment.
        (gappy,) = coverage_report(
            [rec("c1", "ACGTACGT")], [rec("r1", "ACGTCGT")]
        )
        assert gappy.observed_bases == 7
        (strict,) = coverage_report(
            [rec("c1", "ACGTACGT")],
            [rec("r1", "ACGTCGT")],
            gap_open=50,
        )
        assert strict.observed_bases == 4


class TestErrors:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"match_score": 0},
            {"match_score": True},
            {"mismatch_penalty": -1},
            {"gap_open": -1},
            {"gap_extend": -1},
            {"min_score": 0},
            {"min_score": 1.5},
            {"min_score": True},
        ],
    )
    def test_invalid_parameters_raise_value_error(self, kwargs) -> None:
        with pytest.raises(ValueError):
            coverage_report([rec("c1", "ACGT")], [rec("r1", "ACGT")], **kwargs)

    def test_invalid_parameters_raise_before_consuming_inputs(self) -> None:
        def exploding():
            raise AssertionError("input consumed")
            yield

        with pytest.raises(ValueError):
            coverage_report(exploding(), exploding(), min_score=0)

    def test_empty_reference_raises(self) -> None:
        with pytest.raises(MappingReferenceError):
            coverage_report([], [rec("r1", "ACGT")])

    def test_empty_reference_raises_before_consuming_reads(self) -> None:
        def exploding():
            raise AssertionError("reads consumed")
            yield

        with pytest.raises(MappingReferenceError):
            coverage_report([], exploding())

    def test_invalid_reference_symbol_raises_on_iteration(self) -> None:
        report = coverage_report([rec("c1", "ACGT-")], [rec("r1", "ACGT")])
        with pytest.raises(SequenceValidationError):
            list(report)

    def test_invalid_read_symbol_raises_on_iteration(self) -> None:
        report = coverage_report([rec("c1", "ACGT")], [rec("r1", "ACGT-")])
        with pytest.raises(SequenceValidationError):
            list(report)


class TestBatching:
    def test_batching_does_not_change_results(self) -> None:
        references = [rec("c1", "ACGTACGT"), rec("c2", "TTTTACGT")]
        reads = [
            rec("r1", "ACGT"),
            rec("r2", "TTTT"),
            rec("r3", "ACGTCGT"),
            rec("r4", "GGGG"),
        ]
        whole = list(coverage_report(references, reads))
        one_by_one = list(
            coverage_report(iter(references), iter(reads))
        )
        assert whole == one_by_one
