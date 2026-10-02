"""Tests for genome_variant.coverage.coverage_report."""

from __future__ import annotations

import pytest

from genome_variant.coverage import CoverageStats, coverage_report
from genome_variant.mapping import MappingReferenceError
from genome_variant.sequence_io import (
    SequenceRecord,
    SequenceValidationError,
)


def rec(identifier, sequence, quality=None):
    return SequenceRecord(identifier, sequence, quality=quality)


def stats_dict(
    record,
    id_,
    length,
    mapped,
    observed,
    covered,
    coverage,
    mean_depth,
    concordant,
    concordance,
):
    return {
        "reference_record": record,
        "reference": id_,
        "length": length,
        "mapped_reads": mapped,
        "observed_bases": observed,
        "covered_bases": covered,
        "coverage_fraction": coverage,
        "mean_depth": mean_depth,
        "concordant_bases": concordant,
        "concordance_fraction": concordance,
    }


EXPECTED_KEYS = [
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


def test_basic_per_reference_stats() -> None:
    refs = [rec("c1", "ACGTACGT"), rec("c2", "TTTTGGGG")]
    # b carries an interior mismatch: 4=1X3= (8 observations, 7 concordant).
    reads = [rec("a", "ACGT"), rec("b", "ACGTTCGT")]
    result = list(coverage_report(refs, reads))
    assert [list(item.to_dict()) for item in result] == [EXPECTED_KEYS] * 2
    first = result[0].to_dict()
    assert first == stats_dict(
        0, "c1", 8, 2, 12, 8, "1.000000", "1.500000", 11, "0.916667"
    )
    second = result[1].to_dict()
    assert second == stats_dict(
        1, "c2", 8, 0, 0, 0, "0.000000", "0.000000", 0, "0.000000"
    )


def test_unmapped_reads_are_excluded() -> None:
    # Neither the read nor its reverse complement (also GCGC) shares a
    # base with the T-only reference, so no positive fragment exists.
    refs = [rec("c1", "TTTT")]
    stats = list(coverage_report(refs, [rec("lost", "GCGC")]))[0]
    assert stats.mapped_reads == 0
    assert stats.observed_bases == 0
    assert stats.covered_bases == 0
    assert stats.concordant_bases == 0


def test_one_line_per_reference_even_without_hits() -> None:
    # The read's A/T bases never occur in the G/C-only references, so no
    # positive fragment exists against any record.
    refs = [rec("c1", "GGGG"), rec("c2", "CCCC"), rec("c3", "GGCC")]
    result = list(coverage_report(refs, [rec("lost", "TATA")]))
    assert [item.reference_record for item in result] == [0, 1, 2]
    assert all(item.mapped_reads == 0 for item in result)


def test_duplicate_reference_identifiers_not_merged() -> None:
    refs = [rec("x", "ACGT"), rec("x", "ACGT")]
    result = list(coverage_report(refs, [rec("a", "ACGT")]))
    assert len(result) == 2
    assert result[0].mapped_reads == 1
    assert result[1].mapped_reads == 0


def test_depth_accumulates_across_reads() -> None:
    refs = [rec("c1", "ACGT")]
    reads = [rec("a", "ACGT"), rec("b", "ACGT"), rec("c", "ACGT")]
    stats = list(coverage_report(refs, reads))[0]
    assert stats.observed_bases == 12
    assert stats.covered_bases == 4
    assert stats.concordant_bases == 12
    payload = stats.to_dict()
    assert payload["coverage_fraction"] == "1.000000"
    assert payload["mean_depth"] == "3.000000"
    assert payload["concordance_fraction"] == "1.000000"


def test_reverse_strand_read_counts_same() -> None:
    # The read's reverse complement AAAAACGT matches the reference over its
    # full length, so the winning hit is on the "-" strand with 8 =.
    refs = [rec("c1", "AAAAACGT")]
    forward = list(coverage_report(refs, [rec("f", "AAAAACGT")]))[0]
    reverse = list(coverage_report(refs, [rec("r", "ACGTTTTT")]))[0]
    assert reverse == forward


def test_deletion_adds_no_depth() -> None:
    # ACGACGT aligns to ACGTACGT as 3=1D4= with free gaps.
    refs = [rec("c1", "ACGTACGT")]
    stats = list(
        coverage_report(
            refs, [rec("del", "ACGACGT")], gap_open=0, gap_extend=0
        )
    )[0]
    assert stats.observed_bases == 7
    assert stats.covered_bases == 7
    assert stats.concordant_bases == 7


def test_insertion_has_no_reference_position() -> None:
    # Cigar 1=1D1I1=1D1= spanning reference positions 0, 2 and 4.
    refs = [rec("c1", "ACGTACGT")]
    stats = list(
        coverage_report(
            refs, [rec("ins", "AAGAA")], gap_open=0, gap_extend=0
        )
    )[0]
    assert stats.observed_bases == 3
    assert stats.covered_bases == 3
    assert stats.concordant_bases == 3


def test_iupac_verbatim_equality_is_concordant() -> None:
    # R == R is a match, even though R is an ambiguity code.
    refs = [rec("c1", "ARGT")]
    stats = list(coverage_report(refs, [rec("a", "ARGT")]))[0]
    assert stats.observed_bases == 4
    assert stats.concordant_bases == 4
    assert stats.to_dict()["concordance_fraction"] == "1.000000"


def test_iupac_ambiguity_pair_is_discordant() -> None:
    # An interior R versus T are different characters: one X column adds
    # depth but no concordance; the flanking matches keep the fragment
    # positive.
    refs = [rec("c1", "ACGTACGT")]
    stats = list(coverage_report(refs, [rec("a", "ACGRACGT")]))[0]
    assert stats.observed_bases == 8
    assert stats.covered_bases == 8
    assert stats.concordant_bases == 7
    assert stats.to_dict()["concordance_fraction"] == f"{7/8:.6f}"


def test_fastq_quality_values_ignored() -> None:
    refs = [rec("c1", "ACGT")]
    fasta = list(coverage_report(refs, [rec("a", "ACGT")]))[0]
    fastq = list(
        coverage_report(refs, [rec("a", "ACGT", quality=(0, 0, 0, 0))])
    )[0]
    assert fasta == fastq


def test_partial_overlap_covers_only_aligned_columns() -> None:
    refs = [rec("c1", "TTTTACGTTTT")]
    stats = list(coverage_report(refs, [rec("a", "ACGT")]))[0]
    assert stats.observed_bases == 4
    assert stats.covered_bases == 4
    assert stats.length == 11
    payload = stats.to_dict()
    assert payload["coverage_fraction"] == f"{4/11:.6f}"
    assert payload["mean_depth"] == f"{4/11:.6f}"


def test_zero_denominators_use_zero() -> None:
    stats = CoverageStats(0, "c1", 0, 0, 0, 0, 0)
    payload = stats.to_dict()
    assert payload["coverage_fraction"] == "0.000000"
    assert payload["mean_depth"] == "0.000000"
    assert payload["concordance_fraction"] == "0.000000"


def test_min_score_filters_reads() -> None:
    refs = [rec("c1", "ACGTACGT")]
    stats = list(
        coverage_report(refs, [rec("a", "ACGT")], min_score=9)
    )[0]
    assert stats.mapped_reads == 0
    assert stats.observed_bases == 0


def test_empty_reference_raises_before_reads_consumed() -> None:
    def reads():
        raise AssertionError("reads must not be consumed")
        yield  # pragma: no cover

    with pytest.raises(MappingReferenceError):
        list(coverage_report([], reads()))


def test_invalid_parameters_raise_value_error() -> None:
    refs = [rec("c1", "ACGT")]
    reads = [rec("a", "ACGT")]
    with pytest.raises(ValueError):
        list(coverage_report(refs, reads, match_score=0))
    with pytest.raises(ValueError):
        list(coverage_report(refs, reads, mismatch_penalty=-1))
    with pytest.raises(ValueError):
        list(coverage_report(refs, reads, min_score=0))
    with pytest.raises(ValueError):
        list(coverage_report(refs, reads, min_score=True))


def test_invalid_read_base_raises_validation_error() -> None:
    refs = [rec("c1", "ACGT")]
    with pytest.raises(SequenceValidationError):
        list(coverage_report(refs, [rec("a", "AZ")]))


def test_batching_does_not_change_results() -> None:
    refs = [rec("c1", "ACGTACGT"), rec("c2", "TTTTGGGG")]
    reads = [
        rec("a", "ACGT"),
        rec("b", "AAAA"),
        rec("c", "ACGTACGT"),
        rec("d", "GCGC"),
    ]
    once = [item.to_dict() for item in coverage_report(refs, iter(reads))]

    class Batched:
        def __init__(self, items):
            self.items = list(items)

        def __iter__(self):
            for item in self.items:
                yield item

    batched = [
        item.to_dict() for item in coverage_report(refs, Batched(reads))
    ]
    assert once == batched


def test_results_are_lazy_over_reads() -> None:
    refs = [rec("c1", "ACGT")]
    consumed = []

    def reads():
        consumed.append(True)
        yield rec("a", "ACGT")

    iterator = coverage_report(refs, reads())
    assert consumed == []
    list(iterator)
    assert consumed == [True]
