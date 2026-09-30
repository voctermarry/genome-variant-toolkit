"""Tests for genome_variant.quality read trimming and filtering."""

from __future__ import annotations

import pytest

from genome_variant.quality import ReadQualityError, trim_and_filter_reads
from genome_variant.sequence_io import SequenceRecord


def record(identifier, sequence, quality, description=""):
    return SequenceRecord(identifier, sequence, description, tuple(quality))


def test_end_trimming_drops_bases_below_threshold() -> None:
    # Phred 10 at both ends, 20 in the middle; threshold 20 keeps the core.
    rec = record("r1", "AACCGG", [10, 10, 20, 20, 10, 10])
    (kept,) = trim_and_filter_reads(
        [rec], min_end_quality=20, min_mean_quality=0, min_length=1
    )
    assert kept.sequence == "CC"
    assert kept.quality == (20, 20)
    assert kept.identifier == "r1"


def test_quality_equal_to_threshold_is_kept() -> None:
    rec = record("r1", "ACGT", [20, 20, 20, 20])
    (kept,) = trim_and_filter_reads(
        [rec], min_end_quality=20, min_mean_quality=20, min_length=1
    )
    assert kept.sequence == "ACGT"


def test_internal_low_quality_bases_do_not_split() -> None:
    rec = record("r1", "ACGT", [30, 5, 5, 30])
    (kept,) = trim_and_filter_reads(
        [rec], min_end_quality=20, min_mean_quality=0, min_length=1
    )
    assert kept.sequence == "ACGT"
    assert kept.quality == (30, 5, 5, 30)


def test_all_bases_trimmed_drops_record() -> None:
    rec = record("r1", "ACGT", [1, 2, 3, 4])
    assert list(trim_and_filter_reads([rec], min_end_quality=20, min_length=1)) == []


def test_short_after_trimming_drops_record() -> None:
    rec = record("r1", "AACCGGTT", [5, 5, 40, 40, 40, 40, 5, 5])
    assert (
        list(
            trim_and_filter_reads(
                [rec], min_end_quality=20, min_mean_quality=0, min_length=5
            )
        )
        == []
    )


def test_mean_quality_boundary() -> None:
    # Sum 60 over 4 bases: mean 15 < 20 drops, == 15 keeps.
    rec = record("r1", "ACGT", [10, 10, 20, 20])
    assert (
        list(
            trim_and_filter_reads(
                [rec], min_end_quality=0, min_mean_quality=20, min_length=1
            )
        )
        == []
    )
    (kept,) = trim_and_filter_reads(
        [rec], min_end_quality=0, min_mean_quality=15, min_length=1
    )
    assert kept.sequence == "ACGT"


def test_identifier_and_description_preserved_and_order_kept() -> None:
    records = [
        record("b", "ACGT", [40] * 4, "second read"),
        record("a", "TTTT", [1] * 4),
        record("c", "GG", [40, 40]),
    ]
    kept = list(
        trim_and_filter_reads(
            records, min_end_quality=20, min_mean_quality=20, min_length=1
        )
    )
    assert [(r.identifier, r.description) for r in kept] == [
        ("b", "second read"),
        ("c", ""),
    ]


def test_defaults_are_20_20_30() -> None:
    good = record("good", "A" * 30, [20] * 30)
    short = record("short", "A" * 29, [40] * 29)
    kept = list(trim_and_filter_reads([good, short]))
    assert [r.identifier for r in kept] == ["good"]


def test_missing_quality_raises_with_identifier() -> None:
    rec = SequenceRecord("noid", "ACGT")
    with pytest.raises(ReadQualityError, match="noid"):
        list(trim_and_filter_reads([rec], min_length=1))


def test_quality_length_mismatch_raises() -> None:
    rec = SequenceRecord("bad", "ACGT", "", (10, 10))
    with pytest.raises(ReadQualityError, match="bad"):
        list(trim_and_filter_reads([rec], min_length=1))


def test_quality_out_of_range_raises() -> None:
    rec = SequenceRecord("hot", "ACGT", "", (10, 94, 10, 10))
    with pytest.raises(ReadQualityError, match="hot"):
        list(trim_and_filter_reads([rec], min_length=1))


def test_record_errors_are_lazy() -> None:
    good = record("good", "ACGT", [40] * 4)
    bad = SequenceRecord("bad", "ACGT")
    iterator = trim_and_filter_reads([good, bad], min_length=1)
    assert next(iterator).identifier == "good"
    with pytest.raises(ReadQualityError):
        next(iterator)


@pytest.mark.parametrize("value", [True, 1.5, "20", None])
def test_non_integer_thresholds_rejected(value) -> None:
    with pytest.raises(ValueError):
        trim_and_filter_reads([], min_end_quality=value)
    with pytest.raises(ValueError):
        trim_and_filter_reads([], min_mean_quality=value)
    with pytest.raises(ValueError):
        trim_and_filter_reads([], min_length=value)


@pytest.mark.parametrize("value", [-1, 94, 100])
def test_quality_threshold_range(value) -> None:
    with pytest.raises(ValueError):
        trim_and_filter_reads([], min_end_quality=value)
    with pytest.raises(ValueError):
        trim_and_filter_reads([], min_mean_quality=value)


@pytest.mark.parametrize("value", [0, -5])
def test_min_length_at_least_one(value) -> None:
    with pytest.raises(ValueError):
        trim_and_filter_reads([], min_length=value)


def test_threshold_validation_is_eager() -> None:
    def explode():
        raise AssertionError("records must not be touched")
        yield

    with pytest.raises(ValueError):
        trim_and_filter_reads(explode(), min_length=0)
