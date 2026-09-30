"""Tests for genome_variant.quality."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from genome_variant.quality import ReadQualityError, filter_reads
from genome_variant.sequence_io import SequenceRecord


def rec(
    identifier: str,
    sequence: str,
    quality,
    description: str = "",
) -> SequenceRecord:
    return SequenceRecord(identifier, sequence, description, tuple(quality))


def apply(records, **kwargs) -> list[SequenceRecord]:
    return list(filter_reads(iter(records), **kwargs))


class TestEndTrimming:
    def test_trims_strictly_below_threshold_from_both_ends(self) -> None:
        # Scores 19 and 20 at the ends; 20 is kept, 19 removed.
        record = rec(
            "r1",
            "ACGTACGT",
            (19, 20, 30, 40, 40, 30, 20, 19),
        )
        result = apply([record], min_end_quality=20, min_length=1, min_mean_quality=0)
        assert result == [
            rec("r1", "CGTACG", (20, 30, 40, 40, 30, 20))
        ]

    def test_threshold_bases_are_kept(self) -> None:
        record = rec("r", "AAAA", (0, 20, 20, 0))
        result = apply([record], min_end_quality=20, min_length=1, min_mean_quality=0)
        assert result == [rec("r", "AA", (20, 20))]

    def test_internal_low_quality_bases_are_kept(self) -> None:
        record = rec("r", "ACGT", (40, 2, 3, 40))
        result = apply([record], min_end_quality=20, min_length=1, min_mean_quality=0)
        assert result == [rec("r", "ACGT", (40, 2, 3, 40))]

    def test_zero_threshold_trims_nothing(self) -> None:
        record = rec("r", "ACGT", (0, 0, 0, 0))
        result = apply([record], min_end_quality=0, min_length=1, min_mean_quality=0)
        assert result == [record]

    def test_all_bases_trimmed_drops_record(self) -> None:
        record = rec("r", "ACGT", (5, 6, 7, 8))
        assert apply([record], min_end_quality=20, min_length=1) == []

    def test_runs_of_low_quality_ends(self) -> None:
        record = rec("r", "AAAAACGTAAAA", (1, 1, 1, 1, 1, 40, 40, 40, 19, 18, 17, 16))
        result = apply([record], min_end_quality=20, min_length=1, min_mean_quality=0)
        assert result == [rec("r", "CGT", (40, 40, 40))]


class TestFiltering:
    def test_too_short_after_trimming_is_dropped(self) -> None:
        record = rec("r", "ACGTACGT", (40, 40, 40, 40, 40, 40, 5, 5))
        assert apply([record], min_end_quality=20, min_length=7) == []
        # Length equal to the threshold survives.
        result = apply([record], min_end_quality=20, min_length=6, min_mean_quality=0)
        assert len(result) == 1
        assert len(result[0].sequence) == 6

    def test_mean_quality_filter_uses_trimmed_segment(self) -> None:
        # End trimming removes the trailing 5; the remaining segment
        # (19, 20) has sum 39 < 20 * 2 and the read is dropped.
        record = rec("r", "ACG", (19, 20, 5))
        assert apply(
            [record], min_end_quality=10, min_length=1, min_mean_quality=20
        ) == []
        # A higher-quality tail survives the same mean threshold.
        kept = rec("r2", "ACG", (19, 20, 21))
        result = apply(
            [kept], min_end_quality=10, min_length=1, min_mean_quality=20
        )
        assert len(result) == 1

    def test_mean_equal_to_threshold_is_kept(self) -> None:
        record = rec("r", "ACG", (20, 20, 20))
        result = apply([record], min_end_quality=20, min_length=3, min_mean_quality=20)
        assert len(result) == 1

    def test_mean_just_below_threshold_is_dropped(self) -> None:
        record = rec("r", "ACG", (21, 20, 18))
        assert apply(
            [record], min_end_quality=0, min_length=1, min_mean_quality=20
        ) == []

    def test_surviving_record_keeps_identifier_and_description(self) -> None:
        record = rec(
            "read/1",
            "NACGTN",
            (5, 30, 30, 30, 30, 5),
            description="some description",
        )
        result = apply(
            [record], min_end_quality=20, min_mean_quality=0, min_length=1
        )
        assert result[0].identifier == "read/1"
        assert result[0].description == "some description"
        assert result[0].sequence == "ACGT"
        assert result[0].quality == (30, 30, 30, 30)

    def test_order_is_preserved_and_records_are_not_merged(self) -> None:
        records = [
            rec("a", "A" * 30, (40,) * 30),
            rec("b", "A" * 30, (1,) * 30),
            rec("a", "C" * 30, (40,) * 30),
        ]
        result = apply(records)
        assert [r.identifier for r in result] == ["a", "a"]
        assert result[0].sequence == "A" * 30
        assert result[1].sequence == "C" * 30

    def test_all_filtered_produces_empty_output_iterable(self) -> None:
        records = [rec("a", "ACGT", (1, 1, 1, 1)), rec("b", "ACGT", (2, 2, 2, 2))]
        assert apply(records, min_end_quality=20, min_length=1) == []


class TestRecordErrors:
    def test_missing_quality_raises_with_identifier(self) -> None:
        records = [SequenceRecord("fasta-record", "ACGT")]
        with pytest.raises(ReadQualityError) as info:
            apply(records)
        assert "fasta-record" in str(info.value)

    def test_quality_length_mismatch_raises_with_identifier(self) -> None:
        records = [rec("x", "ACGT", (30, 30))]
        with pytest.raises(ReadQualityError) as info:
            apply(records)
        assert "x" in str(info.value)

    @pytest.mark.parametrize("score", [-1, 94, 100])
    def test_quality_out_of_range_raises(self, score: int) -> None:
        records = [rec("bad", "ACG", (30, score, 30))]
        with pytest.raises(ReadQualityError) as info:
            apply(records)
        assert "bad" in str(info.value)

    def test_boundary_scores_0_and_93_are_accepted(self) -> None:
        record = rec("edge", "AC", (0, 93))
        result = apply(
            [record], min_end_quality=0, min_length=2, min_mean_quality=0
        )
        assert result == [record]

    def test_error_is_lazy_until_record_reached(self) -> None:
        good = rec("good", "A" * 30, (40,) * 30)
        bad = SequenceRecord("bad", "ACGT", None)
        iterator: Iterator[SequenceRecord] = filter_reads(iter([good, bad]))
        assert next(iterator).identifier == "good"
        with pytest.raises(ReadQualityError):
            next(iterator)


class TestParameterValidation:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"min_end_quality": -1},
            {"min_end_quality": 94},
            {"min_mean_quality": -1},
            {"min_mean_quality": 94},
            {"min_length": 0},
            {"min_length": -5},
        ],
    )
    def test_out_of_range_thresholds_raise_value_error(self, kwargs) -> None:
        with pytest.raises(ValueError):
            filter_reads(iter([]), **kwargs)

    @pytest.mark.parametrize(
        "value",
        [True, False, 1.0, 20.5, "20", None],
    )
    def test_non_integer_or_boolean_thresholds_raise_value_error(self, value) -> None:
        for name in ("min_end_quality", "min_mean_quality", "min_length"):
            with pytest.raises(ValueError):
                filter_reads(iter([]), **{name: value})

    def test_boundary_thresholds_accepted(self) -> None:
        # Should not raise; iteration is empty.
        list(
            filter_reads(
                iter([]),
                min_end_quality=0,
                min_mean_quality=93,
                min_length=1,
            )
        )

    def test_invalid_parameter_raises_before_iteration(self) -> None:
        def records() -> Iterator[SequenceRecord]:
            raise AssertionError("input must not be consumed")
            yield  # pragma: no cover

        with pytest.raises(ValueError):
            filter_reads(records(), min_length=0)

    def test_read_quality_error_is_value_error(self) -> None:
        assert issubclass(ReadQualityError, ValueError)
