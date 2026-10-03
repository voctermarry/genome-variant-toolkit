"""Tests for genome_variant.deduplication.deduplicate_reads."""

from __future__ import annotations

import pytest

from genome_variant.deduplication import deduplicate_reads
from genome_variant.quality import ReadQualityError
from genome_variant.sequence_io import SequenceRecord, SequenceValidationError


def rec(identifier, sequence, description="", quality=None):
    return SequenceRecord(identifier, sequence, description, quality)


def test_empty_iterable_returns_empty_tuple() -> None:
    assert deduplicate_reads([]) == ()
    assert deduplicate_reads(iter([]), canonical=True) == ()


def test_identical_sequences_keep_earliest_record() -> None:
    records = [
        rec("r1", "ACGT", "first"),
        rec("r2", "TTTT"),
        rec("r3", "ACGT", "second"),
    ]
    assert deduplicate_reads(records) == (records[0], records[1])


def test_same_identifier_different_sequences_not_merged() -> None:
    records = [rec("r1", "ACGT"), rec("r1", "TGCA")]
    assert deduplicate_reads(records) == tuple(records)


def test_result_ordered_by_first_occurrence() -> None:
    records = [
        rec("r1", "GGGG"),
        rec("r2", "AAAA"),
        rec("r3", "GGGG"),
        rec("r4", "CCCC"),
        rec("r5", "AAAA"),
    ]
    result = deduplicate_reads(records)
    assert [r.identifier for r in result] == ["r1", "r2", "r4"]


def test_chunking_does_not_change_result() -> None:
    records = [
        rec("r1", "ACGT"),
        rec("r2", "TTTT"),
        rec("r3", "ACGT"),
        rec("r4", "GGGG"),
        rec("r5", "TTTT"),
    ]
    expected = deduplicate_reads(records)
    for size in range(1, len(records) + 1):
        chunks = (
            records[start : start + size]
            for start in range(0, len(records), size)
        )
        assert deduplicate_reads(r for chunk in chunks for r in chunk) == expected


def test_canonical_groups_reverse_complements() -> None:
    # The reverse complement of AACGTA is TACGTT.
    records = [
        rec("r1", "AACGTA", "forward"),
        rec("r2", "TACGTT", "reverse"),
    ]
    result = deduplicate_reads(records, canonical=True)
    # The representative keeps its original orientation and description.
    assert result == (records[0],)
    # Without canonical the two orientations are distinct groups.
    assert deduplicate_reads(records) == tuple(records)


def test_canonical_uses_full_iupac_reverse_complement() -> None:
    # The reverse complement of ACGTRYSWKMBDHVN is NBDHVKMWSRYACGT.
    forward = rec("r1", "ACGTRYSWKMBDHVN")
    reverse = rec("r2", "NBDHVKMWSRYACGT")
    assert deduplicate_reads([forward, reverse], canonical=True) == (forward,)
    assert deduplicate_reads([forward, reverse]) == (forward, reverse)


def test_canonical_palindrome_groups_with_itself() -> None:
    records = [rec("r1", "ACGT"), rec("r2", "ACGT")]
    assert deduplicate_reads(records, canonical=True) == (records[0],)


def test_fasta_group_representative_is_unchanged() -> None:
    # The grouping key is the reverse complement AACGTA, but the
    # representative keeps its original sequence, orientation and
    # description instead of being rewritten to the key.
    records = [rec("r1", "TACGTT", "some description"), rec("r2", "AACGTA")]
    result = deduplicate_reads(records, canonical=True)
    assert result[0] is records[0]
    assert result[0].sequence == "TACGTT"
    assert result[0].description == "some description"


def test_fastq_group_keeps_highest_quality_sum() -> None:
    low = rec("r1", "ACGT", quality=(10, 10, 10, 10))
    high = rec("r2", "ACGT", quality=(40, 40, 40, 40))
    mid = rec("r3", "ACGT", quality=(20, 20, 20, 20))
    assert deduplicate_reads([low, high, mid]) == (high,)
    assert deduplicate_reads([high, low, mid]) == (high,)


def test_fastq_quality_tie_keeps_earliest() -> None:
    first = rec("r1", "ACGT", quality=(30, 10, 30, 10))
    second = rec("r2", "ACGT", quality=(10, 30, 10, 30))
    assert deduplicate_reads([first, second]) == (first,)
    assert deduplicate_reads([second, first]) == (second,)


def test_fastq_representative_preserves_record() -> None:
    best = rec("r2", "ACGT", "desc", (40, 40, 40, 40))
    result = deduplicate_reads([rec("r1", "ACGT", quality=(1, 1, 1, 1)), best])
    assert result[0] is best
    assert result[0].quality == (40, 40, 40, 40)


def test_canonical_fastq_groups_across_orientations() -> None:
    forward = rec("r1", "AACGTA", quality=(10, 10, 10, 10, 10, 10))
    reverse = rec("r2", "TACGTT", quality=(30, 30, 30, 30, 30, 30))
    result = deduplicate_reads([forward, reverse], canonical=True)
    # Highest quality sum wins; orientation is preserved as-is.
    assert result == (reverse,)


def test_canonical_must_be_boolean() -> None:
    for bad in (0, 1, "true", None, 1.0):
        with pytest.raises(ValueError):
            deduplicate_reads([], canonical=bad)


def test_canonical_checked_before_input_consumed() -> None:
    def exploding():
        raise AssertionError("input must not be consumed")
        yield

    with pytest.raises(ValueError):
        deduplicate_reads(exploding(), canonical=1)


def test_invalid_base_raises_sequence_validation_error() -> None:
    with pytest.raises(SequenceValidationError):
        deduplicate_reads([rec("r1", "ACGX")])
    with pytest.raises(SequenceValidationError):
        deduplicate_reads([rec("r1", "acgt")])


def test_quality_length_mismatch_raises_read_quality_error() -> None:
    with pytest.raises(ReadQualityError):
        deduplicate_reads([rec("r1", "ACGT", quality=(10, 10))])


def test_quality_score_out_of_range_raises_read_quality_error() -> None:
    with pytest.raises(ReadQualityError):
        deduplicate_reads([rec("r1", "ACGT", quality=(10, 10, 10, 94))])
    with pytest.raises(ReadQualityError):
        deduplicate_reads([rec("r1", "ACGT", quality=(10, 10, 10, -1))])


def test_missing_quality_score_raises_read_quality_error() -> None:
    with pytest.raises(ReadQualityError):
        deduplicate_reads([rec("r1", "ACGT", quality=(10, None, 10, 10))])


def test_mixed_quality_and_plain_records_raise_value_error() -> None:
    with pytest.raises(ValueError):
        deduplicate_reads(
            [rec("r1", "ACGT", quality=(10, 10, 10, 10)), rec("r2", "ACGT")]
        )
    with pytest.raises(ValueError):
        deduplicate_reads(
            [rec("r1", "ACGT"), rec("r2", "ACGT", quality=(10, 10, 10, 10))]
        )


def test_returns_tuple() -> None:
    result = deduplicate_reads([rec("r1", "ACGT")])
    assert isinstance(result, tuple)
