"""Tests for genome_variant.kmer.build_kmer_index."""

from __future__ import annotations

import pytest

from genome_variant.kmer import KmerOccurrence, build_kmer_index
from genome_variant.sequence_io import SequenceRecord, SequenceValidationError


def rec(identifier, sequence, quality=None):
    return SequenceRecord(identifier, sequence, quality=quality)


def test_forward_windows_and_positions() -> None:
    index = build_kmer_index([rec("r1", "ACGAC")], 2, canonical=False)
    assert list(index) == ["AC", "CG", "GA"]
    assert index["AC"] == [
        KmerOccurrence(0, "r1", 0, "+"),
        KmerOccurrence(0, "r1", 3, "+"),
    ]
    assert index["CG"] == [KmerOccurrence(0, "r1", 1, "+")]
    assert index["GA"] == [KmerOccurrence(0, "r1", 2, "+")]


def test_canonical_picks_lexicographically_smaller() -> None:
    index = build_kmer_index([rec("r1", "ACG")], 3)
    # reverse complement of ACG is CGT; ACG < CGT, so the forward window wins
    assert list(index) == ["ACG"]
    assert index["ACG"] == [KmerOccurrence(0, "r1", 0, "+")]

    index = build_kmer_index([rec("r1", "CGT")], 3)
    # reverse complement of CGT is ACG; ACG < CGT, so the key is the
    # reverse complement and the strand is '-'
    assert list(index) == ["ACG"]
    assert index["ACG"] == [KmerOccurrence(0, "r1", 0, "-")]


def test_canonical_palindrome_stays_forward() -> None:
    index = build_kmer_index([rec("r1", "ACGT")], 4)
    assert index == {"ACGT": [KmerOccurrence(0, "r1", 0, "+")]}


def test_no_canonical_uses_forward_windows_only() -> None:
    index = build_kmer_index([rec("r1", "CGT")], 3, canonical=False)
    assert index == {"CGT": [KmerOccurrence(0, "r1", 0, "+")]}


def test_windows_with_non_acgt_symbols_are_skipped() -> None:
    index = build_kmer_index([rec("r1", "ACNGT")], 2, canonical=False)
    assert list(index) == ["AC", "GT"]


def test_records_shorter_than_k_contribute_nothing() -> None:
    index = build_kmer_index([rec("r1", "AC"), rec("r2", "ACGT")], 3, canonical=False)
    assert list(index) == ["ACG", "CGT"]
    assert all(occurrence.id == "r2" for occurrences in index.values() for occurrence in occurrences)


def test_duplicate_identifiers_are_not_merged() -> None:
    index = build_kmer_index([rec("x", "AAA"), rec("x", "AAA")], 2, canonical=False)
    assert index["AA"] == [
        KmerOccurrence(0, "x", 0, "+"),
        KmerOccurrence(0, "x", 1, "+"),
        KmerOccurrence(1, "x", 0, "+"),
        KmerOccurrence(1, "x", 1, "+"),
    ]


def test_keys_iterate_in_lexicographic_order() -> None:
    index = build_kmer_index([rec("r1", "TTTACGAAA")], 3, canonical=False)
    assert list(index) == sorted(index)


def test_batching_does_not_change_keys_or_positions() -> None:
    records = [rec("a", "ACGTAC"), rec("b", "TTT"), rec("c", "GG")]
    whole = build_kmer_index(iter(records), 2)
    piecewise = build_kmer_index(iter(records), 2)
    assert whole == piecewise
    assert list(whole) == list(piecewise)


def test_quality_values_are_ignored() -> None:
    with_quality = build_kmer_index(
        [rec("r1", "ACGT", quality=(0, 0, 0, 0))], 2, canonical=False
    )
    without = build_kmer_index([rec("r1", "ACGT")], 2, canonical=False)
    assert with_quality == without


def _fail_if_consumed():
    raise AssertionError("records were consumed")
    yield


@pytest.mark.parametrize("bad_k", [0, -1, True, False, 2.0, "3", None])
def test_invalid_k_raises_before_consuming_records(bad_k) -> None:
    with pytest.raises(ValueError):
        build_kmer_index(_fail_if_consumed(), bad_k)


@pytest.mark.parametrize("bad_canonical", [0, 1, "yes", None, 1.0])
def test_invalid_canonical_raises_before_consuming_records(bad_canonical) -> None:
    with pytest.raises(ValueError):
        build_kmer_index(_fail_if_consumed(), 3, canonical=bad_canonical)


def test_invalid_base_raises_when_record_is_reached() -> None:
    records = iter([rec("ok", "ACGT"), rec("bad", "ACXT")])
    with pytest.raises(SequenceValidationError) as excinfo:
        build_kmer_index(records, 2)
    message = str(excinfo.value)
    assert "bad" in message
    assert "X" in message
    assert "3" in message  # one-based position


def test_empty_input_gives_empty_index() -> None:
    assert build_kmer_index([], 3) == {}
