"""Tests for genome_variant.kmer."""

from __future__ import annotations

from collections import UserList
from itertools import chain

import pytest

from genome_variant.kmer import (
    KmerIndex,
    KmerOccurrence,
    build_kmer_index,
)
from genome_variant.sequence_io import SequenceRecord, SequenceValidationError


def records(*pairs) -> list[SequenceRecord]:
    return [SequenceRecord(identifier, sequence) for identifier, sequence in pairs]


def occurrence_map(index: KmerIndex) -> dict[str, list[KmerOccurrence]]:
    return {kmer: index.occurrences(kmer) for kmer in index.kmers()}


class TestBuildIndex:
    def test_canonical_basic(self) -> None:
        index = build_kmer_index(records(("r", "ACGT")), 2)
        # AC@0 stays +; CG@1 is a palindrome; GT@2 canonicalizes to AC(-).
        assert list(index) == ["AC", "CG"]
        assert occurrence_map(index) == {
            "AC": [KmerOccurrence(0, "r", 0, "+"), KmerOccurrence(0, "r", 2, "-")],
            "CG": [KmerOccurrence(0, "r", 1, "+")],
        }

    def test_non_canonical_only_forward(self) -> None:
        index = build_kmer_index(records(("r", "ACGT")), 2, canonical=False)
        assert list(index) == ["AC", "CG", "GT"]
        assert all(
            occurrence.strand == "+"
            for occurrences in occurrence_map(index).values()
            for occurrence in occurrences
        )

    def test_palindrome_is_forward(self) -> None:
        # ACGT is its own reverse complement; only one + occurrence.
        index = build_kmer_index(records(("p", "ACGT")), 4)
        assert index.occurrences("ACGT") == [KmerOccurrence(0, "p", 0, "+")]
        assert index.count("ACGT") == 1

    def test_reverse_choice_uses_minus(self) -> None:
        # Forward window "TT" has reverse complement "AA" (AA < TT), so
        # the canonical key is AA with a minus strand.
        index = build_kmer_index(records(("r", "TT")), 2)
        assert index.occurrences("AA") == [KmerOccurrence(0, "r", 0, "-")]
        assert "TT" not in index

    def test_k_one_keys_every_base(self) -> None:
        index = build_kmer_index(records(("r", "ACGT")), 1)
        # Canonical: A<->T -> A(+/-), C<->G -> C(+/-).
        assert index.occurrences("A") == [
            KmerOccurrence(0, "r", 0, "+"),
            KmerOccurrence(0, "r", 3, "-"),
        ]
        assert index.occurrences("C") == [
            KmerOccurrence(0, "r", 1, "+"),
            KmerOccurrence(0, "r", 2, "-"),
        ]

    def test_duplicate_identifiers_are_not_merged(self) -> None:
        index = build_kmer_index(
            records(("x", "AA"), ("x", "AA")), 2, canonical=False
        )
        assert index.occurrences("AA") == [
            KmerOccurrence(0, "x", 0, "+"),
            KmerOccurrence(1, "x", 0, "+"),
        ]

    def test_record_numbers_and_positions(self) -> None:
        index = build_kmer_index(
            records(("a", "AAAA"), ("b", "CCC"), ("c", "GGG")), 3, canonical=False
        )
        assert [o.record for o in index.occurrences("CCC")] == [1]
        assert [o.position for o in index.occurrences("AAA")] == [0, 1]

    def test_short_records_have_no_entries(self) -> None:
        index = build_kmer_index(records(("a", "AC"), ("b", ""), ("c", "G")), 3)
        assert len(index) == 0
        assert list(index.kmers()) == []

    def test_windows_with_ambiguity_symbols_are_skipped(self) -> None:
        index = build_kmer_index(records(("r", "ACGNTAC")), 2, canonical=False)
        assert list(index) == ["AC", "CG", "TA"]
        # The AC windows are at 0 and 5.
        assert [o.position for o in index.occurrences("AC")] == [0, 5]

    def test_iteration_is_lexicographic(self) -> None:
        index = build_kmer_index(records(("r", "TTTTGGGGAAAACCCC")), 2)
        keys = list(index.kmers())
        assert keys == sorted(keys)

    def test_batching_does_not_change_order(self) -> None:
        all_records = records(
            ("r0", "ACGT"), ("r1", "GGGG"), ("r2", "TTAA"), ("r3", "ACGNAC")
        )

        def as_pairs(index: KmerIndex):
            return [
                (kmer, [(o.record, o.id, o.position, o.strand) for o in occurrences])
                for kmer, occurrences in index.items()
            ]

        # The same logical record stream, chunked at the source in
        # different ways, must produce identical keys and positions.
        def chained(chunk_size: int):
            chunks = [
                all_records[i : i + chunk_size]
                for i in range(0, len(all_records), chunk_size)
            ]
            return chain.from_iterable(chunks)

        expected = as_pairs(build_kmer_index(iter(all_records), 3))
        assert as_pairs(build_kmer_index(chained(1), 3)) == expected
        assert as_pairs(build_kmer_index(chained(2), 3)) == expected
        assert as_pairs(build_kmer_index(chained(3), 3)) == expected
        assert as_pairs(build_kmer_index(chained(10), 3)) == expected

    def test_lowercase_is_invalid_in_handbuilt_records(self) -> None:
        # The reader upper-cases input; a directly constructed record with
        # lowercase symbols is treated like any other invalid symbol.
        with pytest.raises(SequenceValidationError, match="'b'"):
            build_kmer_index(records(("a", "AC"), ("b", "acgt")), 2)


class TestInvalidParameters:
    @pytest.mark.parametrize("value", [0, -1, -100, 1.0, "3", None])
    def test_bad_k_raises_value_error(self, value) -> None:
        with pytest.raises(ValueError):
            build_kmer_index(records(("r", "ACGT")), value)

    def test_boolean_k_rejected(self) -> None:
        with pytest.raises(ValueError, match="non-boolean"):
            build_kmer_index(records(("r", "ACGT")), True)

    @pytest.mark.parametrize("value", [0, 1, "true", None])
    def test_non_boolean_canonical_rejected(self, value) -> None:
        with pytest.raises(ValueError, match="canonical"):
            build_kmer_index(records(("r", "ACGT")), 2, canonical=value)

    def test_invalid_parameters_raised_before_consuming_records(self) -> None:
        consumed = []

        def source():
            for item in records(("r", "ACGT")):
                consumed.append(item)
                yield item

        with pytest.raises(ValueError):
            build_kmer_index(source(), 0)
        with pytest.raises(ValueError):
            build_kmer_index(source(), 2, canonical="yes")  # type: ignore[arg-type]
        assert consumed == []

    def test_iterator_source_is_consumed(self) -> None:
        def source():
            yield SequenceRecord("a", "ACGT")
            yield SequenceRecord("b", "AC")

        index = build_kmer_index(source(), 2, canonical=False)
        assert index.occurrences("AC") == [
            KmerOccurrence(0, "a", 0, "+"),
            KmerOccurrence(1, "b", 0, "+"),
        ]

    def test_userlist_source_accepted(self) -> None:
        index = build_kmer_index(UserList(records(("r", "ACGT"))), 2)
        assert index.count("AC") == 2


class TestLazyRecordErrors:
    def test_invalid_symbol_raises_when_record_reached(self) -> None:
        reached: list[str] = []

        def source():
            reached.append("good")
            yield SequenceRecord("good", "ACGT")
            reached.append("bad")
            yield SequenceRecord("bad", "ACQ")
            raise AssertionError("iteration must stop at the bad record")

        with pytest.raises(SequenceValidationError) as excinfo:
            build_kmer_index(source(), 2, canonical=False)
        # The first record was processed and the error points at the second.
        assert reached == ["good", "bad"]
        message = str(excinfo.value)
        assert "'bad'" in message
        assert "'Q'" in message
        assert "position 3" in message

    def test_error_identifies_symbol_and_one_based_position(self) -> None:
        with pytest.raises(SequenceValidationError) as excinfo:
            build_kmer_index(records(("z", "AC!GT")), 2)
        message = str(excinfo.value)
        assert "'z'" in message
        assert "'!'" in message
        assert "position 3" in message

    def test_first_record_error_propagates(self) -> None:
        def source():
            yield SequenceRecord("x", "AZ")

        with pytest.raises(SequenceValidationError, match="position 2"):
            build_kmer_index(source(), 2)

    def test_records_after_bad_record_are_not_consumed(self) -> None:
        def source():
            yield SequenceRecord("ok", "ACGT")
            yield SequenceRecord("bad", "Z")
            yield SequenceRecord("unreached", "ACGT")

        with pytest.raises(SequenceValidationError):
            build_kmer_index(source(), 2)
