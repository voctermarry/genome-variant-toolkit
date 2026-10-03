"""Tests for genome_variant.deduplication."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from genome_variant.deduplication import deduplicate_reads
from genome_variant.quality import ReadQualityError
from genome_variant.sequence_io import (
    SequenceRecord,
    SequenceValidationError,
)


def fa(identifier: str, sequence: str, description: str = "") -> SequenceRecord:
    return SequenceRecord(identifier, sequence, description)


def fq(
    identifier: str,
    sequence: str,
    quality,
    description: str = "",
) -> SequenceRecord:
    return SequenceRecord(identifier, sequence, description, tuple(quality))


class TestBasicGrouping:
    def test_empty_iterable_returns_empty_tuple(self) -> None:
        result = deduplicate_reads(iter(()))
        assert result == ()
        assert isinstance(result, tuple)

    def test_empty_iterable_canonical_returns_empty_tuple(self) -> None:
        assert deduplicate_reads(iter(()), canonical=True) == ()

    def test_identical_fasta_sequences_are_grouped(self) -> None:
        records = [fa("a", "ACGT"), fa("b", "ACGT", "dup"), fa("c", "TTTT")]
        result = deduplicate_reads(iter(records))
        assert result == (fa("a", "ACGT"), fa("c", "TTTT"))

    def test_fasta_keeps_earliest_record_with_its_metadata(self) -> None:
        records = [fa("first", "ACGT", "kept desc"), fa("second", "ACGT")]
        (representative,) = deduplicate_reads(iter(records))
        assert representative is records[0]
        assert representative.identifier == "first"
        assert representative.description == "kept desc"

    def test_equal_identifiers_with_different_sequences_are_not_merged(self) -> None:
        records = [fa("x", "ACGT"), fa("x", "ACGA")]
        result = deduplicate_reads(iter(records))
        assert len(result) == 2

    def test_reverse_complements_are_distinct_without_canonical(self) -> None:
        records = [fa("a", "AAAA"), fa("b", "TTTT")]
        result = deduplicate_reads(iter(records))
        assert len(result) == 2

    def test_result_order_tracks_first_appearance(self) -> None:
        records = [
            fa("a", "GG"),
            fa("b", "AA"),
            fa("c", "GG"),
            fa("d", "CC"),
            fa("e", "AA"),
        ]
        result = deduplicate_reads(iter(records))
        assert [r.sequence for r in result] == ["GG", "AA", "CC"]

    def test_result_is_independent_of_batching(self) -> None:
        records = [
            fa("a", "GG"), fa("b", "AA"), fa("c", "GG"),
            fa("d", "CC"), fa("e", "AA"),
        ]

        def chunks(data, size) -> Iterator[SequenceRecord]:
            for start in range(0, len(data), size):
                yield from data[start : start + size]

        once = deduplicate_reads(iter(records))
        for size in (1, 2, 3, 5):
            assert deduplicate_reads(chunks(records, size)) == once


class TestFastqRepresentative:
    def test_highest_quality_sum_wins(self) -> None:
        records = [
            fq("low", "ACGT", (1, 1, 1, 1)),
            fq("high", "ACGT", (40, 40, 40, 40)),
            fq("mid", "ACGT", (10, 10, 10, 10)),
        ]
        (representative,) = deduplicate_reads(iter(records))
        assert representative is records[1]

    def test_tie_keeps_earliest(self) -> None:
        records = [
            fq("first", "ACGT", (10, 20, 30, 40)),
            fq("second", "ACGT", (40, 30, 20, 10)),
        ]
        (representative,) = deduplicate_reads(iter(records))
        assert representative is records[0]

    def test_winner_keeps_own_quality_and_metadata(self) -> None:
        records = [
            fq("low", "ACGT", (0, 0, 0, 0), "low desc"),
            fq("high", "ACGT", (40, 40, 40, 40), "high desc"),
        ]
        (representative,) = deduplicate_reads(iter(records))
        assert representative.identifier == "high"
        assert representative.description == "high desc"
        assert representative.quality == (40, 40, 40, 40)
        assert representative.sequence == "ACGT"

    def test_distinct_sequences_are_kept(self) -> None:
        records = [
            fq("a", "ACGT", (40, 40, 40, 40)),
            fq("b", "TTTT", (1, 1, 1, 1)),
        ]
        result = deduplicate_reads(iter(records))
        assert [r.identifier for r in result] == ["a", "b"]


class TestCanonical:
    def test_sequence_and_reverse_complement_are_grouped(self) -> None:
        records = [fa("a", "AAAA"), fa("b", "TTTT")]
        (representative,) = deduplicate_reads(iter(records), canonical=True)
        assert representative is records[0]

    def test_reverse_complement_seen_first_is_kept_unchanged(self) -> None:
        records = [fq("rc", "TTTT", (5, 6, 7, 8)), fq("fwd", "AAAA", (1, 2, 3, 4))]
        (representative,) = deduplicate_reads(iter(records), canonical=True)
        # Highest quality wins regardless of orientation; the record is
        # emitted exactly as stored (sequence/quality not rewritten).
        assert representative.identifier == "rc"
        assert representative.sequence == "TTTT"
        assert representative.quality == (5, 6, 7, 8)

    def test_representative_sequence_is_not_rewritten_to_group_key(self) -> None:
        # ACGT's reverse complement is also ACGT? No: rc(ACGT) = ACGT.
        # Use a sequence whose reverse complement is lexicographically
        # smaller but differs.
        forward = "CCAA"            # rc -> TTGG (TTGG > CCAA, key CCAA)
        other = "TTGG"
        records = [fa("fwd", forward), fa("rev", other)]
        (representative,) = deduplicate_reads(iter(records), canonical=True)
        assert representative.sequence == forward

        forward2 = "TTGG"           # rc -> CCAA (CCAA < TTGG, key CCAA)
        records = [fa("rev", forward2), fa("fwd", "CCAA")]
        result = deduplicate_reads(iter(records), canonical=True)
        assert len(result) == 1
        assert result[0].sequence == "TTGG"  # original orientation kept

    def test_full_iupac_reverse_complement(self) -> None:
        # R (A/G) complements to Y (C/T); rc("RY") = "RY" (palindrome).
        records = [fa("a", "RY"), fa("b", "RY")]
        assert len(deduplicate_reads(iter(records), canonical=True)) == 1
        # rc("R") = "Y"; R < Y so the canonical key is R.
        records = [fa("a", "R"), fa("b", "Y")]
        (representative,) = deduplicate_reads(iter(records), canonical=True)
        assert representative.sequence == "R"

    def test_canonical_groups_across_multiple_groups(self) -> None:
        records = [
            fa("a", "AAAA"),
            fa("b", "CCCC"),
            fa("c", "TTTT"),  # rc(AAAA) = TTTT
            fa("d", "GGGG"),  # rc(CCCC) = GGGG
        ]
        result = deduplicate_reads(iter(records), canonical=True)
        assert [r.sequence for r in result] == ["AAAA", "CCCC"]

    def test_canonical_quality_winner_across_orientation(self) -> None:
        records = [
            fq("fwd", "AAAA", (1, 1, 1, 1)),
            fq("rev", "TTTT", (9, 9, 9, 9)),
        ]
        (representative,) = deduplicate_reads(iter(records), canonical=True)
        assert representative.identifier == "rev"
        assert representative.sequence == "TTTT"
        assert representative.quality == (9, 9, 9, 9)


class TestValidation:
    def test_non_boolean_canonical_raises_before_consuming_input(self) -> None:
        for value in (1, 0, "true", None, 1.0):
            consumed = []

            def records() -> Iterator[SequenceRecord]:
                consumed.append(True)
                yield fa("a", "ACGT")

            with pytest.raises(ValueError):
                deduplicate_reads(records(), canonical=value)
            assert consumed == []

    def test_boolean_subclass_value_is_accepted(self) -> None:
        assert deduplicate_reads(iter([fa("a", "ACGT")]), canonical=False) == (
            fa("a", "ACGT"),
        )

    def test_invalid_iupac_base_raises(self) -> None:
        with pytest.raises(SequenceValidationError):
            deduplicate_reads(iter([fa("a", "ACGZ")]))

    def test_invalid_base_in_later_record_raises(self) -> None:
        with pytest.raises(SequenceValidationError):
            deduplicate_reads(
                iter([fa("a", "ACGT"), fa("b", "TTTZ")]), canonical=True
            )

    def test_quality_length_mismatch_raises(self) -> None:
        with pytest.raises(ReadQualityError):
            deduplicate_reads(iter([fq("a", "ACGT", (40, 40, 40))]))

    def test_quality_out_of_range_raises(self) -> None:
        with pytest.raises(ReadQualityError):
            deduplicate_reads(iter([fq("a", "ACGT", (40, 40, 94, 40))]))
        with pytest.raises(ReadQualityError):
            deduplicate_reads(iter([fq("a", "ACGT", (40, -1, 40, 40))]))

    def test_non_integer_quality_raises(self) -> None:
        with pytest.raises(ReadQualityError):
            deduplicate_reads(iter([fq("a", "ACGT", (40, 40, "x", 40))]))

    def test_mixed_quality_and_quality_free_records_raise(self) -> None:
        with pytest.raises(ValueError):
            deduplicate_reads(
                iter([fa("a", "ACGT"), fq("b", "ACGT", (1, 2, 3, 4))])
            )
        with pytest.raises(ValueError):
            deduplicate_reads(
                iter([fq("a", "ACGT", (1, 2, 3, 4)), fa("b", "ACGT")])
            )

    def test_mixed_error_is_not_quality_or_validation_error(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            deduplicate_reads(
                iter([fa("a", "ACGT"), fq("b", "ACGT", (1, 2, 3, 4))])
            )
        assert not isinstance(exc_info.value, ReadQualityError)
        assert not isinstance(exc_info.value, SequenceValidationError)
