"""Deterministic read-level deduplication.

Public API:

- :func:`deduplicate_reads` — group identical records (optionally across
  strand via canonical IUPAC reverse complementation) and return one
  representative per group.
"""

from __future__ import annotations

from collections.abc import Iterable

from .quality import ReadQualityError
from .sequence_io import SequenceRecord, _validate_sequence

__all__ = [
    "deduplicate_reads",
]

_PHRED_MIN = 0
_PHRED_MAX = 93

# Full IUPAC complementation; the table maps a base to its complement and
# the reversed translation is the reverse complement.
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _check_quality(record: SequenceRecord) -> None:
    """Validate the quality values of a quality-bearing record."""
    identifier = record.identifier
    quality = record.quality
    if len(quality) != len(record.sequence):
        raise ReadQualityError(
            f"record {identifier!r}: quality length {len(quality)} does not "
            f"match sequence length {len(record.sequence)}"
        )
    for score in quality:
        if not isinstance(score, int):
            raise ReadQualityError(
                f"record {identifier!r}: quality score {score!r} is not an integer"
            )
        if not _PHRED_MIN <= score <= _PHRED_MAX:
            raise ReadQualityError(
                f"record {identifier!r}: quality score {score} is outside the "
                "Phred+33 range 0-93"
            )


def deduplicate_reads(
    records: Iterable[SequenceRecord],
    canonical: bool = False,
) -> tuple[SequenceRecord, ...]:
    """Deduplicate *records* deterministically and return representatives.

    By default records whose complete upper-case sequences are identical
    form one group.  With *canonical* true the group key is the
    lexicographically smaller of the sequence and its full IUPAC reverse
    complement, so a sequence and its reverse complement are treated as
    duplicates.  Each group yields exactly one record:

    - records without quality values (FASTA): the earliest record seen;
    - records with quality values (FASTQ): the record with the greatest
      Phred quality sum, falling back to the earliest on a tie.

    Representative records are returned unchanged (identifier,
    description, original orientation, sequence and quality are never
    rewritten), in the order in which each group first appears in the
    input, so the result is independent of how the iterable is batched.
    Equal identifiers never merge records with different sequences.

    *canonical* must be a boolean; any other type raises
    :class:`ValueError` before *records* is consumed.  An empty iterable
    returns an empty tuple.  Invalid IUPAC symbols raise
    :class:`SequenceValidationError`; missing, mismatched or out-of-range
    quality values raise :class:`ReadQualityError`; mixing records with
    and without quality values in the same iterable raises
    :class:`ValueError`.
    """
    if not isinstance(canonical, bool):
        kind = type(canonical).__name__
        raise ValueError(f"canonical must be a boolean, not {kind}")

    representatives: dict[str, SequenceRecord] = {}
    quality_sums: dict[str, int] = {}
    quality_records: bool | None = None

    for record in records:
        _validate_sequence(record.identifier, record.sequence)

        has_quality = record.quality is not None
        if quality_records is None:
            quality_records = has_quality
        elif has_quality != quality_records:
            if has_quality:
                raise ValueError(
                    f"record {record.identifier!r}: records with quality values "
                    "cannot be mixed with records without quality values"
                )
            raise ValueError(
                f"record {record.identifier!r}: records without quality values "
                "cannot be mixed with records with quality values"
            )

        if has_quality:
            _check_quality(record)

        key = record.sequence
        if canonical:
            reverse_complement = _reverse_complement(key)
            if reverse_complement < key:
                key = reverse_complement

        if key not in representatives:
            representatives[key] = record
            if has_quality:
                quality_sums[key] = sum(record.quality)
        elif has_quality:
            score = sum(record.quality)
            # Strictly greater: a tie keeps the earlier record.
            if score > quality_sums[key]:
                representatives[key] = record
                quality_sums[key] = score

    return tuple(representatives.values())
