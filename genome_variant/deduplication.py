"""Read-level deduplication.

Public API:

- :func:`deduplicate_reads` — group identical (or canonically equivalent)
  reads and return one representative record per group, ordered by each
  group's first occurrence in the input.
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

_IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_IUPAC_COMPLEMENT)[::-1]


def _quality_sum(record: SequenceRecord) -> int:
    identifier = record.identifier
    quality = record.quality
    assert quality is not None  # caller only handles quality-bearing records
    if len(quality) != len(record.sequence):
        raise ReadQualityError(
            f"record {identifier!r}: quality length {len(quality)} does not "
            f"match sequence length {len(record.sequence)}"
        )
    total = 0
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
        total += score
    return total


def deduplicate_reads(
    records: Iterable[SequenceRecord],
    canonical: bool = False,
) -> tuple[SequenceRecord, ...]:
    """Deduplicate *records* and return one representative per group.

    Records are grouped by their full upper-case sequence; with
    *canonical* true the grouping key is the lexicographically smaller of
    the sequence and its full IUPAC reverse complement, so the two
    orientations of the same read are duplicates of each other.  Records
    that merely share an identifier are never merged.

    Each group contributes exactly one representative, and the
    representatives are ordered by the input position of their group's
    first occurrence, so the result does not depend on how the record
    iterator is chunked:

    - Groups of records without quality values (FASTA) keep the earliest
      record.
    - Groups of records with quality values (FASTQ) keep the record with
      the highest Phred score sum; ties keep the earliest record.

    Representative records are returned unchanged: identifier,
    description, original orientation, sequence and quality values are
    preserved as-is and the sequence is never rewritten to the grouping
    key.  An empty *records* iterable yields an empty tuple.

    *canonical* must be a boolean; other values raise :class:`ValueError`
    before *records* is consumed.  A record whose sequence holds a symbol
    outside the IUPAC set raises :class:`SequenceValidationError`.  A
    record with quality values whose scores are missing, mismatched in
    number or outside the Phred+33 range 0-93 raises
    :class:`ReadQualityError`.  Mixing records with and without quality
    values in one input raises :class:`ValueError`.
    """
    if not isinstance(canonical, bool):
        kind = type(canonical).__name__
        raise ValueError(f"canonical must be a boolean, not {kind}")

    # key -> [representative, representative score (None without quality)]
    groups: dict[str, list[object]] = {}
    saw_quality = False
    saw_plain = False

    for record in records:
        _validate_sequence(record.identifier, record.sequence)
        quality = record.quality
        if quality is None:
            saw_plain = True
            score = None
        else:
            saw_quality = True
            score = _quality_sum(record)
        if saw_quality and saw_plain:
            raise ValueError(
                "input mixes records with and without quality values"
            )

        sequence = record.sequence
        key = sequence
        if canonical:
            complement = _reverse_complement(sequence)
            if complement < sequence:
                key = complement

        group = groups.get(key)
        if group is None:
            groups[key] = [record, score]
        elif score is not None and score > group[1]:
            # Strictly greater: quality ties keep the earliest record.
            group[0] = record
            group[1] = score

    # The dict preserves insertion order, which is exactly the order of
    # each group's first occurrence in the input.
    return tuple(group[0] for group in groups.values())
