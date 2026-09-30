"""Read quality trimming and filtering.

Public API:

- :class:`ReadQualityError` — a record lacks quality values, carries the
  wrong number of them, or holds scores outside Phred+33 range 0-93.
- :func:`filter_reads` — lazily trim low-quality ends and drop records by
  mean quality and minimum length.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from .sequence_io import SequenceRecord

__all__ = [
    "ReadQualityError",
    "filter_reads",
]

_PHRED_MIN = 0
_PHRED_MAX = 93


class ReadQualityError(ValueError):
    """A record has no quality values, a quality/sequence length mismatch,
    or a quality score outside the Phred+33 range 0-93.

    The message always contains the offending record's identifier.
    """


def _threshold(name: str, value: object, *, allow_zero: bool) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if allow_zero:
        if not _PHRED_MIN <= value <= _PHRED_MAX:
            raise ValueError(f"{name} must be between {_PHRED_MIN} and {_PHRED_MAX}")
    elif value < 1:
        raise ValueError(f"{name} must be at least 1")
    return value


def filter_reads(
    records: Iterable[SequenceRecord],
    min_end_quality: int = 20,
    min_mean_quality: int = 20,
    min_length: int = 30,
) -> Iterator[SequenceRecord]:
    """Yield records after end trimming and quality/length filtering.

    Each record is processed lazily, in iteration order:

    1. Bases whose Phred score is strictly below *min_end_quality* are
       removed consecutively from both ends; bases equal to the threshold
       are kept and low-quality bases inside the read never split it.
    2. A record is dropped when every base was trimmed, when the trimmed
       length is below *min_length*, or when the mean quality (sum of
       scores divided by trimmed length) is below *min_mean_quality*.
    3. Surviving records keep their identifier and description; sequence
       and quality are sliced together.

    *min_end_quality* and *min_mean_quality* are non-boolean integers in
    the Phred range 0-93 (defaults 20 and 20); *min_length* is a
    non-boolean integer of at least 1 (default 30).  Invalid parameters
    raise :class:`ValueError`; records without quality values, with a
    quality/sequence length mismatch or with scores outside 0-93 raise
    :class:`ReadQualityError` (mentioning the record identifier) when
    that record is reached during iteration.
    """
    end_threshold = _threshold(
        "min_end_quality", min_end_quality, allow_zero=True
    )
    mean_threshold = _threshold(
        "min_mean_quality", min_mean_quality, allow_zero=True
    )
    length_threshold = _threshold("min_length", min_length, allow_zero=False)

    def _iter() -> Iterator[SequenceRecord]:
        # Record problems surface lazily, when each record is reached.
        for record in records:
            trimmed = _filter_one(
                record, end_threshold, mean_threshold, length_threshold
            )
            if trimmed is not None:
                yield trimmed

    return _iter()


def _filter_one(
    record: SequenceRecord,
    end_threshold: int,
    mean_threshold: int,
    length_threshold: int,
) -> SequenceRecord | None:
    identifier = record.identifier
    quality = record.quality
    sequence = record.sequence

    if quality is None:
        raise ReadQualityError(
            f"record {identifier!r}: has no quality values"
        )
    if len(quality) != len(sequence):
        raise ReadQualityError(
            f"record {identifier!r}: quality length {len(quality)} does not "
            f"match sequence length {len(sequence)}"
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

    left = 0
    right = len(quality)
    while left < right and quality[left] < end_threshold:
        left += 1
    while right > left and quality[right - 1] < end_threshold:
        right -= 1

    trimmed_length = right - left
    if trimmed_length < length_threshold:
        return None
    if sum(quality[left:right]) < mean_threshold * trimmed_length:
        return None

    return SequenceRecord(
        record.identifier,
        sequence[left:right],
        record.description,
        quality[left:right],
    )
