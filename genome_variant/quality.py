"""Read quality trimming and filtering.

Public API:

- :class:`ReadQualityError` — a record's quality values are missing,
  misaligned with the sequence, or outside the Phred range 0–93.
- :func:`trim_and_filter_reads` — lazily trim low-quality bases from both
  ends of each read and drop reads that fail the length or mean-quality
  thresholds.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from .sequence_io import SequenceRecord

__all__ = [
    "ReadQualityError",
    "trim_and_filter_reads",
]

_MIN_QUALITY = 0
_MAX_QUALITY = 93


class ReadQualityError(ValueError):
    """A record's quality values are missing, misaligned, or out of range."""


def _check_quality_threshold(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer, got {value!r}")
    if not _MIN_QUALITY <= value <= _MAX_QUALITY:
        raise ValueError(
            f"{name} must be between {_MIN_QUALITY} and {_MAX_QUALITY}, got {value}"
        )


def _check_min_length(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"min_length must be an integer, got {value!r}")
    if value < 1:
        raise ValueError(f"min_length must be at least 1, got {value}")


def trim_and_filter_reads(
    records: Iterable[SequenceRecord],
    min_end_quality: int = 20,
    min_mean_quality: int = 20,
    min_length: int = 30,
) -> Iterator[SequenceRecord]:
    """Lazily trim and filter reads by base quality.

    For each record, bases whose Phred score is strictly below
    *min_end_quality* are removed from both ends (scores equal to the
    threshold are kept; low-quality bases inside the read never split it).
    A record is dropped when every base is trimmed away, when the trimmed
    length is below *min_length*, or when the sum of the remaining quality
    scores is less than ``min_mean_quality`` times the trimmed length.
    Surviving records keep their identifier and description with sequence
    and quality sliced in step; record order is preserved.

    The thresholds must be non-boolean integers; the quality thresholds
    are limited to 0–93 and *min_length* must be at least 1, otherwise
    :class:`ValueError` is raised immediately.  Records without quality
    values, with a quality length different from the sequence length, or
    with scores outside 0–93 raise :class:`ReadQualityError` when the
    iterator reaches them.
    """
    _check_quality_threshold("min_end_quality", min_end_quality)
    _check_quality_threshold("min_mean_quality", min_mean_quality)
    _check_min_length(min_length)
    return _trimmed_records(records, min_end_quality, min_mean_quality, min_length)


def _trimmed_records(
    records: Iterable[SequenceRecord],
    min_end_quality: int,
    min_mean_quality: int,
    min_length: int,
) -> Iterator[SequenceRecord]:
    for record in records:
        trimmed = _trim_ends(record, min_end_quality)
        if trimmed is None:
            continue
        sequence, quality = trimmed
        if len(sequence) < min_length:
            continue
        if sum(quality) < min_mean_quality * len(sequence):
            continue
        yield SequenceRecord(
            record.identifier, sequence, record.description, quality
        )


def _trim_ends(
    record: SequenceRecord, min_end_quality: int
) -> tuple[str, tuple[int, ...]] | None:
    identifier = record.identifier
    quality = record.quality
    if quality is None:
        raise ReadQualityError(f"record {identifier!r}: no quality values")
    if len(quality) != len(record.sequence):
        raise ReadQualityError(
            f"record {identifier!r}: quality length {len(quality)} does not "
            f"match sequence length {len(record.sequence)}"
        )
    for value in quality:
        if not _MIN_QUALITY <= value <= _MAX_QUALITY:
            raise ReadQualityError(
                f"record {identifier!r}: quality score {value} is outside the "
                f"range {_MIN_QUALITY}-{_MAX_QUALITY}"
            )

    sequence = record.sequence
    left = 0
    right = len(sequence)
    while left < right and quality[left] < min_end_quality:
        left += 1
    while right > left and quality[right - 1] < min_end_quality:
        right -= 1
    if left == right:
        return None
    return sequence[left:right], tuple(quality[left:right])
