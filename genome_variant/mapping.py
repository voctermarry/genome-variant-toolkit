"""Deterministic multi-reference mapping of sequencing reads.

Public API:

- :class:`MappingReferenceError` — the reference collection contained no
  records.
- :class:`ReadMapping` — the mapping result for one read: mapped or
  unmapped, with zero-based half-open coordinates and a CIGAR string.
- :func:`map_reads` — align every read against every reference record,
  on both strands, using the local-alignment scoring of
  :func:`~genome_variant.alignment.align_pair`.

Each read is aligned, as given and as its complete IUPAC reverse
complement, against every reference record.  Candidates scoring below
``min_score`` are discarded; the surviving candidate with the highest
score wins, with ties resolved by reference input index, reference
coordinates, strand, original-read query coordinates and CIGAR.  A read
with no qualifying candidate still produces an unmapped result.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from .alignment import align_pair, _score_parameter
from .sequence_io import SequenceRecord

__all__ = [
    "MappingReferenceError",
    "ReadMapping",
    "map_reads",
]

# Full IUPAC complementation; the table maps a base to its complement and
# the reversed translation is the reverse complement.
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class MappingReferenceError(ValueError):
    """The reference collection contained no reference records."""


@dataclass(frozen=True)
class ReadMapping:
    """The mapping of one read.

    For mapped reads ``reference_record`` is the zero-based index of the
    chosen reference in reference input order, ``reference`` its
    identifier, ``reference_start``/``reference_end`` the zero-based
    half-open reference interval, ``query_start``/``query_end`` the
    matching interval on the original read (zero-based, half-open),
    ``strand`` ``"+"`` or ``"-"``, ``score`` the local-alignment score
    and ``cigar`` the CIGAR relative to the aligned strand.  For unmapped
    reads all fields after ``mapped`` except ``record`` and ``id`` are
    ``None``.
    """

    record: int
    id: str
    mapped: bool
    reference_record: int | None
    reference: str | None
    reference_start: int | None
    reference_end: int | None
    query_start: int | None
    query_end: int | None
    strand: str | None
    score: int | None
    cigar: str | None

    def to_dict(self) -> dict[str, object]:
        """Return the fixed-order JSON-object representation."""
        return {
            "record": self.record,
            "id": self.id,
            "mapped": self.mapped,
            "reference_record": self.reference_record,
            "reference": self.reference,
            "reference_start": self.reference_start,
            "reference_end": self.reference_end,
            "query_start": self.query_start,
            "query_end": self.query_end,
            "strand": self.strand,
            "score": self.score,
            "cigar": self.cigar,
        }


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _unmapped(record_index: int, read: SequenceRecord) -> ReadMapping:
    return ReadMapping(
        record=record_index,
        id=read.identifier,
        mapped=False,
        reference_record=None,
        reference=None,
        reference_start=None,
        reference_end=None,
        query_start=None,
        query_end=None,
        strand=None,
        score=None,
        cigar=None,
    )


def map_reads(
    references: Iterable[SequenceRecord],
    reads: Iterable[SequenceRecord],
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
    min_score: int = 1,
) -> Iterator[ReadMapping]:
    """Map each read against every reference on both strands.

    Every read is locally aligned against every reference record twice:
    with its forward sequence and with its complete IUPAC reverse
    complement.  Only candidates whose score is at least *min_score*
    survive.  The winner is the highest-scoring candidate; ties are
    resolved by, in order, the reference input index,
    ``reference_start``, ``reference_end``, the forward strand (``"+"``)
    over the reverse strand (``"-"``), ``query_start`` and ``query_end``
    measured on the original read, and the lexicographically smallest
    CIGAR.  On the reverse strand the query interval is converted back to
    the original read while the CIGAR stays oriented against the
    reference and the reverse-complement read.

    Reads without a qualifying candidate yield an unmapped
    :class:`ReadMapping`; they are never dropped.  Results are produced in
    read input order and are byte-stable for fixed inputs regardless of
    batching.

    Scoring parameters have the same meaning as for
    :func:`~genome_variant.alignment.align_pair` in local mode;
    *min_score* is a non-boolean positive integer (default 1).  Invalid
    parameters raise :class:`ValueError` before either input is consumed;
    an empty reference collection raises :class:`MappingReferenceError`
    before *reads* is consumed.  Invalid sequence symbols raise
    :class:`~genome_variant.sequence_io.SequenceValidationError` when the
    offending reference or read record is reached.
    """
    match = _score_parameter("match_score", match_score, positive=True)
    mismatch = _score_parameter(
        "mismatch_penalty", mismatch_penalty, positive=False
    )
    open_penalty = _score_parameter("gap_open", gap_open, positive=False)
    extend_penalty = _score_parameter(
        "gap_extend", gap_extend, positive=False
    )
    if not isinstance(min_score, int) or isinstance(min_score, bool):
        kind = type(min_score).__name__
        raise ValueError(
            f"min_score must be a non-boolean integer, not {kind}"
        )
    if min_score < 1:
        raise ValueError("min_score must be a positive integer")

    reference_records = list(references)
    if not reference_records:
        raise MappingReferenceError("at least one reference record is required")

    def _iter() -> Iterator[ReadMapping]:
        for record_index, read in enumerate(reads):
            yield _map_one(
                record_index,
                read,
                reference_records,
                match,
                mismatch,
                open_penalty,
                extend_penalty,
                min_score,
            )

    return _iter()


def _map_one(
    record_index: int,
    read: SequenceRecord,
    references: list[SequenceRecord],
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
    min_score: int,
) -> ReadMapping:
    reverse_record = SequenceRecord(
        read.identifier,
        _reverse_complement(read.sequence),
        read.description,
    )
    read_length = len(read.sequence)

    # The tuple is ordered by the deterministic tie-break rules; the
    # smallest tuple wins among equal-score candidates.
    best_key: tuple[
        int, int, int, int, int, int, int, str
    ] | None = None
    best_alignment = None
    best_strand = ""

    for ref_index, reference in enumerate(references):
        for strand, query_record in (("+", read), ("-", reverse_record)):
            alignment = align_pair(
                reference,
                query_record,
                mode="local",
                match_score=match,
                mismatch_penalty=mismatch,
                gap_open=gap_open,
                gap_extend=gap_extend,
            )
            if alignment.score < min_score:
                continue

            if strand == "+":
                query_start = alignment.query_start
                query_end = alignment.query_end
                strand_rank = 0
            else:
                # Convert the interval on the reverse-complement read back
                # to zero-based half-open coordinates on the original read.
                query_start = read_length - alignment.query_end
                query_end = read_length - alignment.query_start
                strand_rank = 1

            key = (
                -alignment.score,
                ref_index,
                alignment.reference_start,
                alignment.reference_end,
                strand_rank,
                query_start,
                query_end,
                alignment.cigar,
            )
            if best_key is None or key < best_key:
                best_key = key
                best_alignment = alignment
                best_strand = strand

    if best_alignment is None or best_key is None:
        return _unmapped(record_index, read)

    ref_index = best_key[1]
    return ReadMapping(
        record=record_index,
        id=read.identifier,
        mapped=True,
        reference_record=ref_index,
        reference=references[ref_index].identifier,
        reference_start=best_alignment.reference_start,
        reference_end=best_alignment.reference_end,
        query_start=best_key[5],
        query_end=best_key[6],
        strand=best_strand,
        score=best_alignment.score,
        cigar=best_alignment.cigar,
    )
