"""Map reads against a multi-record reference.

Public API:

- :class:`ReadMapping` — the mapping outcome of one read: zero-based read
  record index, read identifier, mapped flag and, for mapped reads, the
  chosen reference record, zero-based half-open coordinates, strand,
  score and CIGAR string.
- :class:`MappingReferenceError` — raised when the reference set is empty.
- :func:`map_reads` — lazily map every read against every reference
  record on both strands and yield one :class:`ReadMapping` per read in
  input order.

Each read is aligned locally (the ``align_pair`` local mode and its
scoring parameters) against every reference record twice: once as the
forward sequence and once as the full IUPAC reverse complement.  Only
candidates scoring at least ``min_score`` are kept.  FASTQ quality
values never participate in scoring.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from .alignment import align_pair, _score_parameter
from .sequence_io import SequenceRecord, _validate_sequence

__all__ = [
    "MappingReferenceError",
    "ReadMapping",
    "map_reads",
]

# Full IUPAC DNA complement: A<->T, C<->G, R<->Y, S<->S, W<->W, K<->M,
# B<->V, D<->H, N<->N.
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


class MappingReferenceError(ValueError):
    """The reference set passed to :func:`map_reads` is empty."""


@dataclass(frozen=True)
class ReadMapping:
    """The mapping outcome of a single read.

    ``record`` is the zero-based index of the read in the input stream
    and ``id`` its identifier (identifiers are never merged, so equal
    identifiers from different reads stay distinct outcomes).  For a
    mapped read ``reference_record`` is the zero-based index of the
    chosen reference record, ``reference`` its identifier, the four
    coordinates are zero-based and half-open on the reference and on the
    original read respectively, ``strand`` is ``"+"`` or ``"-"``, and
    ``score`` and ``cigar`` describe the chosen local alignment.  For an
    unmapped read every field after ``mapped`` is ``None``.
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


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def map_reads(
    references: Iterable[SequenceRecord],
    reads: Iterable[SequenceRecord],
    min_score: int = 1,
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
) -> Iterator[ReadMapping]:
    """Map *reads* against *references* and yield one outcome per read.

    Every read is aligned locally against every reference record, both as
    the forward sequence and as its full IUPAC reverse complement, using
    the ``align_pair`` local mode and scoring parameters.  Alignments
    scoring below *min_score* are discarded.  Among the remaining
    candidates the highest score wins; ties are broken by the smallest
    reference input index, ``reference_start`` and ``reference_end``,
    then the forward strand, then the smallest ``query_start`` and
    ``query_end`` converted back onto the original read, and finally the
    lexicographically smallest CIGAR string.  A read without any
    qualifying candidate still yields an unmapped outcome.

    Outcomes are produced in read input order; the same records always
    yield the same outcomes regardless of how the iterables are batched.
    For reverse-strand outcomes the query coordinates are converted to
    the zero-based half-open interval on the original read, while the
    CIGAR keeps its orientation against the reference and the reverse
    complement read.

    *min_score* must be a non-boolean positive integer; the scoring
    parameters follow the same rules as :func:`align_pair`.  Invalid
    parameters raise :class:`ValueError` before *reads* is consumed.  An
    empty reference set raises :class:`MappingReferenceError`; a record
    whose sequence holds a symbol outside the IUPAC set raises
    :class:`~genome_variant.sequence_io.SequenceValidationError`.
    """
    if (
        not isinstance(min_score, int)
        or isinstance(min_score, bool)
        or min_score < 1
    ):
        kind = type(min_score).__name__
        raise ValueError(
            f"min_score must be a non-boolean positive integer, not {kind}"
        )
    match = _score_parameter("match_score", match_score, positive=True)
    mismatch = _score_parameter(
        "mismatch_penalty", mismatch_penalty, positive=False
    )
    open_penalty = _score_parameter("gap_open", gap_open, positive=False)
    extend_penalty = _score_parameter("gap_extend", gap_extend, positive=False)
    return _map_reads(
        references,
        reads,
        min_score,
        match,
        mismatch,
        open_penalty,
        extend_penalty,
    )


def _map_reads(
    references: Iterable[SequenceRecord],
    reads: Iterable[SequenceRecord],
    min_score: int,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
) -> Iterator[ReadMapping]:
    reference_list = list(references)
    if not reference_list:
        raise MappingReferenceError("the reference set is empty")
    for reference in reference_list:
        _validate_sequence(reference.identifier, reference.sequence)
    for record_index, read in enumerate(reads):
        yield _map_read(
            record_index,
            read,
            reference_list,
            min_score,
            match,
            mismatch,
            gap_open,
            gap_extend,
        )


def _map_read(
    record_index: int,
    read: SequenceRecord,
    references: list[SequenceRecord],
    min_score: int,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
) -> ReadMapping:
    _validate_sequence(read.identifier, read.sequence)
    length = len(read.sequence)
    reverse = _reverse_complement(read.sequence)

    best_key: tuple[int, int, int, int, int, int, int, str] | None = None
    best: ReadMapping | None = None
    for reference_index, reference in enumerate(references):
        for strand, query_sequence in (
            ("+", read.sequence),
            ("-", reverse),
        ):
            alignment = align_pair(
                reference,
                SequenceRecord(read.identifier, query_sequence),
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
            else:
                # Convert the reverse-complement interval back onto the
                # original read; the CIGAR keeps its own orientation.
                query_start = length - alignment.query_end
                query_end = length - alignment.query_start
            key = (
                -alignment.score,
                reference_index,
                alignment.reference_start,
                alignment.reference_end,
                0 if strand == "+" else 1,
                query_start,
                query_end,
                alignment.cigar,
            )
            if best_key is None or key < best_key:
                best_key = key
                best = ReadMapping(
                    record=record_index,
                    id=read.identifier,
                    mapped=True,
                    reference_record=reference_index,
                    reference=reference.identifier,
                    reference_start=alignment.reference_start,
                    reference_end=alignment.reference_end,
                    query_start=query_start,
                    query_end=query_end,
                    strand=strand,
                    score=alignment.score,
                    cigar=alignment.cigar,
                )

    if best is None:
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
    return best
