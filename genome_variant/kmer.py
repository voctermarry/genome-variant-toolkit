"""k-mer index construction.

Public API:

- :class:`KmerOccurrence` — one occurrence of a k-mer: zero-based record
  index, record identifier, zero-based sequence offset and strand.
- :func:`build_kmer_index` — index an iterable of sequence records by
  k-mer, optionally canonicalizing each k-mer against its reverse
  complement.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from .sequence_io import SequenceRecord, _validate_sequence

__all__ = [
    "KmerOccurrence",
    "build_kmer_index",
]

_KMER_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class KmerOccurrence:
    """A single k-mer occurrence.

    ``record`` is the zero-based index of the record in the input stream
    and ``id`` its identifier (identifiers are never merged, so equal
    identifiers from different records stay distinct occurrences).
    ``position`` is the zero-based offset of the k-mer window in the
    record's sequence; ``strand`` is ``"+"`` when the indexed key is the
    forward window and ``"-"`` when it is the reverse complement.
    """

    record: int
    id: str
    position: int
    strand: str


def _reverse_complement(window: str) -> str:
    return window.translate(_COMPLEMENT)[::-1]


def build_kmer_index(
    records: Iterable[SequenceRecord],
    k: int,
    canonical: bool = True,
) -> dict[str, list[KmerOccurrence]]:
    """Index *records* by k-mer and return a key-sorted mapping.

    Every k-mer window of length *k* contributes one
    :class:`KmerOccurrence`, in input order, to the list stored under its
    key.  With *canonical* true (the default) the key is the
    lexicographically smaller of the forward window and its reverse
    complement; the strand is ``"+"`` for the forward window (always for
    palindromes) and ``"-"`` for the reverse complement.  With
    *canonical* false only forward windows are indexed, all on strand
    ``"+"``.

    Windows containing any symbol outside ``ACGT`` are skipped and
    records shorter than *k* contribute nothing.  The returned mapping
    iterates keys in lexicographic order; the same record stream always
    yields the same keys and occurrence order regardless of how the
    records are batched.

    *k* must be a non-boolean integer of at least 1 and *canonical* a
    boolean; invalid parameters raise :class:`ValueError` before
    *records* is consumed.  A record whose sequence holds a symbol
    outside the IUPAC set raises :class:`SequenceValidationError` when
    that record is reached.
    """
    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        kind = type(k).__name__
        raise ValueError(f"k must be a non-boolean integer of at least 1, not {kind}")
    if not isinstance(canonical, bool):
        kind = type(canonical).__name__
        raise ValueError(f"canonical must be a boolean, not {kind}")

    index: dict[str, list[KmerOccurrence]] = {}
    for record_index, record in enumerate(records):
        _validate_sequence(record.identifier, record.sequence)
        sequence = record.sequence
        for start in range(0, len(sequence) - k + 1):
            window = sequence[start : start + k]
            if not _KMER_BASES.issuperset(window):
                continue
            key = window
            strand = "+"
            if canonical:
                complement = _reverse_complement(window)
                if complement < window:
                    key = complement
                    strand = "-"
            occurrence = KmerOccurrence(record_index, record.identifier, start, strand)
            index.setdefault(key, []).append(occurrence)

    return dict(sorted(index.items()))
