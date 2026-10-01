"""Reusable k-mer indexing over sequence records.

Public API:

- :class:`KmerOccurrence` — one occurrence of a k-mer: zero-based record
  number, record identifier, zero-based sequence start and strand.
- :class:`KmerIndex` — mapping of k-mer keys to their occurrences,
  iterated in lexicographic key order.
- :func:`build_kmer_index` — build an index from a
  :class:`~genome_variant.sequence_io.SequenceRecord` iterable.

When canonical indexing is enabled a window is keyed by the
lexicographically smaller of the forward fragment and its reverse
complement; the strand is ``"+"`` for the forward choice and ``"-"`` for
the reverse-complement choice (palindromes are always ``"+"``).
Non-canonical indexing keys windows by their forward fragment and marks
every occurrence ``"+"``.  Windows containing symbols other than
``A``, ``C``, ``G`` and ``T`` are skipped.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from .sequence_io import SequenceRecord, SequenceValidationError

__all__ = [
    "KmerOccurrence",
    "KmerIndex",
    "build_kmer_index",
]

# The reader accepts the full IUPAC DNA alphabet; only ACGT windows are
# indexed, the other IUPAC symbols make windows ineligible.
_IUPAC_BASES = frozenset("ACGTRYSWKMBDHVN")
_KMER_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class KmerOccurrence:
    """One occurrence of a k-mer.

    ``record`` is the zero-based index of the record in the input stream,
    ``id`` its identifier, ``position`` the zero-based start of the window
    within the record sequence and ``strand`` ``"+"`` or ``"-"``.
    """

    record: int
    id: str
    position: int
    strand: str


class KmerIndex:
    """K-mer to occurrence-list mapping with sorted iteration.

    Occurrences within each key keep input order: records are processed in
    sequence and windows from left to right, so the same record stream
    produces the same keys and positions no matter how it is batched.
    """

    def __init__(self, k: int, canonical: bool) -> None:
        self._k = k
        self._canonical = canonical
        self._occurrences: dict[str, list[KmerOccurrence]] = {}

    @property
    def k(self) -> int:
        return self._k

    @property
    def canonical(self) -> bool:
        return self._canonical

    def kmers(self) -> Iterator[str]:
        """Yield the indexed k-mer keys in lexicographic order."""
        yield from sorted(self._occurrences)

    def occurrences(self, kmer: str) -> list[KmerOccurrence]:
        """Return the occurrences of *kmer* in input order (empty if absent)."""
        return list(self._occurrences.get(kmer, ()))

    def items(self) -> Iterator[tuple[str, list[KmerOccurrence]]]:
        """Yield ``(kmer, occurrences)`` pairs in lexicographic key order."""
        for kmer in sorted(self._occurrences):
            yield kmer, list(self._occurrences[kmer])

    def count(self, kmer: str) -> int:
        """Return the number of occurrences of *kmer*."""
        return len(self._occurrences.get(kmer, ()))

    def __len__(self) -> int:
        return len(self._occurrences)

    def __contains__(self, kmer: object) -> bool:
        return kmer in self._occurrences

    def __iter__(self) -> Iterator[str]:
        return self.kmers()

    def _add(self, kmer: str, occurrence: KmerOccurrence) -> None:
        self._occurrences.setdefault(kmer, []).append(occurrence)


def _reverse_complement(fragment: str) -> str:
    return fragment.translate(_COMPLEMENT)[::-1]


def _validate_record(record: SequenceRecord) -> None:
    """Check a hand-built record against the reader's IUPAC rule set.

    Records produced by the reader were already validated, but records
    constructed directly may carry anything; report an offending symbol
    the same way the reader does, with identifier and one-based position.
    IUPAC ambiguity symbols (N, R, ...) pass here and simply make every
    window containing them ineligible downstream.
    """
    for position, symbol in enumerate(record.sequence, start=1):
        if symbol not in _IUPAC_BASES:
            raise SequenceValidationError(
                f"record {record.identifier!r}: invalid symbol {symbol!r} at position {position}"
            )


def build_kmer_index(
    records: Iterable[SequenceRecord], k: int, canonical: bool = True
) -> KmerIndex:
    """Build a :class:`KmerIndex` from *records*.

    *k* must be a non-boolean integer of at least 1; *canonical* must be a
    boolean.  Invalid parameters raise :class:`ValueError` before *records*
    is consumed.  Records are processed in input order; duplicate
    identifiers are never merged and records shorter than *k* contribute
    no entries.  A record containing a symbol outside the IUPAC DNA
    alphabet raises :class:`SequenceValidationError` (identifying the
    record, symbol and one-based position) when that record is reached;
    IUPAC ambiguity symbols such as ``N`` are valid in a record but make
    every window containing them ineligible.
    """
    if not isinstance(k, int) or isinstance(k, bool):
        kind = type(k).__name__
        raise ValueError(f"k must be a non-boolean integer, not {kind}")
    if k < 1:
        raise ValueError("k must be at least 1")
    if not isinstance(canonical, bool):
        kind = type(canonical).__name__
        raise ValueError(f"canonical must be a boolean, not {kind}")

    index = KmerIndex(k, canonical)
    for record_number, record in enumerate(records):
        sequence = record.sequence
        _validate_record(record)
        last_start = len(sequence) - k
        for position in range(last_start + 1):
            fragment = sequence[position : position + k]
            if any(symbol not in _KMER_BASES for symbol in fragment):
                continue
            if canonical:
                reverse = _reverse_complement(fragment)
                if reverse < fragment:
                    key, strand = reverse, "-"
                else:
                    # Lexicographic ties are palindromes and stay forward.
                    key, strand = fragment, "+"
            else:
                key, strand = fragment, "+"
            index._add(
                key,
                KmerOccurrence(record_number, record.identifier, position, strand),
            )
    return index
