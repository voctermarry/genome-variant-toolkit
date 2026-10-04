"""Reference coverage and concordance statistics from mapped reads.

Public API:

- :class:`ReferenceCoverage` — the per-reference-record statistics:
  mapped read count, depth of coverage, breadth of coverage and the
  share of observed bases that literally match the reference.
- :func:`coverage_report` — map every read with the same
  local-alignment scoring, strand search and candidate adjudication as
  :func:`~genome_variant.mapping.map_reads` and aggregate one
  :class:`ReferenceCoverage` per reference record, in reference input
  order.

Only the single winning mapping of each read contributes; unmapped
reads contribute nothing.  CIGAR ``=`` and ``X`` operations each add
one to the depth of the corresponding reference position, ``D`` adds
nothing and ``I`` has no reference position.  A base is concordant
only when the observed character is literally identical to the
reference character, including IUPAC ambiguity symbols.  FASTQ quality
values never participate.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from ._columns import alignment_columns
from .alignment import _score_parameter
from .mapping import MappingReferenceError, _best_candidate
from .sequence_io import SequenceRecord

__all__ = [
    "ReferenceCoverage",
    "coverage_report",
]


@dataclass(frozen=True)
class ReferenceCoverage:
    """Coverage and concordance statistics for one reference record.

    ``reference_record`` is the zero-based index of the reference in
    reference input order and ``reference`` its identifier.
    ``mapped_reads`` counts reads whose unique winning mapping landed
    on this reference.  ``observed_bases`` is the total depth (the
    number of ``=``/``X`` columns over all those reads),
    ``covered_bases`` the number of positions with non-zero depth and
    ``concordant_bases`` the number of observed bases literally equal
    to the reference character.  The three fractions are fixed
    six-decimal strings: ``coverage_fraction`` is
    ``covered_bases/length``, ``mean_depth`` is
    ``observed_bases/length`` and ``concordance_fraction`` is
    ``concordant_bases/observed_bases``; a zero denominator yields
    ``"0.000000"``.
    """

    reference_record: int
    reference: str
    length: int
    mapped_reads: int
    observed_bases: int
    covered_bases: int
    coverage_fraction: str
    mean_depth: str
    concordant_bases: int
    concordance_fraction: str

    def to_dict(self) -> dict[str, object]:
        """Return the fixed-order JSON-object representation."""
        return {
            "reference_record": self.reference_record,
            "reference": self.reference,
            "length": self.length,
            "mapped_reads": self.mapped_reads,
            "observed_bases": self.observed_bases,
            "covered_bases": self.covered_bases,
            "coverage_fraction": self.coverage_fraction,
            "mean_depth": self.mean_depth,
            "concordant_bases": self.concordant_bases,
            "concordance_fraction": self.concordance_fraction,
        }


def _six_decimals(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "0.000000"
    return f"{numerator / denominator:.6f}"


def coverage_report(
    references: Iterable[SequenceRecord],
    reads: Iterable[SequenceRecord],
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
    min_score: int = 1,
) -> Iterator[ReferenceCoverage]:
    """Aggregate per-reference coverage and concordance statistics.

    Each read is mapped exactly as in
    :func:`~genome_variant.mapping.map_reads` (same local-alignment
    scoring, both strands, same candidate adjudication) and only the
    unique winning mapping is used; reads without a qualifying
    candidate contribute nothing.  Every ``=`` or ``X`` CIGAR column of
    a winning mapping adds one to the depth of its reference position;
    ``D`` columns add nothing and ``I`` columns have no reference
    position.  An observed base is concordant only when its character
    is literally identical to the reference character, including IUPAC
    ambiguity symbols.  Quality values never participate.

    One :class:`ReferenceCoverage` is produced per reference record, in
    reference input order, including records without any mapped read;
    duplicate identifiers are not merged and are distinguished by
    ``reference_record``.  Results are byte-stable for fixed inputs
    regardless of batching.

    Scoring parameters have the same meaning as for
    :func:`~genome_variant.mapping.map_reads`.  Invalid parameters
    raise :class:`ValueError` before either input is consumed; an
    empty reference collection raises
    :class:`~genome_variant.mapping.MappingReferenceError` before
    *reads* is consumed.  Invalid sequence symbols raise
    :class:`~genome_variant.sequence_io.SequenceValidationError` when
    the offending reference or read record is reached.
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

    def _iter() -> Iterator[ReferenceCoverage]:
        depths = [[0] * len(record.sequence) for record in reference_records]
        mapped_reads = [0] * len(reference_records)
        concordant = [0] * len(reference_records)

        for read in reads:
            winner = _best_candidate(
                read,
                reference_records,
                match,
                mismatch,
                open_penalty,
                extend_penalty,
                min_score,
            )
            if winner is None:
                continue
            ref_index, alignment, _strand = winner
            mapped_reads[ref_index] += 1
            concordant[ref_index] += _collect(
                alignment_columns(alignment), depths[ref_index]
            )

        for ref_index, reference in enumerate(reference_records):
            depth = depths[ref_index]
            length = len(reference.sequence)
            observed = sum(depth)
            covered = sum(1 for position_depth in depth if position_depth)
            yield ReferenceCoverage(
                reference_record=ref_index,
                reference=reference.identifier,
                length=length,
                mapped_reads=mapped_reads[ref_index],
                observed_bases=observed,
                covered_bases=covered,
                coverage_fraction=_six_decimals(covered, length),
                mean_depth=_six_decimals(observed, length),
                concordant_bases=concordant[ref_index],
                concordance_fraction=_six_decimals(
                    concordant[ref_index], observed
                ),
            )

    return _iter()


def _collect(columns, depth: list[int]) -> int:
    """Add one read's winning alignment to *depth*; return concordance.

    Columns share the canonical interpretation used by every other
    post-mapping entry: ``=``/``X`` add one observation at their
    reference position, ``D`` consumes the reference without depth and
    ``I`` has no reference position.  The observed base is already
    expressed on the reference's forward strand on a reverse-strand win.
    """
    concordant = 0
    for column in columns:
        if column.op in ("=", "X"):
            position = column.reference_position
            assert position is not None
            depth[position] += 1
            if column.query_base == column.reference_base:
                concordant += 1
    return concordant
