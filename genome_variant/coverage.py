"""Reference coverage and concordance statistics from mapped reads.

Public API:

- :class:`CoverageStats` — per-reference coverage totals and the fixed
  JSON-object representation rendered by the ``coverage-report``
  subcommand.
- :func:`coverage_report` — map every read with the same local-alignment
  scoring, both-strand search and candidate adjudication as
  :func:`~genome_variant.mapping.map_reads`, then tally per-reference
  depth and verbatim concordance.

Only a read's unique winning mapping contributes; unmapped reads are
ignored.  In the winning CIGAR each ``=`` and ``X`` column adds one depth
observation at its reference position, ``D`` advances the reference
position without adding depth and ``I`` has no reference position.
Concordance is strict character equality (including IUPAC characters),
which is exactly what the ``=`` operation records.  FASTQ quality values
never participate.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from .alignment import _expand_cigar, _score_parameter
from .mapping import MappingReferenceError, _best_candidate
from .sequence_io import SequenceRecord

__all__ = [
    "CoverageStats",
    "MappingReferenceError",
    "coverage_report",
]


@dataclass(frozen=True)
class CoverageStats:
    """Coverage and concordance totals for one reference record.

    ``reference_record`` is the zero-based index in reference input
    order, ``reference`` the identifier and ``length`` the reference
    sequence length.  ``mapped_reads`` counts reads whose unique winning
    mapping lands on this record.  ``observed_bases`` is the total number
    of reference-consuming (``=``/``X``) alignment columns,
    ``covered_bases`` the number of reference positions with non-zero
    depth and ``concordant_bases`` the number of ``=`` columns.  The
    three ratios are rendered as fixed six-decimal strings by
    :meth:`to_dict`.
    """

    reference_record: int
    reference: str
    length: int
    mapped_reads: int
    observed_bases: int
    covered_bases: int
    concordant_bases: int

    def to_dict(self) -> dict[str, object]:
        """Return the fixed-order JSON-object representation.

        ``coverage_fraction`` is ``covered_bases/length``,
        ``mean_depth`` is ``observed_bases/length`` and
        ``concordance_fraction`` is ``concordant_bases/observed_bases``;
        a zero denominator yields ``"0.000000"``.
        """
        coverage = (
            self.covered_bases / self.length if self.length else 0.0
        )
        mean_depth = (
            self.observed_bases / self.length if self.length else 0.0
        )
        concordance = (
            self.concordant_bases / self.observed_bases
            if self.observed_bases
            else 0.0
        )
        return {
            "reference_record": self.reference_record,
            "reference": self.reference,
            "length": self.length,
            "mapped_reads": self.mapped_reads,
            "observed_bases": self.observed_bases,
            "covered_bases": self.covered_bases,
            "coverage_fraction": f"{coverage:.6f}",
            "mean_depth": f"{mean_depth:.6f}",
            "concordant_bases": self.concordant_bases,
            "concordance_fraction": f"{concordance:.6f}",
        }


def coverage_report(
    references: Iterable[SequenceRecord],
    reads: Iterable[SequenceRecord],
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
    min_score: int = 1,
) -> Iterator[CoverageStats]:
    """Yield :class:`CoverageStats` for each reference, in input order.

    Every read is mapped exactly as in
    :func:`~genome_variant.mapping.map_reads`: local alignment against
    every reference record on both strands, with only the unique winning
    candidate (score at least *min_score*) retained.  Reads without a
    winning candidate contribute nothing.  Each ``=`` or ``X`` CIGAR
    column adds one depth observation at its reference position and each
    ``=`` column adds one concordant observation; ``D`` consumes a
    reference position without depth and ``I`` has none.  Quality values
    are ignored.

    One stats object is produced per reference record even when no read
    maps to it; duplicate reference identifiers are never merged.
    Results are byte-stable for fixed inputs regardless of how the caller
    batches read iteration.

    Scoring parameters and *min_score* have the same meaning, validation
    and error timing as :func:`~genome_variant.mapping.map_reads`:
    invalid parameters raise :class:`ValueError` before either input is
    consumed, an empty reference collection raises
    :class:`MappingReferenceError` before *reads* is consumed, and
    invalid sequence symbols raise
    :class:`~genome_variant.sequence_io.SequenceValidationError` when the
    offending record is reached.
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

    def _iter() -> Iterator[CoverageStats]:
        depth: list[list[int]] = [
            [0] * len(record.sequence) for record in reference_records
        ]
        mapped_counts = [0] * len(reference_records)
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
            mapped_counts[ref_index] += 1
            reference_position = alignment.reference_start
            record_depth = depth[ref_index]
            for op in _expand_cigar(alignment.cigar):
                if op == "=":
                    record_depth[reference_position] += 1
                    concordant[ref_index] += 1
                    reference_position += 1
                elif op == "X":
                    record_depth[reference_position] += 1
                    reference_position += 1
                elif op == "D":
                    reference_position += 1
                # "I" consumes no reference position.

        for ref_index, record in enumerate(reference_records):
            record_depth = depth[ref_index]
            yield CoverageStats(
                reference_record=ref_index,
                reference=record.identifier,
                length=len(record.sequence),
                mapped_reads=mapped_counts[ref_index],
                observed_bases=sum(record_depth),
                covered_bases=sum(1 for value in record_depth if value),
                concordant_bases=concordant[ref_index],
            )

    return _iter()
