"""Shared interpretation of one winning pairwise alignment.

Internal to the package — the public entry points remain
:func:`genome_variant.mapping.map_reads`,
:func:`genome_variant.coverage.coverage_report` and
:func:`genome_variant.calling.call_variants`.  Every post-alignment
consumer (mapping output, coverage statistics and variant evidence)
asks this module for the single interpretation of a winning
alignment, so one winning alignment can never be explained twice in
slightly different ways.

A :class:`WinningAlignment` owns the adjudicator result (reference
record index, strand and the pairwise alignment of the reference with
the strand-aligned read) together with the read length.  Its
:attr:`WinningAlignment.columns` sequence is the one place that walks
the CIGAR/gapped strings once:

- each column is one CIGAR operation (``=``/``X``/``I``/``D``);
- ``=`` and ``X`` columns each contribute one observation at a
  reference position; ``I`` columns never consume a reference
  position and ``D`` columns consume one without adding depth;
- query coordinates are always stated on the **original** read; on
  the reverse strand the aligned (reverse-complement) columns are
  mapped back, so quality values keep indexing the original quality
  string exactly;
- observed bases are always expressed on the reference's forward
  strand — on the reverse strand the aligned base is the
  reverse-complement read's base, which is already the
  forward-strand observation;
- gap boundaries follow the CIGAR: an ``I``/``D`` run is a maximal
  contiguous run of columns gapped on the same side, with the
  insertion point/deletion interval taken from the next unconsumed
  reference position.

The module never re-aligns, never re-adjudicates and never changes
tie-breaking; it only gives the three post-alignment features one
identical set of reference coordinates, query coordinates, strand,
base observations and gap semantics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .alignment import PairwiseAlignment

__all__: list[str] = []


@dataclass(frozen=True)
class AlignmentColumn:
    """One CIGAR column of a winning alignment.

    ``op`` is ``"="``, ``"X"`` (aligned base pair), ``"I"`` (a query
    base with no reference position) or ``"D"`` (a reference base the
    read skips).  ``reference_position`` is the zero-based reference
    coordinate, or ``None`` on an insertion.  ``boundary_position`` is
    the reference coordinate at the column for every column: the
    consumed coordinate on an ``=``/``X``/``D`` column and the
    insertion point (the gap between ``boundary_position - 1`` and
    ``boundary_position``) on an ``I`` column.  ``original_query_index``
    is the zero-based coordinate on the original read, or ``None`` on
    a deletion; quality values index the original read there even on
    the reverse strand.  ``reference_base``/``observed_base`` hold the
    reference character and the observed character (``None`` on the
    side that is gapped); observed bases are always oriented on the
    reference forward strand.
    """

    op: str
    reference_position: int | None
    boundary_position: int
    original_query_index: int | None
    reference_base: str | None
    observed_base: str | None

    @property
    def is_aligned_pair(self) -> bool:
        """Whether the column is a depth-contributing ``=``/``X`` pair."""
        return self.op in ("=", "X")

    @property
    def is_insertion(self) -> bool:
        """Whether the column is a reference-gap (``I``) column."""
        return self.op == "I"

    @property
    def is_deletion(self) -> bool:
        """Whether the column is a query-gap (``D``) column."""
        return self.op == "D"


@dataclass(frozen=True)
class WinningAlignment:
    """The adjudicator's choice interpreted on the original read.

    Instances are produced by :func:`genome_variant.mapping._best_candidate`.
    ``reference_record`` is the winning reference's input index,
    ``strand`` ``"+"`` or ``"-"``, ``read_length`` the length of the
    original read and ``alignment`` the pairwise alignment of the
    reference with the read on the chosen strand (the read's reverse
    complement on a reverse-strand win).
    """

    reference_record: int
    strand: str
    read_length: int
    alignment: "PairwiseAlignment"

    def original_query_interval(self) -> tuple[int, int]:
        """Return the aligned interval as zero-based half-open original-read coords."""
        aligned_start = self.alignment.query_start
        aligned_end = self.alignment.query_end
        if self.strand == "+":
            return aligned_start, aligned_end
        # Convert the interval on the reverse-complement read back to
        # zero-based half-open coordinates on the original read.
        return self.read_length - aligned_end, self.read_length - aligned_start

    def columns(self) -> tuple[AlignmentColumn, ...]:
        """Interpret the winning alignment once, column by column.

        Coordinates and bases follow the module semantics: original-read
        coordinates with reverse-strand mapping, forward-strand observed
        bases and the reference position as ``None`` on insertions.
        """
        alignment = self.alignment
        read_length = self.read_length
        reverse = self.strand == "-"
        aligned_reference = alignment.aligned_reference
        aligned_query = alignment.aligned_query

        columns: list[AlignmentColumn] = []
        reference_position = alignment.reference_start
        aligned_query_index = alignment.query_start
        for column in range(len(aligned_reference)):
            # Reference coordinate at this column before consuming it:
            # the gap sits immediately left of this position.
            boundary_position = reference_position
            reference_char = aligned_reference[column]
            query_char = aligned_query[column]
            insertion = reference_char == "-"
            deletion = query_char == "-"

            if insertion:
                op = "I"
                ref_pos: int | None = None
                ref_base: str | None = None
            else:
                ref_pos = reference_position
                ref_base = reference_char
                if deletion:
                    op = "D"
                else:
                    op = "=" if reference_char == query_char else "X"

            if deletion:
                original_index: int | None = None
                observed: str | None = None
            else:
                if reverse:
                    original_index = read_length - 1 - aligned_query_index
                else:
                    original_index = aligned_query_index
                observed = query_char

            columns.append(
                AlignmentColumn(
                    op,
                    ref_pos,
                    boundary_position,
                    original_index,
                    ref_base,
                    observed,
                )
            )

            if not insertion:
                reference_position += 1
            if not deletion:
                aligned_query_index += 1

        return tuple(columns)
