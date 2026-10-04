"""Canonical interpretation of the columns of a winning alignment.

Internal to the package — no public symbol is defined here.  The
mapping adjudication in :mod:`genome_variant.mapping` selects exactly
one :class:`~genome_variant.alignment.PairwiseAlignment` per read; every
post-mapping entry (coverage statistics, SNV evidence and indel
evidence) must interpret that one alignment with one set of semantics.
This module expands the winning CIGAR once into per-column coordinates,
bases and gap boundaries so the entries never re-walk the gapped
alignment strings independently.

Column semantics, shared by every consumer:

- ``=`` and ``X`` each consume one reference position and one aligned
  query position and contribute one base observation at that reference
  position;
- ``I`` consumes the aligned query only: it has no reference position;
- ``D`` consumes one reference position without adding depth;
- reference coordinates are always zero-based on the reference's
  forward strand;
- query coordinates and query bases are on the aligned strand, i.e. for
  a reverse-strand win they index the reverse-complement read and the
  observed bases are already expressed on the reference's forward
  strand;
- :func:`original_read_position` maps a column back to the original
  read coordinate, so quality values always correspond to the original
  read even though the quality string is never reverse complemented;
- ``boundary`` is the reference cursor at the left edge of a column:
  the reference position consumed by a ``=``/``X``/``D`` column, or the
  insertion point of an ``I`` column (the boundary between zero-based
  reference positions ``boundary - 1`` and ``boundary``).
"""

from __future__ import annotations

from dataclasses import dataclass

from .alignment import PairwiseAlignment, _expand_cigar

__all__: list[str] = []


@dataclass(frozen=True)
class AlignmentColumn:
    """One interpreted CIGAR column of a winning alignment.

    ``op`` is ``"="``, ``"X"``, ``"I"`` or ``"D"``.
    ``reference_position`` is the zero-based reference position, or
    ``None`` for an insertion; ``query_position`` is the zero-based
    position on the aligned strand, or ``None`` for a deletion.
    ``reference_base``/``query_base`` carry the observed characters
    (``None`` on the gapped side) and ``boundary`` the left-edge
    reference cursor.
    """

    op: str
    reference_position: int | None
    query_position: int | None
    reference_base: str | None
    query_base: str | None
    boundary: int


def alignment_columns(alignment: PairwiseAlignment) -> tuple[AlignmentColumn, ...]:
    """Expand *alignment* into one :class:`AlignmentColumn` per CIGAR column.

    The CIGAR is the single source of which sequences each column
    consumes; the gapped alignment strings only supply the observed
    characters.  Coordinates start at the alignment's reference/query
    intervals and advance exactly as the CIGAR consumes them.
    """
    aligned_reference = alignment.aligned_reference
    aligned_query = alignment.aligned_query
    reference_position = alignment.reference_start
    query_position = alignment.query_start

    columns: list[AlignmentColumn] = []
    for column_index, op in enumerate(_expand_cigar(alignment.cigar)):
        # The reference cursor before this column consumes anything:
        # the consumed position for =/X/D, the insertion point for I.
        boundary = reference_position
        reference_base = aligned_reference[column_index]
        query_base = aligned_query[column_index]

        if op == "I":
            columns.append(
                AlignmentColumn(
                    op,
                    None,
                    query_position,
                    None,
                    query_base,
                    boundary,
                )
            )
            query_position += 1
        elif op == "D":
            columns.append(
                AlignmentColumn(
                    op,
                    reference_position,
                    None,
                    reference_base,
                    None,
                    boundary,
                )
            )
            reference_position += 1
        else:  # "=" or "X": one reference-position observation.
            columns.append(
                AlignmentColumn(
                    op,
                    reference_position,
                    query_position,
                    reference_base,
                    query_base,
                    boundary,
                )
            )
            reference_position += 1
            query_position += 1

    return tuple(columns)


def original_read_position(
    column: AlignmentColumn, read_length: int, strand: str
) -> int:
    """Map a query-consuming column back to its original-read coordinate.

    On the forward strand the aligned-strand position *is* the original
    read coordinate.  On the reverse strand the aligned query is the
    reverse complement, so the coordinate is mirrored; the quality
    string is never reverse complemented itself.
    """
    aligned_position = column.query_position
    assert aligned_position is not None
    if strand == "+":
        return aligned_position
    return read_length - 1 - aligned_position
