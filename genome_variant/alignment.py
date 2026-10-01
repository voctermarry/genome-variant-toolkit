"""Pairwise sequence alignment with affine gap penalties.

Public API:

- :class:`PairwiseAlignment` — the result of aligning two records:
  scoring mode, score, zero-based half-open coordinates, CIGAR string and
  gapped alignment strings.
- :func:`align_pair` — align two :class:`~genome_variant.sequence_io.
  SequenceRecord` objects using global (Needleman–Wunsch style) or local
  (Smith–Waterman style) dynamic programming with affine gap costs.

Only exactly equal characters count as a match; every other valid IUPAC
combination is a mismatch.  FASTQ quality values never participate in
scoring.  A gap run of length ``L`` costs ``gap_open + gap_extend *
(L - 1)``.  Ties are resolved deterministically by the smallest start and
end coordinates (reference first, then query), then by the smallest
CIGAR string.
"""

from __future__ import annotations

from dataclasses import dataclass

from .sequence_io import SequenceRecord, _validate_sequence

__all__ = [
    "PairwiseAlignment",
    "align_pair",
]


@dataclass(frozen=True)
class PairwiseAlignment:
    """The result of a pairwise alignment.

    Coordinates are zero-based and half-open.  ``aligned_reference`` and
    ``aligned_query`` use ``"-"`` for gap characters and have equal
    lengths.  ``cigar`` uses ``=`` (match), ``X`` (mismatch), ``I``
    (insertion; consumes the query only) and ``D`` (deletion; consumes the
    reference only); adjacent equal operations are merged.
    """

    reference: str
    query: str
    mode: str
    score: int
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    cigar: str
    aligned_reference: str
    aligned_query: str


def _score_parameter(name: str, value: object, *, positive: bool) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if positive:
        if value <= 0:
            raise ValueError(f"{name} must be a positive integer")
    elif value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return value


def align_pair(
    reference: SequenceRecord,
    query: SequenceRecord,
    mode: str = "global",
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
) -> PairwiseAlignment:
    """Align two sequence records and return the best alignment.

    *mode* is ``"global"`` (align the complete sequences end to end) or
    ``"local"`` (return the best-scoring positive fragment).  *match_score*
    is a positive integer awarded for equal characters; *mismatch_penalty*,
    *gap_open* and *gap_extend* are non-negative integer penalties
    (subtracted for mismatches and gaps).  A run of ``L`` gaps costs
    ``gap_open + gap_extend * (L - 1)``.

    Only exactly equal characters are matches; every other valid IUPAC
    pair is a mismatch.  Quality values are ignored.

    Invalid parameters raise :class:`ValueError` before either record is
    examined; invalid sequence symbols raise
    :class:`~genome_variant.sequence_io.SequenceValidationError`.

    In local mode, when no fragment has a positive score the result has a
    zero score, all four coordinates equal to zero and empty alignment
    strings and CIGAR.  Among equal-score alignments the one with the
    smallest reference start, reference end, query start and query end (in
    that order) is chosen; ties remaining are broken by the
    lexicographically smallest CIGAR string.
    """
    if mode not in ("global", "local"):
        raise ValueError(f"mode must be 'global' or 'local', not {mode!r}")
    match = _score_parameter("match_score", match_score, positive=True)
    mismatch = _score_parameter(
        "mismatch_penalty", mismatch_penalty, positive=False
    )
    open_penalty = _score_parameter("gap_open", gap_open, positive=False)
    extend_penalty = _score_parameter("gap_extend", gap_extend, positive=False)

    _validate_sequence(reference.identifier, reference.sequence)
    _validate_sequence(query.identifier, query.sequence)

    if mode == "global":
        score, cigar = _best_path(
            reference.sequence,
            query.sequence,
            0,
            0,
            len(reference.sequence),
            len(query.sequence),
            match,
            mismatch,
            open_penalty,
            extend_penalty,
            leading_gaps=True,
        )
        return _build_result(
            reference.identifier,
            query.identifier,
            "global",
            reference.sequence,
            query.sequence,
            score,
            cigar,
            0,
            len(reference.sequence),
            0,
            len(query.sequence),
        )

    return _align_local(
        reference.identifier,
        query.identifier,
        reference.sequence,
        query.sequence,
        match,
        mismatch,
        open_penalty,
        extend_penalty,
    )


# --- dynamic programming -------------------------------------------------

# Per-cell path entries are (score, merged_cigar); "larger score, then
# smaller cigar" dominates.  At a fixed grid cell two distinct merged
# CIGARs can never be prefixes of one another (every extra token consumes
# more coordinates), and their first difference always lies in an already
# closed token that later columns never mutate, so keeping one winner per
# cell is safe.


def _append_op(cigar: str, op: str) -> str:
    """Append one column operation to a merged CIGAR string."""
    if cigar and cigar[-1] == op:
        cut = len(cigar) - 1
        while cut > 0 and cigar[cut - 1].isdigit():
            cut -= 1
        count = int(cigar[cut:-1]) + 1
        return cigar[:cut] + str(count) + op
    return cigar + "1" + op


def _better(
    winner: tuple[int, str] | None, candidate: tuple[int, str] | None
) -> tuple[int, str] | None:
    if candidate is None:
        return winner
    if winner is None or candidate[0] > winner[0] or (
        candidate[0] == winner[0] and candidate[1] < winner[1]
    ):
        return candidate
    return winner


def _best_path(
    ref: str,
    qry: str,
    ref_start: int,
    query_start: int,
    ref_end: int,
    query_end: int,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
    *,
    leading_gaps: bool,
    trailing_gaps: bool = True,
) -> tuple[int, str]:
    """Best score and smallest merged CIGAR over paths inside the box.

    The match/mismatch state holds the three classic affine-gap states.
    With *leading_gaps* false the path must start with a diagonal step
    (local alignments only begin at a restart column); with it true pure
    gap runs along the top and left edges are allowed (global mode).
    """
    width = query_end - query_start

    def empty_row() -> list[tuple[int, str] | None]:
        return [None] * (width + 1)

    prev_match = empty_row()
    prev_delete = empty_row()
    prev_insert = empty_row()
    prev_match[0] = (0, "")
    if leading_gaps and width >= 1:
        # Top edge: the first insertion opens a gap from the match state;
        # the rest extend it.  Match-state cells along the edge stay
        # unreachable, so no later run can re-open there.
        seed = prev_match[0]
        assert seed is not None
        prev_insert[1] = (seed[0] - gap_open, _append_op(seed[1], "I"))
        for offset in range(2, width + 1):
            earlier = prev_insert[offset - 1]
            assert earlier is not None
            prev_insert[offset] = (
                earlier[0] - gap_extend,
                _append_op(earlier[1], "I"),
            )

    for i in range(ref_start + 1, ref_end + 1):
        cur_match = empty_row()
        cur_delete = empty_row()
        cur_insert = empty_row()

        if leading_gaps:
            depth = i - ref_start
            if depth == 1:
                cur_delete[0] = (-gap_open, "1D")
            else:
                earlier = prev_delete[0]
                assert earlier is not None
                cur_delete[0] = (
                    earlier[0] - gap_extend,
                    _append_op(earlier[1], "D"),
                )

        for offset in range(1, width + 1):
            j = query_start + offset
            substitution = match if ref[i - 1] == qry[j - 1] else -mismatch
            op = "=" if substitution == match and ref[i - 1] == qry[j - 1] else "X"

            diagonal: tuple[int, str] | None = None
            for predecessor in (
                prev_match[offset - 1],
                prev_delete[offset - 1],
                prev_insert[offset - 1],
            ):
                if predecessor is None:
                    continue
                diagonal = _better(
                    diagonal,
                    (predecessor[0] + substitution, _append_op(predecessor[1], op)),
                )
            cur_match[offset] = diagonal

            deletion: tuple[int, str] | None = None
            for predecessor, cost in (
                (prev_match[offset], -gap_open),
                (prev_delete[offset], -gap_extend),
                (prev_insert[offset], -gap_open),
            ):
                if predecessor is None:
                    continue
                deletion = _better(
                    deletion,
                    (predecessor[0] + cost, _append_op(predecessor[1], "D")),
                )
            cur_delete[offset] = deletion

            insertion: tuple[int, str] | None = None
            for predecessor, cost in (
                (cur_match[offset - 1], -gap_open),
                (cur_delete[offset - 1], -gap_open),
                (cur_insert[offset - 1], -gap_extend),
            ):
                if predecessor is None:
                    continue
                insertion = _better(
                    insertion,
                    (predecessor[0] + cost, _append_op(predecessor[1], "I")),
                )
            cur_insert[offset] = insertion

        prev_match, prev_delete, prev_insert = cur_match, cur_delete, cur_insert

    ending: tuple[int, str] | None
    if trailing_gaps:
        ending = _better(
            _better(prev_match[width], prev_delete[width]), prev_insert[width]
        )
    else:
        # Local alignments must finish with a match/mismatch column; gap
        # states past the final diagonal are not valid endpoints.
        ending = prev_match[width]
    assert ending is not None
    return ending


def _local_scan(
    ref: str,
    qry: str,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
) -> tuple[int, int, int, int] | None:
    """Scalar affine Smith-Waterman scan.

    Returns ``(ref_start, ref_end, query_start, query_end)`` of the
    preferred positive-score endpoint (highest score, then the four
    coordinates in spec order) or ``None`` when no positive fragment
    exists.  Only the match/mismatch state may restart an alignment.
    """
    m, n = len(ref), len(qry)
    # Each reachable state cell holds (score, ref_start, query_start).
    match_grid: list[list[tuple[int, int, int] | None]] = [
        [None] * (n + 1) for _ in range(m + 1)
    ]
    delete_grid = [[None] * (n + 1) for _ in range(m + 1)]
    insert_grid = [[None] * (n + 1) for _ in range(m + 1)]

    best: tuple[int, int, int, int, int] | None = None  # score + coords

    def earlier_start(
        winner: tuple[int, int, int] | None, candidate: tuple[int, int, int] | None
    ) -> tuple[int, int, int] | None:
        if candidate is None:
            return winner
        if winner is None or candidate[0] > winner[0] or (
            candidate[0] == winner[0] and candidate[1:] < winner[1:]
        ):
            return candidate
        return winner

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            substitution = match if ref[i - 1] == qry[j - 1] else -mismatch

            diagonal: tuple[int, int, int] | None = None
            for grid in (match_grid, delete_grid, insert_grid):
                pred = grid[i - 1][j - 1]
                if pred is not None:
                    diagonal = earlier_start(
                        diagonal,
                        (pred[0] + substitution, pred[1], pred[2]),
                    )
            restart: tuple[int, int, int] = (
                substitution,
                i - 1,
                j - 1,
            )
            match_grid[i][j] = earlier_start(diagonal, restart)

            deletion: tuple[int, int, int] | None = None
            for grid, cost in (
                (match_grid, gap_open),
                (insert_grid, gap_open),
            ):
                pred = grid[i - 1][j]
                if pred is not None:
                    deletion = earlier_start(
                        deletion, (pred[0] - cost, pred[1], pred[2])
                    )
            pred = delete_grid[i - 1][j]
            if pred is not None:
                deletion = earlier_start(
                    deletion, (pred[0] - gap_extend, pred[1], pred[2])
                )
            delete_grid[i][j] = deletion

            insertion: tuple[int, int, int] | None = None
            for grid, cost in (
                (match_grid, gap_open),
                (delete_grid, gap_open),
            ):
                pred = grid[i][j - 1]
                if pred is not None:
                    insertion = earlier_start(
                        insertion, (pred[0] - cost, pred[1], pred[2])
                    )
            pred = insert_grid[i][j - 1]
            if pred is not None:
                insertion = earlier_start(
                    insertion, (pred[0] - gap_extend, pred[1], pred[2])
                )
            insert_grid[i][j] = insertion

            # Only match/mismatch state cells are valid endpoints: a path
            # ending in a gap state is dominated by trimming its trailing
            # gap run (non-positive cost, smaller end coordinates).
            entry = match_grid[i][j]
            if entry is not None and entry[0] > 0:
                key = (entry[0], entry[1], i, entry[2], j)
                if best is None or (
                    -key[0], key[1], key[2], key[3], key[4]
                ) < (-best[0], best[1], best[2], best[3], best[4]):
                    best = key

    if best is None:
        return None
    return best[1], best[2], best[3], best[4]


def _expand_cigar(cigar: str) -> list[str]:
    """Expand a merged CIGAR string into one operation character per column."""
    ops: list[str] = []
    digits = ""
    for char in cigar:
        if char.isdigit():
            digits += char
        else:
            ops.extend([char] * int(digits))
            digits = ""
    return ops


def _build_result(
    reference_id: str,
    query_id: str,
    mode: str,
    ref: str,
    qry: str,
    score: int,
    cigar: str,
    reference_start: int,
    reference_end: int,
    query_start: int,
    query_end: int,
) -> PairwiseAlignment:
    aligned_reference: list[str] = []
    aligned_query: list[str] = []
    i, j = reference_start, query_start
    for op in _expand_cigar(cigar):
        if op in ("=", "X"):
            aligned_reference.append(ref[i])
            aligned_query.append(qry[j])
            i += 1
            j += 1
        elif op == "D":
            aligned_reference.append(ref[i])
            aligned_query.append("-")
            i += 1
        else:  # "I"
            aligned_reference.append("-")
            aligned_query.append(qry[j])
            j += 1

    return PairwiseAlignment(
        reference=reference_id,
        query=query_id,
        mode=mode,
        score=score,
        reference_start=reference_start,
        reference_end=reference_end,
        query_start=query_start,
        query_end=query_end,
        cigar=cigar,
        aligned_reference="".join(aligned_reference),
        aligned_query="".join(aligned_query),
    )


def _align_local(
    reference_id: str,
    query_id: str,
    ref: str,
    qry: str,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
) -> PairwiseAlignment:
    box = _local_scan(ref, qry, match, mismatch, gap_open, gap_extend)
    if box is None:
        return PairwiseAlignment(
            reference=reference_id,
            query=query_id,
            mode="local",
            score=0,
            reference_start=0,
            reference_end=0,
            query_start=0,
            query_end=0,
            cigar="",
            aligned_reference="",
            aligned_query="",
        )

    ref_start, ref_end, query_start, query_end = box
    score, cigar = _best_path(
        ref,
        qry,
        ref_start,
        query_start,
        ref_end,
        query_end,
        match,
        mismatch,
        gap_open,
        gap_extend,
        leading_gaps=False,
        trailing_gaps=False,
    )
    return _build_result(
        reference_id,
        query_id,
        "local",
        ref,
        qry,
        score,
        cigar,
        ref_start,
        ref_end,
        query_start,
        query_end,
    )
