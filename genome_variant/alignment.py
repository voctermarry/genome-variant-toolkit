"""Pairwise sequence alignment with affine gap penalties.

Public API:

- :func:`align_pair` — align two :class:`~genome_variant.sequence_io.SequenceRecord`
  values in global (Needleman–Wunsch) or local (Smith–Waterman) mode.

Scoring uses four integer parameters: *match_score* (positive) for equal
characters, *mismatch_penalty*, *gap_open_penalty* and
*gap_extend_penalty* (all non-negative, subtracted from the score).  A
gap of length L costs ``gap_open + gap_extend * (L - 1)``.  Only
identical characters count as a match; every other IUPAC combination,
ambiguity symbols included, is a mismatch.  FASTQ quality values never
participate in scoring.
"""

from __future__ import annotations

from .sequence_io import SequenceRecord, _validate_sequence

__all__ = [
    "align_pair",
]

# State slots per DP cell: match/mismatch column, insertion (gap in the
# reference, consumes query), deletion (gap in the query, consumes
# reference).
_MATCH = 0
_INSERTION = 1
_DELETION = 2

# A stored cell candidate is
# ``(score, reference_start, query_start, merged_cigar)``.
# All paths meeting in one cell share its end coordinates, so only the
# starts and the CIGAR are needed for tie-breaking there.
_Candidate = tuple[int, int, int, str]


def align_pair(
    reference: SequenceRecord,
    query: SequenceRecord,
    mode: str = "global",
    match_score: int = 2,
    mismatch_penalty: int = 2,
    gap_open_penalty: int = 2,
    gap_extend_penalty: int = 1,
) -> dict[str, object]:
    """Align *query* against *reference*.

    *mode* is ``"global"`` (align both full sequences, Needleman–Wunsch)
    or ``"local"`` (return the best positive-scoring segment,
    Smith–Waterman).  In local mode a run without a positive-scoring
    segment returns score 0 with all four coordinates at 0 and empty
    alignment strings and CIGAR.

    The result is a dict with keys in the fixed order ``reference``,
    ``query``, ``mode``, ``score``, ``reference_start``,
    ``reference_end``, ``query_start``, ``query_end``, ``cigar``,
    ``aligned_reference``, ``aligned_query``.  Coordinates are
    zero-based and half-open; gaps in an alignment string are ``-``;
    the CIGAR merges adjacent identical operations and uses ``=``
    (match), ``X`` (mismatch), ``I`` (insertion, consumes query only)
    and ``D`` (deletion, consumes reference only).

    Ties resolve by the smallest reference_start, query_start,
    reference_end, query_end in that order, then by the lexicographically
    smallest CIGAR.  *match_score* must be a non-boolean positive integer
    and the three penalties non-boolean non-negative integers; invalid
    parameters raise :class:`ValueError` before the records are examined.
    A sequence holding a symbol outside the IUPAC set raises
    :class:`~genome_variant.sequence_io.SequenceValidationError`.
    """
    if mode not in ("global", "local"):
        raise ValueError(f"mode must be 'global' or 'local', not {mode!r}")
    match = _score_value("match_score", match_score, minimum=1)
    mismatch = _score_value("mismatch_penalty", mismatch_penalty, minimum=0)
    gap_open = _score_value("gap_open_penalty", gap_open_penalty, minimum=0)
    gap_extend = _score_value("gap_extend_penalty", gap_extend_penalty, minimum=0)

    _validate_sequence(reference.identifier, reference.sequence)
    _validate_sequence(query.identifier, query.sequence)

    ref = reference.sequence
    qry = query.sequence
    if mode == "global":
        score, cigar, ref_start, query_start = _solve_global(
            ref, qry, match, mismatch, gap_open, gap_extend
        )
    else:
        score, cigar, ref_start, query_start = _solve_local(
            ref, qry, match, mismatch, gap_open, gap_extend
        )

    aligned_reference, aligned_query, ref_end, query_end = _build_strings(
        ref, qry, ref_start, query_start, cigar
    )
    return {
        "reference": reference.identifier,
        "query": query.identifier,
        "mode": mode,
        "score": score,
        "reference_start": ref_start,
        "reference_end": ref_end,
        "query_start": query_start,
        "query_end": query_end,
        "cigar": cigar,
        "aligned_reference": aligned_reference,
        "aligned_query": aligned_query,
    }


def _score_value(name: str, value: object, *, minimum: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if value < minimum:
        if minimum == 1:
            raise ValueError(f"{name} must be a positive integer")
        raise ValueError(f"{name} must be a non-negative integer")
    return value


def _append_operation(cigar: str, operation: str) -> str:
    """Extend a merged CIGAR with one operation, merging a same-op run."""
    if cigar and cigar[-1] == operation:
        split = len(cigar) - 1
        while split > 0 and cigar[split - 1].isdigit():
            split -= 1
        count = int(cigar[split:-1])
        return cigar[:split] + str(count + 1) + operation
    return cigar + "1" + operation


def _extend(candidate: _Candidate, delta: int, operation: str) -> _Candidate:
    """Append one alignment column to a cell candidate."""
    score, ref_start, query_start, cigar = candidate
    return (
        score + delta,
        ref_start,
        query_start,
        _append_operation(cigar, operation),
    )


def _better(candidate: tuple, best: tuple | None) -> bool:
    """Higher score wins; then smaller coordinates in order; then a
    lexicographically smaller CIGAR (the last tuple element)."""
    if best is None:
        return True
    if candidate[0] != best[0]:
        return candidate[0] > best[0]
    for index in range(1, len(candidate) - 1):
        if candidate[index] != best[index]:
            return candidate[index] < best[index]
    return candidate[-1] < best[-1]


def _choose(candidates: list[_Candidate | None]) -> _Candidate | None:
    best: _Candidate | None = None
    for candidate in candidates:
        if candidate is not None and _better(candidate, best):
            best = candidate
    return best


def _empty_row(width: int) -> list[list[_Candidate | None]]:
    return [[None, None, None] for _ in range(width + 1)]


def _column_score(
    ref_char: str, query_char: str, match: int, mismatch: int
) -> int:
    return match if ref_char == query_char else -mismatch


def _gap_candidates(
    same_row_opening: _Candidate | None,
    extending: _Candidate | None,
    other_state_opening: _Candidate | None,
    gap_open: int,
    gap_extend: int,
    operation: str,
) -> list[_Candidate | None]:
    """Build the predecessors of a gap cell.

    A gap opens from a match/mismatch column or straight after a gap of
    the other kind (each switch pays the opening penalty), and extends
    within the same gap state at the extension penalty.  Allowing the
    cross-state switch keeps every legal operation sequence reachable,
    which the CIGAR tie-break relies on when penalties are zero.
    """
    candidates: list[_Candidate | None] = []
    if same_row_opening is not None:
        candidates.append(_extend(same_row_opening, -gap_open, operation))
    if extending is not None:
        candidates.append(_extend(extending, -gap_extend, operation))
    if other_state_opening is not None:
        candidates.append(_extend(other_state_opening, -gap_open, operation))
    return candidates


def _solve_global(
    ref: str,
    qry: str,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
) -> tuple[int, str, int, int]:
    n, m = len(ref), len(qry)

    # Row 0: the empty alignment and leading insertions are reachable.
    previous = _empty_row(m)
    previous[0][_MATCH] = (0, 0, 0, "")
    for j in range(1, m + 1):
        previous[j][_INSERTION] = (
            -(gap_open + gap_extend * (j - 1)),
            0,
            0,
            f"{j}I",
        )

    for i in range(1, n + 1):
        current = _empty_row(m)
        # Column 0: only leading deletions are reachable.
        current[0][_DELETION] = (
            -(gap_open + gap_extend * (i - 1)),
            0,
            0,
            f"{i}D",
        )
        ref_char = ref[i - 1]
        for j in range(1, m + 1):
            query_char = qry[j - 1]
            operation = "=" if ref_char == query_char else "X"
            delta = _column_score(ref_char, query_char, match, mismatch)
            current[j][_MATCH] = _choose(
                [
                    _extend(previous[j - 1][state], delta, operation)
                    for state in (_MATCH, _INSERTION, _DELETION)
                    if previous[j - 1][state] is not None
                ]
            )
            current[j][_INSERTION] = _choose(
                _gap_candidates(
                    current[j - 1][_MATCH],
                    current[j - 1][_INSERTION],
                    current[j - 1][_DELETION],
                    gap_open,
                    gap_extend,
                    "I",
                )
            )
            current[j][_DELETION] = _choose(
                _gap_candidates(
                    previous[j][_MATCH],
                    previous[j][_DELETION],
                    previous[j][_INSERTION],
                    gap_open,
                    gap_extend,
                    "D",
                )
            )
        previous = current

    ending = _choose(
        [previous[m][_MATCH], previous[m][_INSERTION], previous[m][_DELETION]]
    )
    assert ending is not None
    return ending[0], ending[3], 0, 0


def _solve_local(
    ref: str,
    qry: str,
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
) -> tuple[int, str, int, int]:
    n, m = len(ref), len(qry)
    previous = _empty_row(m)
    # Global best layout: (score, ref_start, query_start, ref_end,
    # query_end, cigar) — exactly the user-facing tie-break order.
    best: tuple | None = None

    # Row 0 holds leading insertions.  At every cell a gap may also open
    # fresh from the virtual zero-score boundary: with a zero opening
    # penalty such a leading gap ties the trimmed segment while carrying
    # a smaller start coordinate, which the tie-break is required to
    # prefer.
    for j in range(1, m + 1):
        candidates: list[_Candidate | None] = []
        if previous[j - 1][_INSERTION] is not None:
            candidates.append(
                _extend(previous[j - 1][_INSERTION], -gap_extend, "I")
            )
        candidates.append((-gap_open, 0, j - 1, "1I"))
        previous[j][_INSERTION] = _choose(candidates)

    for i in range(1, n + 1):
        current = _empty_row(m)
        # Column 0 holds leading deletions, with the same fresh-restart
        # choice at every row.
        current[0][_DELETION] = _choose(
            [
                _extend(previous[0][_DELETION], -gap_extend, "D")
                if previous[0][_DELETION] is not None
                else None,
                (-gap_open, i - 1, 0, "1D"),
            ]
        )
        ref_char = ref[i - 1]
        for j in range(1, m + 1):
            query_char = qry[j - 1]
            operation = "=" if ref_char == query_char else "X"
            delta = _column_score(ref_char, query_char, match, mismatch)

            # A match/mismatch column either extends a diagonal
            # predecessor or starts a fresh segment here.  A fresh
            # mismatch seed has score -mismatch and can only reach the
            # global best when that penalty is zero, where it ties the
            # trimmed segment with a smaller start coordinate.
            column_candidates: list[_Candidate | None] = [
                _extend(previous[j - 1][state], delta, operation)
                for state in (_MATCH, _INSERTION, _DELETION)
                if previous[j - 1][state] is not None
            ]
            column_candidates.append((delta, i - 1, j - 1, "1" + operation))
            # Negative or zero intermediates are retained: after more
            # positive columns they may finish positive.
            current[j][_MATCH] = _choose(column_candidates)

            # An insertion opens from a match column, straight after a
            # deletion (paying the opening penalty), or from the virtual
            # zero boundary of the cell.
            insertion_candidates = _gap_candidates(
                current[j - 1][_MATCH],
                current[j - 1][_INSERTION],
                current[j - 1][_DELETION],
                gap_open,
                gap_extend,
                "I",
            )
            insertion_candidates.append((-gap_open, i, j - 1, "1I"))
            current[j][_INSERTION] = _choose(insertion_candidates)

            deletion_candidates = _gap_candidates(
                previous[j][_MATCH],
                previous[j][_DELETION],
                previous[j][_INSERTION],
                gap_open,
                gap_extend,
                "D",
            )
            deletion_candidates.append((-gap_open, i - 1, j, "1D"))
            current[j][_DELETION] = _choose(deletion_candidates)

            # Every path stored at cell (i, j) ends at coordinates
            # (i, j); only strictly positive complete paths compete, and
            # gap-ending paths lose the end-coordinate tie-break against
            # the same alignment trimmed of its zero-cost trailing gaps.
            for state in (_MATCH, _INSERTION, _DELETION):
                candidate = current[j][state]
                if candidate is None or candidate[0] <= 0:
                    continue
                scored = (
                    candidate[0],
                    candidate[1],
                    candidate[2],
                    i,
                    j,
                    candidate[3],
                )
                if _better(scored, best):
                    best = scored
        previous = current

    if best is None:
        return 0, "", 0, 0
    return best[0], best[5], best[1], best[2]


def _build_strings(
    ref: str,
    qry: str,
    ref_start: int,
    query_start: int,
    cigar: str,
) -> tuple[str, str, int, int]:
    """Reconstruct the gapped alignment strings from a merged CIGAR."""
    aligned_reference: list[str] = []
    aligned_query: list[str] = []
    ref_index = ref_start
    query_index = query_start

    run_length = 0
    for char in cigar:
        if char.isdigit():
            run_length = run_length * 10 + int(char)
            continue
        for _ in range(run_length):
            if char in ("=", "X"):
                aligned_reference.append(ref[ref_index])
                aligned_query.append(qry[query_index])
                ref_index += 1
                query_index += 1
            elif char == "I":
                aligned_reference.append("-")
                aligned_query.append(qry[query_index])
                query_index += 1
            else:  # "D"
                aligned_reference.append(ref[ref_index])
                aligned_query.append("-")
                ref_index += 1
        run_length = 0

    return "".join(aligned_reference), "".join(aligned_query), ref_index, query_index
