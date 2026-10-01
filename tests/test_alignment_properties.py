"""Property tests comparing the DP with exhaustive path enumeration."""

from __future__ import annotations

import itertools

import pytest

from genome_variant.alignment import align_pair
from genome_variant.sequence_io import SequenceRecord

PARAM_SETS = [
    (1, 0, 0, 0),
    (2, 1, 2, 1),
    (2, 3, 5, 2),
]
BASES = "ACG"


def _merge(ops: list[str]) -> str:
    out: list[list] = []
    for op in ops:
        if out and out[-1][1] == op:
            out[-1][0] += 1
        else:
            out.append([1, op])
    return "".join(f"{c}{o}" for c, o in out)


def _score(ops: list[str], match: int, mismatch: int, gap_open: int, gap_extend: int) -> int:
    total = 0
    previous = None
    for op in ops:
        if op in ("=", "X"):
            total += match if op == "=" else -mismatch
        else:
            total -= gap_extend if previous == op else gap_open
        previous = op
    return total


def _paths(m: int, n: int) -> list[tuple[str, ...]]:
    found: list[tuple[str, ...]] = []

    def walk(i: int, j: int, steps: list[str]) -> None:
        if i == m and j == n:
            found.append(tuple(steps))
            return
        if i < m and j < n:
            walk(i + 1, j + 1, steps + ["S"])
        if i < m:
            walk(i + 1, j, steps + ["D"])
        if j < n:
            walk(i, j + 1, steps + ["I"])

    walk(0, 0, [])
    return found


def _columns(path: tuple[str, ...], ref: str, qry: str) -> list[str]:
    ops: list[str] = []
    i = j = 0
    for step in path:
        if step == "S":
            ops.append("=" if ref[i] == qry[j] else "X")
            i += 1
            j += 1
        elif step == "D":
            ops.append("D")
            i += 1
        else:
            ops.append("I")
            j += 1
    return ops


def _best_global(ref: str, qry: str, params: tuple[int, ...]):
    best = None
    for path in _paths(len(ref), len(qry)):
        ops = _columns(path, ref, qry)
        cigar = _merge(ops)
        score = _score(ops, *params)
        if best is None or (-score, cigar) < (-best[0], best[1]):
            best = (score, cigar)
    return best


def _best_local(ref: str, qry: str, params: tuple[int, ...]):
    best = None
    for a in range(len(ref)):
        for b in range(a + 1, len(ref) + 1):
            for c in range(len(qry)):
                for d in range(c + 1, len(qry) + 1):
                    rseg, qseg = ref[a:b], qry[c:d]
                    for path in _paths(len(rseg), len(qseg)):
                        if path[0] != "S":
                            continue
                        ops = _columns(path, rseg, qseg)
                        score = _score(ops, *params)
                        if score <= 0:
                            continue
                        cigar = _merge(ops)
                        key = (-score, a, b, c, d, cigar)
                        if best is None or key < best:
                            best = key
    return best


def _all_sequences(length: int):
    for letters in itertools.product(BASES, repeat=length):
        yield "".join(letters)


def _cases():
    for m in range(4):
        for n in range(4):
            for ref in _all_sequences(m):
                for qry in _all_sequences(n):
                    for params in PARAM_SETS:
                        yield ref, qry, params


@pytest.mark.parametrize("ref,qry,params", list(_cases()))
def test_global_matches_exhaustive_enumeration(ref, qry, params):
    result = align_pair(
        SequenceRecord("r", ref),
        SequenceRecord("q", qry),
        mode="global",
        match_score=params[0],
        mismatch_penalty=params[1],
        gap_open=params[2],
        gap_extend=params[3],
    )
    score, cigar = _best_global(ref, qry, params)
    assert result.score == score
    assert result.cigar == cigar
    assert (result.reference_start, result.reference_end) == (0, len(ref))
    assert (result.query_start, result.query_end) == (0, len(qry))


@pytest.mark.parametrize("ref,qry,params", list(_cases()))
def test_local_matches_exhaustive_enumeration(ref, qry, params):
    result = align_pair(
        SequenceRecord("r", ref),
        SequenceRecord("q", qry),
        mode="local",
        match_score=params[0],
        mismatch_penalty=params[1],
        gap_open=params[2],
        gap_extend=params[3],
    )
    expected = _best_local(ref, qry, params)
    if expected is None:
        assert result.score == 0
        assert result.cigar == ""
        assert (
            result.reference_start,
            result.reference_end,
            result.query_start,
            result.query_end,
        ) == (0, 0, 0, 0)
    else:
        neg_score, a, b, c, d, cigar = expected
        assert result.score == -neg_score
        assert result.cigar == cigar
        assert (
            result.reference_start,
            result.reference_end,
            result.query_start,
            result.query_end,
        ) == (a, b, c, d)
