"""Single-sample SNV calling from reads mapped against a reference.

Public API:

- :class:`VariantCallingError` — the reference collection was empty,
  contained duplicate identifiers, or the sample name was unusable.
- :func:`call_variants` — map reads exactly as
  :func:`~genome_variant.mapping.map_reads` does, pile up quality-filtered
  base evidence and emit single-nucleotide variants as a VCFv4.2 document.

Only substitutions between ``A``, ``C``, ``G`` and ``T`` are called.
Insertions, deletions, unaligned read portions, unmapped reads and
ambiguous (IUPAC) bases never contribute evidence.  Quality values never
participate in alignment scoring; they only filter base evidence.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

from .alignment import _score_parameter
from .mapping import _best_candidate
from .quality import ReadQualityError
from .sequence_io import SequenceRecord
from .vcf import VcfFile, VcfHeader, VcfRecord

__all__ = [
    "VariantCallingError",
    "call_variants",
]

_BASES = ("A", "C", "G", "T")
_PHRED_MIN = 0
_PHRED_MAX = 93
_DEFAULT_SAMPLE = "SAMPLE"

# VCFv4.2 meta-information describing the produced fields.
_META_LINES = (
    '##fileformat=VCFv4.2',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele fraction">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth after quality filtering">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
)


class VariantCallingError(ValueError):
    """The reference collection was empty, identifiers were duplicated, or
    the sample name was empty, whitespace-only or contained a tab.
    """


def _quality_threshold(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if not _PHRED_MIN <= value <= _PHRED_MAX:
        raise ValueError(f"{name} must be between {_PHRED_MIN} and {_PHRED_MAX}")
    return value


def _positive_count(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _ratio(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a finite number, not {kind}")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    if not 0.0 <= result <= 1.0:
        raise ValueError(f"{name} must be between 0.0 and 1.0")
    return result


def _validate_sample_name(sample_name: object) -> str:
    if not isinstance(sample_name, str):
        kind = type(sample_name).__name__
        raise VariantCallingError(
            f"sample_name must be a string, not {kind}"
        )
    if not sample_name or sample_name.isspace() or "\t" in sample_name:
        raise VariantCallingError(
            "sample name must not be empty, whitespace-only or contain a tab"
        )
    return sample_name


def call_variants(
    references: Iterable[SequenceRecord],
    reads: Iterable[SequenceRecord],
    *,
    min_base_quality: int = 20,
    min_alt_count: int = 2,
    min_alt_fraction: float = 0.2,
    homozygous_fraction: float = 0.8,
    sample_name: str = _DEFAULT_SAMPLE,
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
    min_score: int = 1,
) -> VcfFile:
    """Call single-sample SNVs and return a readable VCFv4.2 document.

    Reads are mapped with the same local-alignment scoring, forward/reverse
    strand search and candidate adjudication as
    :func:`~genome_variant.mapping.map_reads` (the scoring keyword
    arguments mirror it); quality values never affect mapping.

    Each mapped read contributes at most one observation per reference
    position, and only for columns where both the reference and the read
    base are among ``A``, ``C``, ``G`` and ``T``; insertions, deletions,
    unaligned portions, unmapped reads and ambiguous bases are ignored.
    On the reverse strand both the base and its Phred quality are taken
    from the original read position.  An observation counts only when its
    Phred quality is at least *min_base_quality*.

    ``DP`` is the number of passing ``ACGT`` observations at a position.
    The most frequent non-reference base is the sole ALT (ties resolved in
    ``A, C, G, T`` order); a record is emitted only when the ALT count is
    at least *min_alt_count* and ``AC/DP`` is at least *min_alt_fraction*.
    The genotype is ``1/1`` when the ALT fraction reaches
    *homozygous_fraction*, otherwise ``0/1``.

    Thresholds are validated before either input is consumed:
    *min_base_quality* is an integer from 0 to 93, *min_alt_count* a
    positive integer, *min_alt_fraction* and *homozygous_fraction* finite
    numbers between 0 and 1 with the latter no smaller than the former, and
    *sample_name* a non-empty, non-whitespace-only string without tabs.
    An empty reference collection or duplicate reference identifiers raise
    :class:`VariantCallingError`; reads lacking quality values or with a
    quality length mismatch raise :class:`ReadQualityError`, and malformed
    sequences raise the usual sequence exceptions, when reached.
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

    base_quality = _quality_threshold("min_base_quality", min_base_quality)
    alt_count = _positive_count("min_alt_count", min_alt_count)
    alt_fraction = _ratio("min_alt_fraction", min_alt_fraction)
    hom_fraction = _ratio("homozygous_fraction", homozygous_fraction)
    if hom_fraction < alt_fraction:
        raise ValueError(
            "homozygous_fraction must be at least min_alt_fraction"
        )
    sample = _validate_sample_name(sample_name)

    reference_records = list(references)
    if not reference_records:
        raise VariantCallingError("at least one reference record is required")
    seen: set[str] = set()
    for reference in reference_records:
        if reference.identifier in seen:
            raise VariantCallingError(
                f"duplicate reference identifier {reference.identifier!r}"
            )
        seen.add(reference.identifier)

    # One position -> A/C/G/T count table per reference record.  A read
    # touches at most one column per position, so plain counters already
    # enforce the once-per-position rule.
    counts: list[dict[int, dict[str, int]]] = [
        {} for _ in reference_records
    ]

    for read in reads:
        _pile_up_one(
            read,
            reference_records,
            counts,
            match,
            mismatch,
            open_penalty,
            extend_penalty,
            min_score,
            base_quality,
        )

    records = []
    for ref_index, reference in enumerate(reference_records):
        chromosome = reference.identifier
        for pos_zero in sorted(counts[ref_index]):
            table = counts[ref_index][pos_zero]
            depth = sum(table.get(base, 0) for base in _BASES)
            reference_base = reference.sequence[pos_zero]
            if reference_base not in _BASES:
                # Every recorded observation matches a reference ACGT
                # column, so no non-reference call is possible here.
                continue
            alt_base = ""
            alt_total = 0
            for base in _BASES:
                if base == reference_base:
                    continue
                if table.get(base, 0) > alt_total:
                    alt_total = table[base]
                    alt_base = base
            if alt_total < alt_count:
                continue
            if depth == 0 or alt_total / depth < alt_fraction:
                continue
            fraction = alt_total / depth
            genotype = "1/1" if fraction >= hom_fraction else "0/1"
            ref_depth = table.get(reference_base, 0)
            records.append(
                VcfRecord(
                    chrom=chromosome,
                    pos=pos_zero + 1,
                    id=".",
                    ref=reference_base,
                    alt=(alt_base,),
                    qual=".",
                    filter="PASS",
                    info=(
                        f"DP={depth};AC={alt_total};"
                        f"AF={fraction:.6f}"
                    ),
                    format_text="GT:DP:AD",
                    sample_text=(
                        f"{genotype}:{depth}:{ref_depth},{alt_total}",
                    ),
                )
            )

    header = VcfHeader(_META_LINES, (sample,))
    return VcfFile(header, tuple(records))


def _pile_up_one(
    read: SequenceRecord,
    references: list[SequenceRecord],
    counts: list[dict[int, dict[str, int]]],
    match: int,
    mismatch: int,
    gap_open: int,
    gap_extend: int,
    min_score: int,
    min_base_quality: int,
) -> None:
    """Map one read and add its quality-filtered observations."""
    quality = read.quality
    if quality is None:
        raise ReadQualityError(
            f"record {read.identifier!r}: has no quality values"
        )
    if len(quality) != len(read.sequence):
        raise ReadQualityError(
            f"record {read.identifier!r}: quality length {len(quality)} does "
            f"not match sequence length {len(read.sequence)}"
        )
    for score in quality:
        if not isinstance(score, int) or not _PHRED_MIN <= score <= _PHRED_MAX:
            raise ReadQualityError(
                f"record {read.identifier!r}: quality score {score!r} is "
                "outside the Phred+33 range 0-93"
            )

    winner = _best_candidate(
        read,
        references,
        match,
        mismatch,
        gap_open,
        gap_extend,
        min_score,
    )
    if winner is None:
        return
    alignment, ref_index, strand = winner
    table = counts[ref_index]

    # Walk the gapped strings in their native orientation against the
    # reference.  On the reverse strand the aligned bases already face the
    # forward reference (they come from the reverse-complemented read), so
    # only the quality index is translated back to the original read:
    # RC-read column j corresponds to original index L-1-j.
    aligned_reference = alignment.aligned_reference
    aligned_query = alignment.aligned_query
    read_length = len(read.sequence)

    ref_position = alignment.reference_start
    if strand == "+":
        query_position = alignment.query_start
        query_step = 1
    else:
        query_position = read_length - 1 - alignment.query_start
        query_step = -1

    for ref_char, query_char in zip(aligned_reference, aligned_query):
        consumes_reference = ref_char != "-"
        consumes_query = query_char != "-"

        if consumes_reference and consumes_query:
            if (
                quality[query_position] >= min_base_quality
                and query_char in _BASES
                and ref_char in _BASES
            ):
                position_table = table.setdefault(ref_position, {})
                position_table[query_char] = (
                    position_table.get(query_char, 0) + 1
                )
            ref_position += 1
            query_position += query_step
        elif consumes_reference:
            ref_position += 1
        elif consumes_query:
            query_position += query_step
