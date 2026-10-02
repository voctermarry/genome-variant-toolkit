"""Single-sample SNV calling from reads mapped against a reference.

Public API:

- :class:`VariantCallingError` — the reference or sample configuration is
  unusable for calling (empty reference, duplicate reference identifier
  or invalid sample name).
- :func:`call_variants` — map FASTQ reads with the same local-alignment
  scoring, strand search and candidate adjudication as
  :func:`~genome_variant.mapping.map_reads`, pile up the quality-filtered
  substitution evidence and return a VCFv4.2 document.

Only substitutions between the bases ``A``, ``C``, ``G`` and ``T`` are
called.  Insertions, deletions, unaligned portions, unmapped reads and
ambiguous bases never contribute evidence.  Quality values never
participate in mapping; they only filter base evidence.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

from .mapping import _best_candidate
from .quality import ReadQualityError
from .sequence_io import SequenceRecord
from .vcf import VcfFile, VcfHeader, VcfRecord

__all__ = [
    "VariantCallingError",
    "call_variants",
]

_ACGT = ("A", "C", "G", "T")
_ACGT_SET = frozenset(_ACGT)
_PHRED_MIN = 0
_PHRED_MAX = 93

_META_LINES = (
    "##fileformat=VCFv4.2",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
)


class VariantCallingError(ValueError):
    """The reference or sample configuration is unusable for calling."""


def _quality_threshold(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if not _PHRED_MIN <= value <= _PHRED_MAX:
        raise ValueError(f"{name} must be between {_PHRED_MIN} and {_PHRED_MAX}")
    return value


def _positive_threshold(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a non-boolean integer, not {kind}")
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _fraction_threshold(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        kind = type(value).__name__
        raise ValueError(f"{name} must be a finite number, not {kind}")
    fraction = float(value)
    if not math.isfinite(fraction):
        raise ValueError(f"{name} must be finite")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1")
    return fraction


def _validate_sample_name(sample_name: object) -> str:
    if not isinstance(sample_name, str):
        raise VariantCallingError("sample name must be a string")
    if not sample_name or sample_name.isspace() or "\t" in sample_name:
        raise VariantCallingError(
            "sample name must not be empty, whitespace-only or contain tabs"
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
    sample_name: str = "SAMPLE",
) -> VcfFile:
    """Call single-sample SNVs and return them as a :class:`~genome_variant.vcf.VcfFile`.

    Each read is mapped exactly as in :func:`~genome_variant.mapping.map_reads`
    (same local-alignment scoring, both strands, same candidate
    adjudication); unmapped reads contribute nothing.  Every mapped read
    contributes at most one observation per covered reference position.
    A base is admitted as evidence only when its Phred quality is at least
    *min_base_quality* and both the reference base and the observed base
    are among ``A``, ``C``, ``G`` and ``T``; on the reverse strand the
    quality and base are taken from the corresponding original-read
    position.  Insertions, deletions and unaligned alignment columns are
    ignored.

    ``DP`` at a position is the total number of admitted ``ACGT``
    observations.  The single ALT is the non-reference base with the
    highest count (ties resolved in ``A``, ``C``, ``G``, ``T`` order); a
    record is emitted only when its count is at least *min_alt_count* and
    its share ``AC/DP`` is at least *min_alt_fraction*.  ``GT`` is
    ``1/1`` when that share is at least *homozygous_fraction*, otherwise
    ``0/1``.  Records follow reference input order and ascending POS.

    *min_base_quality* is a non-boolean integer in 0-93 (default 20);
    *min_alt_count* a non-boolean positive integer (default 2);
    *min_alt_fraction* and *homozygous_fraction* finite numbers in 0-1
    (defaults 0.2 and 0.8), with the homozygous fraction no smaller than
    the minimum ALT fraction.  Invalid thresholds raise :class:`ValueError`
    before either input is consumed.  Reads must carry quality values; a
    read without them or with a quality/sequence length mismatch raises
    :class:`~genome_variant.quality.ReadQualityError` when that read is
    reached.  An empty reference collection, duplicate reference
    identifiers or an empty/whitespace-only/tab-containing *sample_name*
    raise :class:`VariantCallingError`.
    """
    base_quality = _quality_threshold("min_base_quality", min_base_quality)
    alt_count = _positive_threshold("min_alt_count", min_alt_count)
    alt_fraction = _fraction_threshold("min_alt_fraction", min_alt_fraction)
    homo_fraction = _fraction_threshold(
        "homozygous_fraction", homozygous_fraction
    )
    if homo_fraction < alt_fraction:
        raise ValueError(
            "homozygous_fraction must not be smaller than min_alt_fraction"
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

    # Mapping uses the same local-alignment scoring and candidate
    # adjudication as map-reads; callers cannot tune them here.
    match, mismatch = 2, 3
    open_penalty, extend_penalty, min_score = 5, 2, 1

    # One counter table per reference record: position -> [A, C, G, T].
    counts_per_reference: list[dict[int, list[int]]] = [
        {} for _ in reference_records
    ]

    for read in reads:
        quality = read.quality
        if quality is None:
            raise ReadQualityError(
                f"record {read.identifier!r}: has no quality values"
            )
        if len(quality) != len(read.sequence):
            raise ReadQualityError(
                f"record {read.identifier!r}: quality length {len(quality)} "
                f"does not match sequence length {len(read.sequence)}"
            )
        for score in quality:
            if not isinstance(score, int) or not _PHRED_MIN <= score <= _PHRED_MAX:
                raise ReadQualityError(
                    f"record {read.identifier!r}: quality score {score} is "
                    "outside the Phred+33 range 0-93"
                )

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
        ref_index, alignment, strand = winner
        _collect_evidence(
            read,
            quality,
            reference_records[ref_index].sequence,
            alignment,
            strand,
            base_quality,
            counts_per_reference[ref_index],
        )

    records = []
    for ref_index, reference in enumerate(reference_records):
        chromosome = reference.identifier
        reference_sequence = reference.sequence
        for position in sorted(counts_per_reference[ref_index]):
            counts = counts_per_reference[ref_index][position]
            record = _call_position(
                chromosome,
                reference_sequence,
                position,
                counts,
                alt_count,
                alt_fraction,
                homo_fraction,
            )
            if record is not None:
                records.append(record)

    header = VcfHeader(_META_LINES, (sample,))
    return VcfFile(header, tuple(records))


def _collect_evidence(
    read: SequenceRecord,
    quality: tuple[int, ...],
    reference_sequence: str,
    alignment,
    strand: str,
    min_base_quality: int,
    counts: dict[int, list[int]],
) -> None:
    """Accumulate admitted substitution observations for one read."""
    read_length = len(read.sequence)
    aligned_reference = alignment.aligned_reference
    aligned_query = alignment.aligned_query
    query_start = alignment.query_start
    reference_position = alignment.reference_start
    query_column = 0

    for column in range(len(aligned_reference)):
        reference_char = aligned_reference[column]
        # The aligned query is oriented against the reference; on the
        # reverse strand it is a base of the reverse-complement read, so
        # it is already expressed on the reference's forward strand.
        observed = aligned_query[column]
        consumes_reference = reference_char != "-"
        consumes_query = observed != "-"

        if consumes_reference and consumes_query:
            if strand == "+":
                quality_position = query_start + query_column
            else:
                # Only the quality is mapped back to the original read:
                # the quality string is never reverse complemented.
                quality_position = read_length - 1 - (query_start + query_column)

            if quality[quality_position] >= min_base_quality:
                reference_base = reference_sequence[reference_position]
                if observed in _ACGT_SET and reference_base in _ACGT_SET:
                    table = counts.setdefault(reference_position, [0, 0, 0, 0])
                    table[_ACGT.index(observed)] += 1

        if consumes_reference:
            reference_position += 1
        if consumes_query:
            query_column += 1


def _call_position(
    chromosome: str,
    reference_sequence: str,
    position: int,
    counts: list[int],
    min_alt_count: int,
    min_alt_fraction: float,
    homo_fraction: float,
) -> VcfRecord | None:
    """Build the VCF record for one piled-up position, or ``None``."""
    depth = sum(counts)
    if depth == 0:
        return None

    reference_base = reference_sequence[position]
    reference_index = _ACGT.index(reference_base)

    # Highest non-reference count wins; iterating A, C, G, T in order and
    # replacing only on a strict increase resolves ties by that order.
    alt_index: int | None = None
    alt_observations = 0
    for base_index, observed in enumerate(counts):
        if base_index == reference_index:
            continue
        if observed > alt_observations:
            alt_observations = observed
            alt_index = base_index
    if alt_index is None or alt_observations < min_alt_count:
        return None

    fraction = alt_observations / depth
    if fraction < min_alt_fraction:
        return None

    genotype = "1/1" if fraction >= homo_fraction else "0/1"
    reference_observations = counts[reference_index]
    info = (
        f"DP={depth};AC={alt_observations};AF={fraction:.6f}"
    )
    sample = (
        f"{genotype}:{depth}:{reference_observations},{alt_observations}"
    )
    return VcfRecord(
        chrom=chromosome,
        pos=position + 1,
        id=".",
        ref=reference_base,
        alt=(_ACGT[alt_index],),
        qual=".",
        filter="PASS",
        info=info,
        format_text="GT:DP:AD",
        sample_text=(sample,),
    )
