"""Single-sample SNV and optional short-indel calling from mapped reads.

Public API:

- :class:`VariantCallingError` — the reference or sample configuration is
  unusable for calling (empty reference, duplicate reference identifier
  or invalid sample name).
- :func:`call_variants` — map FASTQ reads with the same local-alignment
  scoring, strand search and candidate adjudication as
  :func:`~genome_variant.mapping.map_reads`, pile up the quality-filtered
  substitution evidence and return a VCFv4.2 document.  Short insertions
  and deletions are called additionally when *call_indels* is enabled.

By default only substitutions between the bases ``A``, ``C``, ``G`` and
``T`` are called.  Insertions, deletions, unaligned portions, unmapped
reads and ambiguous bases never contribute substitution evidence.
Quality values never participate in mapping; they only filter base
evidence.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

from .mapping import _best_candidate
from .quality import ReadQualityError
from .sequence_io import SequenceRecord
from .vcf import VcfFile, VcfHeader, VcfRecord, normalize_record

__all__ = [
    "VariantCallingError",
    "call_variants",
]

_ACGT = ("A", "C", "G", "T")
_ACGT_SET = frozenset(_ACGT)
_PHRED_MIN = 0
_PHRED_MAX = 93

# One read's normalized indel observations: (pos, ref, alt) mapped to the
# reference gaps the event touches and whether its evidence passes the
# quality floor.
_IndelEventMap = dict[tuple[int, str, str], tuple[frozenset[int], bool]]
# Per mapped read: reference index, aligned-base qualities, indel events.
_IndelObservation = tuple[int, dict[int, int], _IndelEventMap]

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
    call_indels: bool = False,
    max_indel_length: int = 50,
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
    ignored for SNV evidence.

    ``DP`` at a position is the total number of admitted ``ACGT``
    observations.  The single ALT is the non-reference base with the
    highest count (ties resolved in ``A``, ``C``, ``G``, ``T`` order); a
    record is emitted only when its count is at least *min_alt_count* and
    its share ``AC/DP`` is at least *min_alt_fraction*.  ``GT`` is
    ``1/1`` when that share is at least *homozygous_fraction*, otherwise
    ``0/1``.  Records follow reference input order and ascending POS.

    When *call_indels* is true, short insertions and deletions are called
    in addition from the contiguous ``I`` and ``D`` CIGAR segments of each
    winning alignment.  An insertion is reported as the left reference
    base plus the inserted sequence, a deletion as the left anchor base
    plus the deleted reference fragment.  Events without a usable left
    anchor, without aligned read bases on both sides, with any relevant
    reference or inserted base outside ``ACGT``, or longer than
    *max_indel_length* are ignored.  Insertion evidence requires the
    flanking and all inserted read bases to meet *min_base_quality*;
    deletion evidence requires the flanking read bases to meet it.  Raw
    events are minimized and left-aligned with the reference-aware rules
    of :func:`~genome_variant.vcf.normalize_record` and equivalent events
    are merged by ``CHROM``, ``POS``, ``REF`` and ``ALT``.  An event's
    ``DP`` is the number of mapped reads spanning both its flanks with
    qualifying flank qualities and no other insertion or deletion at its
    boundaries; ``AC`` is the number of those reads supporting the event.
    The same *min_alt_count*, *min_alt_fraction* and *homozygous_fraction*
    thresholds apply, ``AD`` is ``DP-AC,AC`` and SNV and indel records
    are ordered by reference input order, ``POS``, ``REF`` and ``ALT``.

    *min_base_quality* is a non-boolean integer in 0-93 (default 20);
    *min_alt_count* a non-boolean positive integer (default 2);
    *min_alt_fraction* and *homozygous_fraction* finite numbers in 0-1
    (defaults 0.2 and 0.8), with the homozygous fraction no smaller than
    the minimum ALT fraction; *max_indel_length* a non-boolean positive
    integer (default 50).  Invalid thresholds raise :class:`ValueError`
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
    indel_length = _positive_threshold("max_indel_length", max_indel_length)

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

    # Per mapped read indel observations, only collected when enabled:
    # (reference index, aligned qualities, normalized events).
    sequences = {
        reference.identifier: reference.sequence
        for reference in reference_records
    }
    indel_observations: list[_IndelObservation] = []

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
        if call_indels:
            aligned_at, events = _observe_indels(
                read,
                quality,
                reference_records[ref_index],
                alignment,
                strand,
                base_quality,
                indel_length,
                sequences,
            )
            indel_observations.append((ref_index, aligned_at, events))

    records: list[tuple[int, VcfRecord]] = []
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
                records.append((ref_index, record))

    if call_indels:
        records.extend(
            _call_indels(
                reference_records,
                indel_observations,
                base_quality,
                alt_count,
                alt_fraction,
                homo_fraction,
            )
        )
    # SNV records alone are already produced in this order; the sort only
    # interleaves indel records and is stable for equal keys.
    records.sort(
        key=lambda item: (item[0], item[1].pos, item[1].ref, item[1].alt_text)
    )

    header = VcfHeader(_META_LINES, (sample,))
    return VcfFile(header, tuple(record for _, record in records))


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


def _event_gaps(pos: int, ref: str, alt: str) -> frozenset[int]:
    """The reference gap positions a normalized event touches.

    A gap numbered ``g`` is the boundary between the zero-based reference
    positions ``g - 1`` and ``g``.  An insertion touches only its
    insertion point; a deletion touches both edges of the deleted
    fragment.
    """
    if len(alt) > len(ref):
        return frozenset((pos,))
    return frozenset((pos, pos + len(ref) - 1))


def _observe_indels(
    read: SequenceRecord,
    quality: tuple[int, ...],
    reference: SequenceRecord,
    alignment,
    strand: str,
    min_base_quality: int,
    max_indel_length: int,
    sequences: dict[str, str],
) -> tuple[dict[int, int], _IndelEventMap]:
    """Extract one mapped read's indel observations from its alignment.

    Returns ``(aligned_at, events)``.  ``aligned_at`` maps a zero-based
    reference position to the Phred quality of the read base aligned
    there.  ``events`` maps each normalized event key ``(pos, ref, alt)``
    (1-based POS) to ``(gaps, quality_ok)``: the reference gap positions
    the event touches and whether this read's supporting evidence meets
    the quality floor.  A read contributes at most one entry per
    normalized event.
    """
    read_length = len(read.sequence)
    aligned_reference = alignment.aligned_reference
    aligned_query = alignment.aligned_query
    columns = len(aligned_reference)

    # Per column: the zero-based reference position (None on a reference
    # gap), the query index on the aligned strand (None on a query gap)
    # and the next unconsumed zero-based reference position.
    ref_positions: list[int | None] = []
    query_indices: list[int | None] = []
    next_reference: list[int] = []
    reference_position = alignment.reference_start
    query_index = alignment.query_start
    for column in range(columns):
        next_reference.append(reference_position)
        reference_char = aligned_reference[column]
        query_char = aligned_query[column]
        ref_positions.append(
            reference_position if reference_char != "-" else None
        )
        query_indices.append(query_index if query_char != "-" else None)
        if reference_char != "-":
            reference_position += 1
        if query_char != "-":
            query_index += 1

    def column_quality(column: int) -> int:
        index = query_indices[column]
        assert index is not None
        if strand == "+":
            return quality[index]
        # Map the quality back to the original read coordinates; the
        # quality string is never reverse complemented.
        return quality[read_length - 1 - index]

    aligned_at: dict[int, int] = {}
    for column in range(columns):
        ref_position = ref_positions[column]
        if ref_position is not None and query_indices[column] is not None:
            aligned_at[ref_position] = column_quality(column)

    events: _IndelEventMap = {}
    reference_sequence = reference.sequence
    column = 0
    while column < columns:
        insertion = aligned_reference[column] == "-"
        deletion = not insertion and aligned_query[column] == "-"
        if not insertion and not deletion:
            column += 1
            continue
        # One contiguous I or D CIGAR segment: a maximal run of columns
        # gapped on the same side.
        start = column
        while column < columns:
            reference_gap = aligned_reference[column] == "-"
            query_gap = aligned_query[column] == "-"
            if reference_gap != insertion or query_gap != deletion:
                break
            column += 1
        end = column

        length = end - start
        if length > max_indel_length:
            continue
        # The event needs an aligned read base on both sides.
        if start == 0 or end == columns:
            continue
        if query_indices[start - 1] is None or query_indices[end] is None:
            continue
        # The gap sits immediately left of the next reference position;
        # the base before it is the left anchor.
        gap = next_reference[start]
        anchor = gap - 1
        if anchor < 0:
            continue
        anchor_base = reference_sequence[anchor]
        if anchor_base not in _ACGT_SET:
            continue
        if insertion:
            inserted = aligned_query[start:end]
            if any(base not in _ACGT_SET for base in inserted):
                continue
            ref_allele = anchor_base
            alt_allele = anchor_base + inserted
        else:
            deleted = reference_sequence[gap : gap + length]
            if any(base not in _ACGT_SET for base in deleted):
                continue
            ref_allele = reference_sequence[anchor : gap + length]
            alt_allele = anchor_base

        quality_ok = (
            column_quality(start - 1) >= min_base_quality
            and column_quality(end) >= min_base_quality
        )
        if quality_ok and insertion:
            quality_ok = all(
                column_quality(inserted_column) >= min_base_quality
                for inserted_column in range(start, end)
            )

        # Minimize and left-align with the existing reference-aware
        # rules, then merge equivalent events.
        raw = VcfRecord(
            chrom=reference.identifier,
            pos=anchor + 1,
            id=".",
            ref=ref_allele,
            alt=(alt_allele,),
            qual=".",
            filter="PASS",
            info="",
        )
        normalized = normalize_record(raw, sequences)
        key = (normalized.pos, normalized.ref, normalized.alt[0])
        gaps = _event_gaps(normalized.pos, normalized.ref, normalized.alt[0])
        existing = events.get(key)
        if existing is None:
            events[key] = (gaps, quality_ok)
        else:
            events[key] = (existing[0], existing[1] or quality_ok)

    return aligned_at, events


def _call_indels(
    reference_records: list[SequenceRecord],
    observations: list[_IndelObservation],
    min_base_quality: int,
    min_alt_count: int,
    min_alt_fraction: float,
    homo_fraction: float,
) -> list[tuple[int, VcfRecord]]:
    """Build the VCF records for the observed indel candidates."""
    records: list[tuple[int, VcfRecord]] = []
    for ref_index, reference in enumerate(reference_records):
        reads_here = [
            (aligned_at, events)
            for observed_index, aligned_at, events in observations
            if observed_index == ref_index
        ]
        if not reads_here:
            continue
        candidates: dict[tuple[int, str, str], frozenset[int]] = {}
        for _, events in reads_here:
            for key, (gaps, _) in events.items():
                candidates.setdefault(key, gaps)
        for (pos, ref, alt), gaps in candidates.items():
            # Zero-based reference positions flanking the event.
            left = pos - 1
            right = pos - 1 + len(ref)
            depth = 0
            alt_support = 0
            for aligned_at, events in reads_here:
                left_quality = aligned_at.get(left)
                right_quality = aligned_at.get(right)
                if left_quality is None or right_quality is None:
                    continue
                if (
                    left_quality < min_base_quality
                    or right_quality < min_base_quality
                ):
                    continue
                key = (pos, ref, alt)
                if any(
                    other_key != key
                    and not other_gaps.isdisjoint(gaps)
                    for other_key, (other_gaps, _) in events.items()
                ):
                    continue
                depth += 1
                evidence = events.get(key)
                if evidence is not None and evidence[1]:
                    alt_support += 1
            if alt_support < min_alt_count:
                continue
            fraction = alt_support / depth
            if fraction < min_alt_fraction:
                continue
            genotype = "1/1" if fraction >= homo_fraction else "0/1"
            info = f"DP={depth};AC={alt_support};AF={fraction:.6f}"
            sample = f"{genotype}:{depth}:{depth - alt_support},{alt_support}"
            records.append(
                (
                    ref_index,
                    VcfRecord(
                        chrom=reference.identifier,
                        pos=pos,
                        id=".",
                        ref=ref,
                        alt=(alt,),
                        qual=".",
                        filter="PASS",
                        info=info,
                        format_text="GT:DP:AD",
                        sample_text=(sample,),
                    ),
                )
            )
    return records
