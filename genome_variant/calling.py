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

from ._columns import alignment_columns, original_read_position
from .alignment import _score_parameter
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
# One normalized indel candidate: (reference index, pos, ref, alt).
_IndelKey = tuple[int, int, str, str]

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
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
    min_score: int = 1,
) -> VcfFile:
    """Call single-sample SNVs and return them as a :class:`~genome_variant.vcf.VcfFile`.

    Each read is mapped exactly as in
    :func:`~genome_variant.mapping.map_reads` with the given
    *match_score*, *mismatch_penalty*, *gap_open*, *gap_extend* and
    *min_score* (same local-alignment scoring, both strands, same
    candidate adjudication); unmapped reads contribute nothing.  Every
    mapped read contributes at most one observation per covered reference
    position.
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
    integer (default 50).  The mapping parameters *match_score*,
    *mismatch_penalty*, *gap_open*, *gap_extend* and *min_score* have the
    same meaning, defaults and constraints as for
    :func:`~genome_variant.mapping.map_reads` and select the winning
    candidate of every read.  Invalid thresholds or mapping parameters
    raise :class:`ValueError` before either input is consumed.  Reads
    must carry quality values; a
    read without them or with a quality/sequence length mismatch raises
    :class:`~genome_variant.quality.ReadQualityError` when that read is
    reached.  An empty reference collection, duplicate reference
    identifiers or an empty/whitespace-only/tab-containing *sample_name*
    raise :class:`VariantCallingError`.

    *reads* may be any iterable, including one that can be traversed only
    once: each read is consumed and mapped exactly once, and resident
    memory stays bounded by the reference records, the per-position SNV
    counts and the distinct normalized indel candidates — never by the
    number of consumed reads.
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

    # One counter table per reference record: position -> [A, C, G, T].
    counts_per_reference: list[dict[int, list[int]]] = [
        {} for _ in reference_records
    ]

    # Indel evidence is aggregated read by read (only when enabled), so
    # resident memory stays bounded by the reference, the per-position
    # SNV counts and the distinct normalized candidates — never by the
    # number of consumed reads.
    sequences = {
        reference.identifier: reference.sequence
        for reference in reference_records
    }
    indels = (
        _IndelTracker(
            [len(reference.sequence) for reference in reference_records],
            base_quality,
            indel_length,
        )
        if call_indels
        else None
    )

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
        columns = alignment_columns(alignment)
        _collect_evidence(
            read,
            quality,
            reference_records[ref_index].sequence,
            columns,
            strand,
            base_quality,
            counts_per_reference[ref_index],
        )
        if call_indels:
            aligned_at, events = _observe_indels(
                read,
                quality,
                reference_records[ref_index],
                columns,
                strand,
                base_quality,
                indel_length,
                sequences,
            )
            assert indels is not None
            indels.observe(ref_index, aligned_at, events)

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
        assert indels is not None
        records.extend(
            _call_indels(
                reference_records,
                indels.candidates,
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
    columns,
    strand: str,
    min_base_quality: int,
    counts: dict[int, list[int]],
) -> None:
    """Accumulate admitted substitution observations for one read.

    Columns carry the single shared interpretation of the winning
    alignment; only ``=``/``X`` columns contribute an observation, each
    read at most once per reference position.
    """
    read_length = len(read.sequence)

    for column in columns:
        if column.op not in ("=", "X"):
            # Insertions have no reference position; deletions add no
            # base observation; neither contributes SNV evidence.
            continue

        # The observed base is on the aligned strand, already expressed
        # on the reference's forward strand for a reverse-strand win;
        # only the quality is mapped back to the original read.
        quality_position = original_read_position(column, read_length, strand)
        if quality[quality_position] >= min_base_quality:
            reference_position = column.reference_position
            assert reference_position is not None
            observed = column.query_base
            reference_base = reference_sequence[reference_position]
            assert observed is not None
            if observed in _ACGT_SET and reference_base in _ACGT_SET:
                table = counts.setdefault(reference_position, [0, 0, 0, 0])
                table[_ACGT.index(observed)] += 1


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
    columns,
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

    The columns use the one shared interpretation of the winning
    alignment: their ``boundary`` is the reference gap at the column's
    left edge (the consumed position for ``=``/``X``/``D`` and the
    insertion point for ``I``).
    """
    read_length = len(read.sequence)

    def column_quality(column) -> int:
        index = original_read_position(column, read_length, strand)
        return quality[index]

    aligned_at: dict[int, int] = {}
    for column in columns:
        if column.reference_position is not None and column.query_position is not None:
            aligned_at[column.reference_position] = column_quality(column)

    events: _IndelEventMap = {}
    reference_sequence = reference.sequence
    column_count = len(columns)
    index = 0
    while index < column_count:
        op = columns[index].op
        if op not in ("I", "D"):
            index += 1
            continue
        insertion = op == "I"
        # One contiguous I or D CIGAR segment: a maximal run of columns
        # gapped on the same side.
        start = index
        while index < column_count and columns[index].op == op:
            index += 1
        end = index

        length = end - start
        if length > max_indel_length:
            continue
        # The event needs an aligned read base on both sides.
        if start == 0 or end == column_count:
            continue
        if (
            columns[start - 1].query_position is None
            or columns[end].query_position is None
        ):
            continue
        # The gap sits immediately left of the next reference position;
        # the base before it is the left anchor.
        gap = columns[start].boundary
        anchor = gap - 1
        if anchor < 0:
            continue
        anchor_base = reference_sequence[anchor]
        if anchor_base not in _ACGT_SET:
            continue
        if insertion:
            inserted = "".join(
                column.query_base or "" for column in columns[start:end]
            )
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
            column_quality(columns[start - 1]) >= min_base_quality
            and column_quality(columns[end]) >= min_base_quality
        )
        if quality_ok and insertion:
            quality_ok = all(
                column_quality(inserted_column) >= min_base_quality
                for inserted_column in columns[start:end]
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


class _IndelTracker:
    """Streaming aggregation of one sample's indel evidence.

    Reads are observed one at a time; nothing per read is retained after
    :meth:`observe` returns.  Resident state is only:

    - ``candidates``: the distinct normalized events observed so far,
      mapped to their aggregate ``[gaps, depth, alt_support]``;
    - the per-reference aggregates used to backfill the depth of a
      candidate that is first observed in a later read.

    An event's ``DP`` counts the mapped reads on its reference that span
    both flanks (zero-based ``left = pos - 1`` and
    ``right = pos - 1 + len(ref)``) with flank qualities of at least
    *min_base_quality* and that carry no other event touching the
    candidate's gaps.  Because ``right - left = len(ref)`` never exceeds
    ``max_indel_length + 1``, the number of reads admitting each flank
    pair ``(left, left + distance)`` with ``distance`` up to that bound
    is sufficient to reconstruct the spanning count of any candidate
    registered later.  Reads carrying an event that touches a gap ``b``
    of a future candidate are excluded through three boundary
    statistics: counts of reads with an event touching ``b`` that admit
    the pair ``(b - 1, b - 1 + distance)`` (``b`` as the candidate's
    left gap) or ``(b - distance, b)`` (``b`` as the right gap), and, for
    the two-gap deletion case, counts of reads touching both gaps that
    admit the pair ``(b1 - 1, b2)`` — combined by inclusion-exclusion
    over the at most two gaps of a candidate.  All of these are bounded
    by the reference lengths and *max_indel_length*, never by the number
    of consumed reads.
    """

    def __init__(
        self,
        reference_lengths: list[int],
        min_base_quality: int,
        max_indel_length: int,
    ) -> None:
        self._min_base_quality = min_base_quality
        # (ref_index, pos, ref, alt) -> [gaps, depth, alt_support].
        self.candidates: dict[_IndelKey, list[object]] = {}
        # (ref_index, left flank) -> candidate keys registered there.
        self._by_left: dict[tuple[int, int], list[_IndelKey]] = {}
        # Per reference, left flank -> distance -> reads admitting both.
        self._flank_pairs: list[dict[int, dict[int, int]]] = [
            {} for _ in reference_lengths
        ]
        # Per reference, touched gap -> distance -> reads admitting the
        # pair (gap - 1, gap - 1 + distance) / (gap - distance, gap).
        self._gap_left: list[dict[int, dict[int, int]]] = [
            {} for _ in reference_lengths
        ]
        self._gap_right: list[dict[int, dict[int, int]]] = [
            {} for _ in reference_lengths
        ]
        # Per reference, (first gap, second gap) -> reads touching both
        # and admitting the pair (first - 1, second).
        self._gap_pairs: list[dict[tuple[int, int], int]] = [
            {} for _ in reference_lengths
        ]
        # Per reference, the largest flank distance a candidate can have:
        # len(ref) <= max_indel_length + 1 and both flanks must fit.
        self._max_distance = [
            min(max_indel_length + 1, max(length - 1, 0))
            for length in reference_lengths
        ]

    def observe(
        self,
        ref_index: int,
        aligned_at: dict[int, int],
        events: _IndelEventMap,
    ) -> None:
        """Fold one mapped read's indel observation into the aggregates.

        ``aligned_at`` maps zero-based reference positions to the quality
        of the read base aligned there; ``events`` is the read's
        normalized event map as produced by :func:`_observe_indels`.
        """
        admitted = {
            position
            for position, quality in aligned_at.items()
            if quality >= self._min_base_quality
        }
        max_distance = self._max_distance[ref_index]

        # Register newly observed events.  A candidate first observed by
        # this read cannot be carried by any earlier read, so its depth
        # from earlier reads is the spanning count minus the reads whose
        # own events touch the candidate's gaps.
        for (pos, ref, alt), (gaps, _quality_ok) in events.items():
            key = (ref_index, pos, ref, alt)
            if key in self.candidates:
                continue
            left = pos - 1
            distance = len(ref)
            depth = self._flank_pairs[ref_index].get(left, {}).get(distance, 0)
            touched = sorted(gaps)
            blocked = 0
            for boundary in touched:
                if boundary == left + 1:
                    blocked += (
                        self._gap_left[ref_index]
                        .get(boundary, {})
                        .get(distance, 0)
                    )
                else:  # boundary == left + distance (the right flank gap)
                    blocked += (
                        self._gap_right[ref_index]
                        .get(boundary, {})
                        .get(distance, 0)
                    )
            if len(touched) == 2:
                blocked -= self._gap_pairs[ref_index].get(
                    (touched[0], touched[1]), 0
                )
            self.candidates[key] = [gaps, depth - blocked, 0]
            self._by_left.setdefault((ref_index, left), []).append(key)

        # This read's own contribution to every registered candidate it
        # spans: admitted flanks and no other event of this read touching
        # the candidate's gaps.
        for left in admitted:
            for key in self._by_left.get((ref_index, left), ()):
                _, pos, ref, alt = key
                if pos - 1 + len(ref) not in admitted:
                    continue
                gaps = self.candidates[key][0]
                if any(
                    other != (pos, ref, alt)
                    and not other_gaps.isdisjoint(gaps)
                    for other, (other_gaps, _) in events.items()
                ):
                    continue
                self.candidates[key][1] += 1
                evidence = events.get((pos, ref, alt))
                if evidence is not None and evidence[1]:
                    self.candidates[key][2] += 1

        # Fold the read into the aggregates consulted by later
        # registrations; only then is this read part of the backfill.
        if not max_distance:
            return
        flank_pairs = self._flank_pairs[ref_index]
        for left in admitted:
            row = flank_pairs.setdefault(left, {})
            for distance in range(1, max_distance + 1):
                if left + distance in admitted:
                    row[distance] = row.get(distance, 0) + 1
        touched = sorted(
            {gap for gaps, _ in events.values() for gap in gaps}
        )
        if not touched:
            return
        gap_left = self._gap_left[ref_index]
        gap_right = self._gap_right[ref_index]
        for boundary in touched:
            if boundary - 1 in admitted:
                row = gap_left.setdefault(boundary, {})
                for distance in range(1, max_distance + 1):
                    if boundary - 1 + distance in admitted:
                        row[distance] = row.get(distance, 0) + 1
            if boundary in admitted:
                row = gap_right.setdefault(boundary, {})
                for distance in range(1, max_distance + 1):
                    if boundary - distance in admitted:
                        row[distance] = row.get(distance, 0) + 1
        gap_pairs = self._gap_pairs[ref_index]
        for index, first in enumerate(touched):
            if first - 1 not in admitted:
                continue
            for second in touched[index + 1 :]:
                if second - first >= max_distance:
                    break
                if second in admitted:
                    pair = (first, second)
                    gap_pairs[pair] = gap_pairs.get(pair, 0) + 1


def _call_indels(
    reference_records: list[SequenceRecord],
    candidates: dict[_IndelKey, list[object]],
    min_alt_count: int,
    min_alt_fraction: float,
    homo_fraction: float,
) -> list[tuple[int, VcfRecord]]:
    """Build the VCF records for the aggregated indel candidates."""
    records: list[tuple[int, VcfRecord]] = []
    for (ref_index, pos, ref, alt), (
        _gaps,
        depth,
        alt_support,
    ) in candidates.items():
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
                    chrom=reference_records[ref_index].identifier,
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
