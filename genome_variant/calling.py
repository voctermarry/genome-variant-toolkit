"""Single-sample SNV and optional short-indel calling from mapped reads.

Public API:

- :class:`VariantCallingError` — the reference or sample configuration is
  unusable for calling (empty reference, duplicate reference identifier
  or invalid sample name).
- :func:`call_variants` — map FASTQ reads with the same local-alignment
  scoring, strand search and candidate adjudication as
  :func:`~genome_variant.mapping.map_reads`, pile up the quality-filtered
  substitution evidence (and, when enabled, short insertion/deletion
  evidence) and return a VCFv4.2 document.

Only substitutions between the bases ``A``, ``C``, ``G`` and ``T`` are
called by default.  Short insertions and deletions are additionally
called only when *call_indels* is enabled; they are taken from the
contiguous ``I`` and ``D`` CIGAR runs of the winning alignment only.
Ambiguous bases never contribute evidence.  Quality values never
participate in mapping; they only filter base evidence.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

from .mapping import _best_candidate
from .quality import ReadQualityError
from .sequence_io import SequenceRecord
from .vcf import VcfFile, VcfHeader, VcfRecord, _minimize

__all__ = [
    "VariantCallingError",
    "call_variants",
]

_ACGT = ("A", "C", "G", "T")
_ACGT_SET = frozenset(_ACGT)
_PHRED_MIN = 0
_PHRED_MAX = 93
_DEFAULT_MAX_INDEL_LENGTH = 50

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
    max_indel_length: int = _DEFAULT_MAX_INDEL_LENGTH,
) -> VcfFile:
    """Call single-sample variants and return them as a :class:`~genome_variant.vcf.VcfFile`.

    Each read is mapped exactly as in :func:`~genome_variant.mapping.map_reads`
    (same local-alignment scoring, both strands, same candidate
    adjudication); unmapped reads contribute nothing.  Every mapped read
    contributes at most one observation per covered reference position.
    A substitution base is admitted as evidence only when its Phred
    quality is at least *min_base_quality* and both the reference base and
    the observed base are among ``A``, ``C``, ``G`` and ``T``; on the
    reverse strand the quality and base are taken from the corresponding
    original-read position.  Insertion, deletion and unaligned alignment
    columns never contribute substitution evidence.

    When *call_indels* is false (the default) only substitutions are
    called and the other arguments and output are unchanged.  When it is
    true, contiguous ``I``/``D`` CIGAR runs of the winning alignment also
    yield insertion/deletion candidates, at most one observation per read
    per normalized event.  A run needs an aligned read base on each side;
    every involved reference or inserted base must be among
    ``A``, ``C``, ``G`` and ``T``; its length must not exceed
    *max_indel_length*; and the two flanking read bases, plus every
    inserted base for an insertion, must have Phred quality at least
    *min_base_quality* (reverse-strand qualities are mapped back to the
    original read).  Events are reference-minimized and left-aligned and
    equivalent events are merged.  An event's ``DP`` counts mapped reads
    spanning both of its boundaries with both flanks quality-qualified
    and producing no other insertion or deletion at either boundary;
    ``AC`` is the subset supporting that allele.  A record is emitted
    only when ``AC`` is at least *min_alt_count* and ``AC/DP`` is at
    least *min_alt_fraction*; ``GT`` is ``1/1`` when that share is at
    least *homozygous_fraction*, otherwise ``0/1``.  SNV and indel
    records are all retained, sorted by reference input order, POS, REF
    and ALT; distinct qualifying indels at the same position are written
    on separate rows.

    *min_base_quality* is a non-boolean integer in 0-93 (default 20);
    *min_alt_count* a non-boolean positive integer (default 2);
    *min_alt_fraction* and *homozygous_fraction* finite numbers in 0-1
    (defaults 0.2 and 0.8), with the homozygous fraction no smaller than
    the minimum ALT fraction; *call_indels* a boolean (default false);
    *max_indel_length* a non-boolean positive integer (default 50),
    validated before either input is consumed.  Invalid thresholds raise
    :class:`ValueError` before either input is consumed.  Reads must
    carry quality values; a read without them or with a quality/sequence
    length mismatch raises :class:`~genome_variant.quality.ReadQualityError`
    when that read is reached.  An empty reference collection, duplicate
    reference identifiers or an empty/whitespace-only/tab-containing
    *sample_name* raise :class:`VariantCallingError`.
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
    if not isinstance(call_indels, bool):
        kind = type(call_indels).__name__
        raise ValueError(f"call_indels must be a boolean, not {kind}")
    indel_limit = _positive_threshold("max_indel_length", max_indel_length)

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
    # Per reference record, one descriptor per mapped read; events are
    # only known once every read has been seen, so DP/AC are aggregated
    # afterwards.
    indel_reads_per_reference: list[list["_IndelReadEvidence"]] = [
        [] for _ in reference_records
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
        if call_indels:
            indel_reads_per_reference[ref_index].append(
                _collect_indel_evidence(
                    read,
                    quality,
                    reference_records[ref_index].sequence,
                    alignment,
                    strand,
                    base_quality,
                    indel_limit,
                )
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
        if call_indels:
            for event_key, supporting, depth in _aggregate_indels(
                indel_reads_per_reference[ref_index]
            ):
                position, ref, alt = event_key
                record = _call_indel(
                    chromosome,
                    position,
                    ref,
                    alt,
                    supporting,
                    depth,
                    alt_count,
                    alt_fraction,
                    homo_fraction,
                )
                if record is not None:
                    records.append(record)

    if call_indels:
        # Reference input order first (records were appended per
        # reference), then POS, REF and ALT.
        order = {
            reference.identifier: index
            for index, reference in enumerate(reference_records)
        }
        records.sort(key=lambda record: (order[record.chrom], record.pos,
                                         record.ref, record.alt))

    header = VcfHeader(_META_LINES, (sample,))
    return VcfFile(header, tuple(records))


@dataclass(frozen=True)
class _IndelReadEvidence:
    """Indel-relevant evidence of one mapped read.

    ``covered`` holds the zero-based reference positions carrying an
    aligned (match/mismatch) read base whose quality meets the threshold.
    ``events`` lists the distinct normalized
    ``(1-based pos, ref, alt)`` insertion/deletion events the read
    supports; a read contributes at most one observation to an event.
    """

    covered: frozenset[int]
    events: frozenset[tuple[int, str, str]] = frozenset()


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


def _read_column_quality(
    quality: tuple[int, ...],
    strand: str,
    read_length: int,
    query_start: int,
    query_column: int,
) -> int:
    """Phred quality of an aligned (reference-oriented) query column."""
    if strand == "+":
        quality_position = query_start + query_column
    else:
        quality_position = read_length - 1 - (query_start + query_column)
    return quality[quality_position]


def _collect_indel_evidence(
    read: SequenceRecord,
    read_quality: tuple[int, ...],
    reference_sequence: str,
    alignment,
    strand: str,
    min_base_quality: int,
    max_indel_length: int,
) -> _IndelReadEvidence:
    """Extract indel-relevant evidence for one mapped read.

    Every mapped read yields a descriptor: ``covered`` records the
    reference positions carrying a quality-qualified aligned read base
    (used to decide spanning depth), ``events`` lists the distinct
    normalized insertion/deletion events the read supports.
    """
    aligned_reference = alignment.aligned_reference
    aligned_query = alignment.aligned_query
    read_length = len(read.sequence)
    query_start = alignment.query_start

    # One entry per aligned column: (reference_position_or_None,
    # query_column_or_None, consumes_reference, consumes_query).
    columns: list[tuple[int | None, int | None, bool, bool]] = []
    reference_position = alignment.reference_start
    query_column = 0
    for column in range(len(aligned_reference)):
        consumes_reference = aligned_reference[column] != "-"
        consumes_query = aligned_query[column] != "-"
        columns.append(
            (
                reference_position if consumes_reference else None,
                query_column if consumes_query else None,
                consumes_reference,
                consumes_query,
            )
        )
        if consumes_reference:
            reference_position += 1
        if consumes_query:
            query_column += 1

    covered: set[int] = set()
    for ref_pos, qry_col, consumes_reference, consumes_query in columns:
        if consumes_reference and consumes_query:
            assert ref_pos is not None and qry_col is not None
            flank_quality = _read_column_quality(
                read_quality, strand, read_length, query_start, qry_col
            )
            if flank_quality >= min_base_quality:
                covered.add(ref_pos)

    events: set[tuple[int, str, str]] = set()
    column_count = len(columns)
    column = 0
    while column < column_count:
        consumes_reference = columns[column][2]
        consumes_query = columns[column][3]
        if consumes_reference == consumes_query:
            column += 1
            continue

        gap_op = "D" if consumes_reference else "I"
        end = column + 1
        while end < column_count and columns[end][2] == consumes_reference and (
            columns[end][3] == consumes_query
        ):
            end += 1

        left = column - 1
        right = end
        if (
            left >= 0
            and right < column_count
            and columns[left][2]
            and columns[left][3]
            and columns[right][2]
            and columns[right][3]
            and end - column <= max_indel_length
        ):
            raw_event = _build_indel_event(
                gap_op,
                column,
                end,
                left,
                right,
                columns,
                aligned_reference,
                aligned_query,
                reference_sequence,
                read_quality,
                strand,
                read_length,
                query_start,
                min_base_quality,
            )
            if raw_event is not None:
                normalized = _normalize_event(*raw_event, reference_sequence)
                if normalized is not None:
                    events.add(normalized)
        column = end

    return _IndelReadEvidence(
        covered=frozenset(covered), events=frozenset(events)
    )


def _aggregate_indels(
    descriptors: list[_IndelReadEvidence],
) -> list[tuple[tuple[int, str, str], int, int]]:
    """Merge equivalent events and compute ``(event, AC, DP)`` per event.

    An event's ``DP`` counts reads whose quality-qualified aligned bases
    cover both of its boundary reference positions and that produce no
    other insertion or deletion at either boundary.  ``AC`` is the subset
    that supports the event itself.  A read supports a normalized event
    at most once.
    """
    all_events: set[tuple[int, str, str]] = set()
    for descriptor in descriptors:
        all_events.update(descriptor.events)

    result: list[tuple[tuple[int, str, str], int, int]] = []
    for event in all_events:
        position, ref, alt = event
        left_boundary, right_boundary = _event_boundaries(position, ref)

        depth = 0
        supporting = 0
        for descriptor in descriptors:
            if event in descriptor.events:
                blocked = any(
                    other != event
                    and _touches_boundary(
                        other, left_boundary, right_boundary
                    )
                    for other in descriptor.events
                )
                if not blocked:
                    depth += 1
                    supporting += 1
            elif (
                left_boundary in descriptor.covered
                and right_boundary in descriptor.covered
                and not any(
                    _touches_boundary(other, left_boundary, right_boundary)
                    for other in descriptor.events
                )
            ):
                depth += 1

        result.append((event, supporting, depth))

    result.sort(key=lambda entry: entry[0])
    return result


def _event_boundaries(position: int, ref: str) -> tuple[int, int]:
    """Zero-based reference positions of an event's two flanking bases.

    The left boundary is the anchor base and the right boundary is the
    first reference base past the event: one position past the anchor for
    an insertion (its REF is just the anchor) and one position past the
    deleted segment for a deletion.
    """
    return position - 1, position - 1 + len(ref)


def _touches_boundary(
    event: tuple[int, str, str], left_boundary: int, right_boundary: int
) -> bool:
    """Whether *event*'s own flank positions touch either boundary."""
    position, ref, _alt = event
    event_left, event_right = _event_boundaries(position, ref)
    return (
        event_left == left_boundary
        or event_left == right_boundary
        or event_right == left_boundary
        or event_right == right_boundary
    )

def _build_indel_event(
    gap_op: str,
    gap_start: int,
    gap_end: int,
    left: int,
    right: int,
    columns,
    aligned_reference: str,
    aligned_query: str,
    reference_sequence: str,
    read_quality: tuple[int, ...],
    strand: str,
    read_length: int,
    query_start: int,
    min_base_quality: int,
) -> tuple[int, str, str] | None:
    """Build the raw (pos, ref, alt) anchored on the left aligned column."""
    anchor_reference_position = columns[left][0]
    assert anchor_reference_position is not None
    anchor_base = reference_sequence[anchor_reference_position]
    left_query_column = columns[left][1]
    right_query_column = columns[right][1]
    assert left_query_column is not None and right_query_column is not None
    left_quality = _read_column_quality(
        read_quality, strand, read_length, query_start, left_query_column
    )
    right_quality = _read_column_quality(
        read_quality, strand, read_length, query_start, right_query_column
    )
    if left_quality < min_base_quality or right_quality < min_base_quality:
        return None
    if anchor_base not in _ACGT_SET:
        return None

    if gap_op == "I":
        inserted = aligned_query[gap_start:gap_end]
        if any(base not in _ACGT_SET for base in inserted):
            return None
        for gap_column in range(gap_start, gap_end):
            insert_query_column = columns[gap_column][1]
            assert insert_query_column is not None
            insert_quality = _read_column_quality(
                read_quality,
                strand,
                read_length,
                query_start,
                insert_query_column,
            )
            if insert_quality < min_base_quality:
                return None
        return anchor_reference_position, anchor_base, anchor_base + inserted

    deleted = aligned_reference[gap_start:gap_end]
    if any(base not in _ACGT_SET for base in deleted):
        return None
    return (
        anchor_reference_position,
        anchor_base + deleted,
        anchor_base,
    )


def _normalize_event(
    position: int,
    ref: str,
    alt: str,
    reference_sequence: str,
) -> tuple[int, str, str] | None:
    """Reference-aware minimization and left alignment of one indel.

    Returns the 1-based ``(pos, ref, alt)`` or ``None`` when the event
    cannot be represented within the reference.
    """
    pos, minimized_ref, minimized_alts = _minimize(
        position + 1, ref, (alt,)
    )
    minimized_alt = minimized_alts[0]

    while pos > 1 and len(minimized_ref) != len(minimized_alt):
        preceding = reference_sequence[pos - 2].upper()
        if preceding not in _ACGT_SET:
            return None
        next_pos, next_ref, next_alts = _minimize(
            pos - 1,
            preceding + minimized_ref,
            (preceding + minimized_alt,),
        )
        next_alt = next_alts[0]
        if (next_pos, next_ref, next_alt) == (pos, minimized_ref, minimized_alt):
            break
        pos, minimized_ref, minimized_alt = next_pos, next_ref, next_alt

    if pos < 1 or pos - 1 + len(minimized_ref) > len(reference_sequence):
        return None
    expected = reference_sequence[pos - 1 : pos - 1 + len(minimized_ref)]
    if expected != minimized_ref or any(
        base not in _ACGT_SET for base in minimized_ref + minimized_alt
    ):
        return None
    return pos, minimized_ref, minimized_alt


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


def _call_indel(
    chromosome: str,
    pos: int,
    ref: str,
    alt: str,
    supporting: int,
    depth: int,
    min_alt_count: int,
    min_alt_fraction: float,
    homo_fraction: float,
) -> VcfRecord | None:
    """Build the VCF record for one normalized indel event, or ``None``."""
    if supporting < min_alt_count or depth == 0:
        return None
    fraction = supporting / depth
    if fraction < min_alt_fraction:
        return None

    genotype = "1/1" if fraction >= homo_fraction else "0/1"
    info = f"DP={depth};AC={supporting};AF={fraction:.6f}"
    sample = f"{genotype}:{depth}:{depth - supporting},{supporting}"
    return VcfRecord(
        chrom=chromosome,
        pos=pos,
        id=".",
        ref=ref,
        alt=(alt,),
        qual=".",
        filter="PASS",
        info=info,
        format_text="GT:DP:AD",
        sample_text=(sample,),
    )
