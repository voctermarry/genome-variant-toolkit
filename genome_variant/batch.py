"""Multi-sample variant calling from a JSON Lines reads manifest.

Public API:

- :class:`BatchManifestError` — the batch manifest is malformed (bad
  JSON, a wrong field set or field type, an invalid or duplicate sample
  name, or a ``reads`` path of ``"-"``).
- :class:`BatchManifestEntry` — one parsed manifest row: its 1-based
  line number, unique ``sample`` and resolved ``reads`` path.
- :func:`read_batch_manifest` — parse a JSON Lines manifest into
  :class:`BatchManifestEntry` objects in file order.
- :func:`batch_call_variants` — call every manifest sample's FASTQ reads
  with :func:`~genome_variant.calling.call_variants` semantics and merge
  the per-sample results with the reference-aware normalization rules of
  :func:`~genome_variant.summary.summarize_variants`, without
  intermediate VCF files.

The manifest is JSON Lines: each non-empty line is a JSON object with
exactly the string fields ``sample`` (non-empty, not whitespace-only,
tab-free, globally unique) and ``reads`` (a FASTQ path; ``"-"`` is not
allowed).  Relative ``reads`` paths in a file manifest resolve against
the manifest's own directory; relative paths in a stream manifest
resolve against the current working directory.
"""

from __future__ import annotations

import io
import json
import os
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from .alignment import _score_parameter
from .calling import VariantCallingError, call_variants
from .sequence_io import SequenceRecord, read_sequences
from .summary import SampleCall, VariantSummary
from .vcf import normalize_vcf

__all__ = [
    "BatchManifestError",
    "BatchManifestEntry",
    "read_batch_manifest",
    "batch_call_variants",
]


class BatchManifestError(ValueError):
    """The batch manifest is malformed or inconsistent."""


@dataclass(frozen=True)
class BatchManifestEntry:
    """One batch manifest row: 1-based line, unique ``sample``, ``reads``."""

    line_number: int
    sample: str
    reads: str


def read_batch_manifest(
    source: "str | os.PathLike[str] | io.TextIOBase",
) -> tuple[BatchManifestEntry, ...]:
    """Parse a whole JSON Lines batch manifest from *source*.

    *source* is either a filesystem path (``str`` or :class:`os.PathLike`)
    or an already-open text stream.  Blank lines are ignored; every other
    line must decode to a JSON object containing exactly the string fields
    ``sample`` (non-empty, not whitespace-only, tab-free, unique) and
    ``reads`` (a non-empty path; ``"-"`` is rejected).

    Relative ``reads`` values resolve against the manifest file's
    directory for a path source, or against the current working directory
    for a stream source.

    Raises :class:`BatchManifestError` carrying the 1-based manifest line
    number on any structural or field violation, and :class:`OSError`
    when opening a path fails.
    """
    if isinstance(source, io.IOBase):
        stream: io.TextIOBase = source
        source_name = getattr(stream, "name", None)
        if not isinstance(source_name, str) or not source_name:
            source_name = "<stream>"
        base_dir = os.getcwd()
        close = False
    elif isinstance(source, (str, os.PathLike)):
        path = os.fspath(source)
        stream = open(path, "r", encoding="utf-8", newline="")
        source_name = path
        base_dir = os.path.dirname(os.path.abspath(path))
        close = True
    else:
        raise TypeError("source must be a text path or a text stream")

    try:
        seen: set[str] = set()
        entries: list[BatchManifestEntry] = []

        for number, raw_line in enumerate(stream, start=1):
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise BatchManifestError(
                    f"{source_name}:{number}: invalid JSON: {exc.msg}"
                ) from None
            if not isinstance(value, dict):
                raise BatchManifestError(
                    f"{source_name}:{number}: manifest line must be a JSON object"
                )
            if set(value) != {"sample", "reads"}:
                raise BatchManifestError(
                    f"{source_name}:{number}: manifest line must contain exactly "
                    "the fields sample and reads"
                )
            sample = value["sample"]
            reads = value["reads"]
            if (
                not isinstance(sample, str)
                or not sample
                or sample.isspace()
                or "\t" in sample
            ):
                raise BatchManifestError(
                    f"{source_name}:{number}: sample must be a non-empty string "
                    "that is not whitespace-only and contains no tab characters"
                )
            if not isinstance(reads, str) or not reads:
                raise BatchManifestError(
                    f"{source_name}:{number}: reads must be a non-empty string path"
                )
            if reads == "-":
                raise BatchManifestError(
                    f"{source_name}:{number}: reads must not be '-' "
                    "(standard input)"
                )
            if sample in seen:
                raise BatchManifestError(
                    f"{source_name}:{number}: duplicate sample {sample!r}"
                )
            seen.add(sample)
            resolved = reads if os.path.isabs(reads) else os.path.join(base_dir, reads)
            entries.append(BatchManifestEntry(number, sample, resolved))

        return tuple(entries)
    finally:
        if close:
            stream.close()


def _sample_call(record, sample: str) -> SampleCall:
    """Extract the ``GT``/``DP``/``AD`` call of one caller-produced record."""
    format_fields = record.format_text.split(":")
    values = record.sample_text[0].split(":")
    data = dict(zip(format_fields, values))
    ad_ref, ad_alt = data["AD"].split(",")
    return SampleCall(sample, data["GT"], int(data["DP"]), (int(ad_ref), int(ad_alt)))


def batch_call_variants(
    manifest: "str | os.PathLike[str] | io.TextIOBase | Sequence[BatchManifestEntry]",
    reference: "str | os.PathLike[str] | io.TextIOBase | Iterable[SequenceRecord]",
    *,
    min_base_quality: int = 20,
    min_alt_count: int = 2,
    min_alt_fraction: float = 0.2,
    homozygous_fraction: float = 0.8,
    call_indels: bool = False,
    max_indel_length: int = 50,
    match_score: int = 2,
    mismatch_penalty: int = 3,
    gap_open: int = 5,
    gap_extend: int = 2,
    min_score: int = 1,
) -> tuple[VariantSummary, ...]:
    """Call variants for every manifest sample and merge the results.

    *manifest* is a batch manifest path/stream or an already-parsed
    sequence of :class:`BatchManifestEntry`; *reference* is a FASTA
    path/stream or an iterable of
    :class:`~genome_variant.sequence_io.SequenceRecord`.  Relative
    ``reads`` paths are resolved as described for
    :func:`read_batch_manifest`.

    Each entry's FASTQ reads are called independently with
    :func:`~genome_variant.calling.call_variants` (same quality
    thresholds, allele counting, allele fractions, genotyping and
    optional short-indel semantics), using the entry's ``sample`` as the
    sample name.  The per-sample results are normalized against the
    reference and merged by ``CHROM``, ``POS``, ``REF`` and ``ALT``
    exactly as in :func:`~genome_variant.summary.summarize_variants`;
    no intermediate VCF is produced.

    Returns :class:`~genome_variant.summary.VariantSummary` objects
    ordered by reference record order, then POS, REF and ALT; each
    summary's ``samples`` follow manifest order.  A run without variants
    returns an empty tuple.

    The threshold arguments and their constraints are those of
    :func:`~genome_variant.calling.call_variants`, including the mapping
    arguments *match_score*, *mismatch_penalty*, *gap_open*, *gap_extend*
    and *min_score*; the same mapping settings apply to every manifest
    sample.  Invalid thresholds raise :class:`ValueError` before any
    sample's reads are consumed.
    An empty reference collection or duplicate reference identifiers
    raise :class:`VariantCallingError`.  Manifest problems raise
    :class:`BatchManifestError`; malformed reads raise
    :class:`~genome_variant.sequence_io.SequenceFormatError` /
    :class:`~genome_variant.sequence_io.SequenceValidationError` and
    missing or invalid quality values raise
    :class:`~genome_variant.quality.ReadQualityError`, as in single-sample
    calling.  File access failures propagate as :class:`OSError`.
    """
    if isinstance(manifest, (str, os.PathLike, io.IOBase)):
        entries = read_batch_manifest(manifest)
    else:
        entries = tuple(manifest)

    # Validate the mapping parameters before the reference or any
    # sample's reads are consumed; call_variants re-validates them (and
    # the remaining thresholds) on the zero-read check below.
    _score_parameter("match_score", match_score, positive=True)
    _score_parameter("mismatch_penalty", mismatch_penalty, positive=False)
    _score_parameter("gap_open", gap_open, positive=False)
    _score_parameter("gap_extend", gap_extend, positive=False)
    if not isinstance(min_score, int) or isinstance(min_score, bool):
        kind = type(min_score).__name__
        raise ValueError(
            f"min_score must be a non-boolean integer, not {kind}"
        )
    if min_score < 1:
        raise ValueError("min_score must be a positive integer")

    if isinstance(reference, (str, os.PathLike, io.IOBase)):
        reference_records = tuple(read_sequences(reference, format="fasta"))
    else:
        reference_records = tuple(reference)

    # Validate the thresholds and the reference before any sample's reads
    # are consumed: a zero-read call performs exactly the call_variants
    # threshold and reference checks (and nothing else).
    call_variants(
        reference_records,
        (),
        min_base_quality=min_base_quality,
        min_alt_count=min_alt_count,
        min_alt_fraction=min_alt_fraction,
        homozygous_fraction=homozygous_fraction,
        call_indels=call_indels,
        max_indel_length=max_indel_length,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_open=gap_open,
        gap_extend=gap_extend,
        min_score=min_score,
    )

    # Reference order for sorting; identifiers are unique at this point.
    order = {
        record.identifier: index
        for index, record in enumerate(reference_records)
    }

    # key (chrom, pos, ref, alt) -> list of SampleCall in manifest order
    calls: dict[tuple[str, int, str, str], list[SampleCall]] = {}
    key_order: list[tuple[str, int, str, str]] = []

    for entry in entries:
        reads = read_sequences(entry.reads, format="fastq")
        document = call_variants(
            reference_records,
            reads,
            min_base_quality=min_base_quality,
            min_alt_count=min_alt_count,
            min_alt_fraction=min_alt_fraction,
            homozygous_fraction=homozygous_fraction,
            sample_name=entry.sample,
            call_indels=call_indels,
            max_indel_length=max_indel_length,
            match_score=match_score,
            mismatch_penalty=mismatch_penalty,
            gap_open=gap_open,
            gap_extend=gap_extend,
            min_score=min_score,
        )
        normalized = normalize_vcf(document, reference_records)
        for record in normalized.records:
            key = (record.chrom, record.pos, record.ref, record.alt[0])
            if key not in calls:
                calls[key] = []
                key_order.append(key)
            calls[key].append(_sample_call(record, entry.sample))

    summaries: list[VariantSummary] = []
    for key in key_order:
        chrom, pos, ref, alt = key
        summaries.append(
            VariantSummary(
                ref_index=order[chrom],
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                samples=tuple(calls[key]),
            )
        )

    summaries.sort(key=lambda summary: summary.ref_order_key)
    return tuple(summaries)
