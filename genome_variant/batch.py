"""Batch multi-sample variant calling directly from FASTQ reads.

Public API:

- :class:`BatchManifestError` — malformed manifest text, bad fields, an
  invalid or duplicate sample name, or a ``reads`` value of ``"-"``.
- :class:`BatchManifestEntry` — one parsed manifest row.
- :func:`read_batch_manifest` — parse a JSON Lines manifest into
  :class:`BatchManifestEntry` objects in file order.
- :func:`batch_call_variants` — call each manifest sample's reads exactly
  as :func:`~genome_variant.calling.call_variants` and merge the results
  with the reference-aware normalization rules of
  :func:`~genome_variant.summary.summarize_variants`, without writing
  intermediate VCFs.

The manifest is JSON Lines: each non-empty line is a JSON object with
exactly the string fields ``sample`` (non-empty, not whitespace-only,
tab-free, globally unique) and ``reads`` (a FASTQ path, never ``"-"``).
Relative ``reads`` paths in a file manifest resolve against the
manifest's own directory; relative paths in a stream manifest resolve
against the current working directory.
"""

from __future__ import annotations

import io
import json
import os
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

from .calling import (
    VariantCallingError,
    _fraction_threshold,
    _positive_threshold,
    _quality_threshold,
    call_variants,
)
from .sequence_io import SequenceRecord, read_sequences
from .summary import VariantSummary, _merge_sample_documents
from .vcf import VcfFile

__all__ = [
    "BatchManifestError",
    "BatchManifestEntry",
    "read_batch_manifest",
    "batch_call_variants",
]


class BatchManifestError(ValueError):
    """The batch manifest is invalid."""


@dataclass(frozen=True)
class BatchManifestEntry:
    """One manifest row: its 1-based line, a unique ``sample`` and ``reads``."""

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
    ``reads`` (a FASTQ path, never ``"-"``).

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
                    f"{source_name}:{number}: sample must be a non-empty, not "
                    "whitespace-only string without tab characters"
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
) -> tuple[VariantSummary, ...]:
    """Call every manifest sample from its FASTQ reads and merge the results.

    *manifest* is a manifest path/stream or an already-parsed sequence of
    :class:`BatchManifestEntry`; *reference* is a FASTA path/stream or an
    iterable of :class:`~genome_variant.sequence_io.SequenceRecord`.
    Relative ``reads`` paths are resolved as described for
    :func:`read_batch_manifest`.

    Each entry's reads are read in manifest order and called independently
    with the entry's ``sample`` as the sample name, exactly as
    :func:`~genome_variant.calling.call_variants` with the same threshold
    arguments (including *call_indels* and *max_indel_length*).  The
    per-sample results are then normalized against the reference and
    merged with the rules of
    :func:`~genome_variant.summary.summarize_variants`; no intermediate
    VCF is produced.

    Returns summaries ordered by reference record order, then POS, REF and
    ALT; each summary's ``samples`` follow manifest order.  An empty
    manifest or a batch without calls yields an empty tuple.

    The thresholds follow the same rules as
    :func:`~genome_variant.calling.call_variants` and invalid values raise
    :class:`ValueError` before any input is consumed.  Manifest problems
    raise :class:`BatchManifestError` (message carries the manifest source
    and 1-based line number).  An empty reference collection or duplicate
    reference identifiers raise
    :class:`~genome_variant.calling.VariantCallingError`; FASTQ format,
    base and quality problems raise the corresponding existing exceptions
    (:class:`~genome_variant.sequence_io.SequenceFormatError`,
    :class:`~genome_variant.sequence_io.SequenceValidationError`,
    :class:`~genome_variant.quality.ReadQualityError`); file access
    failures raise :class:`OSError`.
    """
    # Same threshold rules as call_variants, before any input is consumed.
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

    if isinstance(manifest, (str, os.PathLike, io.IOBase)):
        entries = read_batch_manifest(manifest)
    else:
        entries = tuple(manifest)

    if isinstance(reference, (str, os.PathLike, io.IOBase)):
        reference_records = tuple(read_sequences(reference, format="fasta"))
    else:
        reference_records = tuple(reference)

    # Same reference validation as call_variants, before any reads file is
    # opened.
    if not reference_records:
        raise VariantCallingError("at least one reference record is required")
    seen: set[str] = set()
    for record in reference_records:
        if record.identifier in seen:
            raise VariantCallingError(
                f"duplicate reference identifier {record.identifier!r}"
            )
        seen.add(record.identifier)

    def _documents() -> Iterator[tuple[str, str, VcfFile]]:
        for entry in entries:
            with open(entry.reads, "r", encoding="utf-8", newline="") as stream:
                read_records = list(read_sequences(stream, format="fastq"))
            document = call_variants(
                reference_records,
                read_records,
                min_base_quality=base_quality,
                min_alt_count=alt_count,
                min_alt_fraction=alt_fraction,
                homozygous_fraction=homo_fraction,
                sample_name=entry.sample,
                call_indels=call_indels,
                max_indel_length=indel_length,
            )
            yield entry.sample, entry.reads, document

    return _merge_sample_documents(_documents(), reference_records)
