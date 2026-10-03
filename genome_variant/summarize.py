"""Summarize single-sample VCFs into multi-sample per-variant summaries.

Public API:

- :class:`ManifestFormatError` — malformed manifest input (bad JSON,
  missing/extra fields, invalid or duplicate sample names).
- :class:`ManifestEntry` — one parsed manifest row: a sample name, its
  single-sample VCF path and the 1-based manifest line number.
- :func:`parse_manifest` — parse and validate manifest JSON Lines.
- :class:`SampleCall` — one sample's genotype call at a variant.
- :class:`VariantSummary` — one merged variant with its calling samples
  in manifest order.
- :func:`summarize_variants` — validate, normalize and merge per-sample
  calls into sorted summaries.

Each manifest line is a JSON object with exactly the string fields
``sample`` and ``vcf``.  Every VCF holds exactly one sample matching the
manifest name; records carry exactly one plain sequence ALT, FILTER
``PASS`` or ``.``, and a FORMAT with ``GT`` (``0/1`` or ``1/1``), ``DP``
(a non-negative integer) and ``AD`` (two non-negative integers summing
to ``DP``).  Records are minimized and left-aligned against the
reference with the same rules as
:func:`~genome_variant.vcf.normalize_vcf` and merged by
``CHROM``, ``POS``, ``REF`` and ``ALT``.
"""

from __future__ import annotations

import json
import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from .sequence_io import SequenceRecord
from .vcf import (
    VcfFile,
    VcfFormatError,
    VcfRecord,
    _build_reference,
    _is_plain_allele,
    _normalize_against,
    _record_error,
)

__all__ = [
    "ManifestFormatError",
    "ManifestEntry",
    "parse_manifest",
    "SampleCall",
    "VariantSummary",
    "summarize_variants",
]

_MANIFEST_FIELDS = ("sample", "vcf")
_GENOTYPES = ("0/1", "1/1")
_REQUIRED_FORMAT_KEYS = ("GT", "DP", "AD")


class ManifestFormatError(ValueError):
    """Malformed manifest input (bad JSON, fields or duplicate samples)."""


@dataclass(frozen=True)
class ManifestEntry:
    """One manifest row.

    ``vcf`` is the VCF path resolved against the manifest's base
    directory when relative (``-`` for standard input is kept as-is);
    ``line_number`` is the 1-based manifest line the row came from.
    """

    sample: str
    vcf: str
    line_number: int


def parse_manifest(
    lines: Iterable[str],
    *,
    source: str = "<manifest>",
    base_dir: str = "",
) -> list[ManifestEntry]:
    """Parse and validate manifest JSON Lines from *lines*.

    Blank lines are skipped.  Every other line must be a JSON object
    with exactly the string fields ``sample`` and ``vcf``; sample names
    must be non-empty, free of tab characters and globally unique.
    Relative ``vcf`` paths are resolved against *base_dir* (an empty
    *base_dir* leaves them relative to the current working directory).
    Raises :class:`ManifestFormatError` with the 1-based line number.
    """
    entries: list[ManifestEntry] = []
    seen: set[str] = set()
    for number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ManifestFormatError(
                f"{source}:{number}: invalid JSON: {exc}"
            ) from None
        if not isinstance(row, dict) or set(row) != set(_MANIFEST_FIELDS):
            raise ManifestFormatError(
                f"{source}:{number}: expected a JSON object with exactly "
                "the fields 'sample' and 'vcf'"
            )
        sample = row["sample"]
        vcf = row["vcf"]
        if not isinstance(sample, str) or not isinstance(vcf, str):
            raise ManifestFormatError(
                f"{source}:{number}: 'sample' and 'vcf' must be strings"
            )
        if not sample:
            raise ManifestFormatError(f"{source}:{number}: empty sample name")
        if "\t" in sample:
            raise ManifestFormatError(
                f"{source}:{number}: sample name {sample!r} contains a tab"
            )
        if sample in seen:
            raise ManifestFormatError(
                f"{source}:{number}: duplicate sample name {sample!r}"
            )
        seen.add(sample)
        if vcf != "-" and base_dir and not os.path.isabs(vcf):
            vcf = os.path.join(base_dir, vcf)
        entries.append(ManifestEntry(sample, vcf, number))
    return entries


@dataclass(frozen=True)
class SampleCall:
    """One sample's call at a variant: genotype, depth and allelic depths."""

    sample: str
    gt: str
    dp: int
    ad: tuple[int, int]


@dataclass(frozen=True)
class VariantSummary:
    """One merged variant with its calling samples in manifest order."""

    chrom: str
    pos: int
    ref: str
    alt: str
    samples: tuple[SampleCall, ...]

    @property
    def sample_count(self) -> int:
        """The number of calling samples."""
        return len(self.samples)

    @property
    def allele_count(self) -> int:
        """The ALT allele count: one per ``0/1`` call, two per ``1/1``."""
        return sum(1 if call.gt == "0/1" else 2 for call in self.samples)

    @property
    def depth(self) -> int:
        """The summed DP over the calling samples."""
        return sum(call.dp for call in self.samples)

    def to_dict(self) -> dict:
        """The fixed-key JSON-ready mapping for this summary."""
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "sample_count": self.sample_count,
            "allele_count": self.allele_count,
            "depth": self.depth,
            "samples": [
                {
                    "sample": call.sample,
                    "gt": call.gt,
                    "dp": call.dp,
                    "ad": [call.ad[0], call.ad[1]],
                }
                for call in self.samples
            ],
        }


def _header_error(source: str, line_number: int, message: str) -> VcfFormatError:
    if source:
        return VcfFormatError(f"{source}:{line_number}: {message}")
    return VcfFormatError(f"line {line_number}: {message}")


def _nonnegative_int(text: str, source: str, record: VcfRecord, field: str) -> int:
    if not text or any(char < "0" or char > "9" for char in text):
        raise _record_error(
            source, record, f"{field} {text!r} is not a non-negative integer"
        )
    return int(text)


def _extract_call(source: str, record: VcfRecord) -> tuple[str, int, tuple[int, int]]:
    """Validate one record's genotype fields and return ``(gt, dp, ad)``."""
    if len(record.alt) != 1 or not _is_plain_allele(record.alt[0]):
        raise _record_error(
            source, record, "expected exactly one plain sequence ALT allele"
        )
    if record.filter not in ("PASS", "."):
        raise _record_error(
            source, record, f"FILTER {record.filter!r} is not PASS or ."
        )
    if record.format_text is None:
        raise _record_error(source, record, "missing FORMAT column")
    keys = record.format_text.split(":")
    for required in _REQUIRED_FORMAT_KEYS:
        if required not in keys:
            raise _record_error(
                source, record, f"FORMAT does not contain {required}"
            )
    values = record.sample_text[0].split(":") if record.sample_text else []
    deepest = max(keys.index(required) for required in _REQUIRED_FORMAT_KEYS)
    if len(values) <= deepest:
        raise _record_error(
            source, record, "sample column does not provide GT, DP and AD values"
        )

    gt = values[keys.index("GT")]
    if gt not in _GENOTYPES:
        raise _record_error(source, record, f"GT {gt!r} is not 0/1 or 1/1")
    dp = _nonnegative_int(values[keys.index("DP")], source, record, "DP")
    ad_text = values[keys.index("AD")]
    ad_parts = ad_text.split(",")
    if len(ad_parts) != 2:
        raise _record_error(
            source, record, f"AD {ad_text!r} does not have exactly two values"
        )
    ad = (
        _nonnegative_int(ad_parts[0], source, record, "AD"),
        _nonnegative_int(ad_parts[1], source, record, "AD"),
    )
    if ad[0] + ad[1] != dp:
        raise _record_error(
            source,
            record,
            f"AD values {ad[0]},{ad[1]} do not sum to DP {dp}",
        )
    return gt, dp, ad


def summarize_variants(
    documents: Iterable[tuple[str, VcfFile]],
    reference: "Mapping[str, str] | Iterable[SequenceRecord]",
) -> list[VariantSummary]:
    """Merge per-sample calls into sorted variant summaries.

    *documents* yields ``(sample, VcfFile)`` pairs in manifest order;
    each document must declare exactly that one sample.  Every record
    must carry exactly one plain sequence ALT, FILTER ``PASS`` or ``.``
    and a FORMAT with ``GT`` (``0/1`` or ``1/1``), ``DP`` (a
    non-negative integer) and ``AD`` (two non-negative integers summing
    to ``DP``).  Records are minimized and left-aligned against
    *reference* (a CHROM-to-sequence mapping or an iterable of
    :class:`~genome_variant.sequence_io.SequenceRecord`) with the same
    rules as :func:`~genome_variant.vcf.normalize_vcf`, then merged by
    ``CHROM``, ``POS``, ``REF`` and ``ALT``; two records of one sample
    normalizing to the same key are a data error.  Summaries are sorted
    by reference record order, then ``POS``, ``REF`` and ``ALT``.

    Raises :class:`ManifestFormatError` for a repeated sample name,
    :class:`~genome_variant.vcf.VcfFormatError` for record or genotype
    violations and :class:`~genome_variant.vcf.ReferenceMismatchError`
    for an absent/duplicate CHROM or a REF that disagrees with the
    reference.
    """
    if isinstance(reference, Mapping):
        sequences = dict(reference)
        duplicates: set[str] = set()
    else:
        sequences, duplicates = _build_reference(reference)
    chrom_order = {name: index for index, name in enumerate(sequences)}

    merged: dict[tuple[str, int, str, str], list[SampleCall]] = {}
    seen_samples: set[str] = set()
    for sample, document in documents:
        if sample in seen_samples:
            raise ManifestFormatError(f"duplicate sample name {sample!r}")
        seen_samples.add(sample)
        source = document.source
        if document.header.samples != (sample,):
            found = (
                ", ".join(repr(name) for name in document.header.samples)
                if document.header.samples
                else "none"
            )
            raise _header_error(
                source,
                len(document.header.meta_lines) + 1,
                f"expected exactly one sample named {sample!r}, found {found}",
            )

        calls = [_extract_call(source, record) for record in document.records]
        normalized = [
            _normalize_against(record, sequences, duplicates, source)
            for record in document.records
        ]
        sample_keys: set[tuple[str, int, str, str]] = set()
        for record, (gt, dp, ad) in zip(normalized, calls):
            key = (record.chrom, record.pos, record.ref, record.alt[0])
            if key in sample_keys:
                raise _record_error(
                    source,
                    record,
                    f"duplicate normalized variant CHROM {record.chrom} "
                    f"POS {record.pos} REF {record.ref} ALT {record.alt[0]} "
                    f"in sample {sample!r}",
                )
            sample_keys.add(key)
            merged.setdefault(key, []).append(SampleCall(sample, gt, dp, ad))

    summaries = [
        VariantSummary(chrom, pos, ref, alt, tuple(sample_calls))
        for (chrom, pos, ref, alt), sample_calls in merged.items()
    ]
    summaries.sort(
        key=lambda summary: (
            chrom_order[summary.chrom],
            summary.pos,
            summary.ref,
            summary.alt,
        )
    )
    return summaries
