"""Multi-sample variant summarization.

Public API:

- :class:`ManifestError` — malformed manifest text, bad fields, duplicate
  sample names, or a sample VCF that fails structural/genotype validation.
- :func:`read_manifest` — parse a JSON Lines manifest into
  :class:`ManifestEntry` objects in file order.
- :func:`summarize_variants` — read one single-sample VCF per manifest
  entry, normalize the calls against a reference and merge equivalent
  calls into :class:`VariantSummary` objects.

The manifest is JSON Lines: each non-empty line is a JSON object with a
non-empty, tab-free, globally unique ``sample`` and a ``vcf`` path to a
single-sample VCF.  Relative VCF paths in a file manifest resolve
against the manifest's own directory; relative paths in a stream manifest
resolve against the current working directory.

Each VCF must declare exactly one sample, whose name matches its manifest
``sample``; every record must carry a single ordinary sequence ALT, a
``PASS``/``.`` FILTER and the fields ``GT``, ``DP`` and ``AD``.  Genotypes
are restricted to ``0/1`` and ``1/1``; DP is a non-negative integer and
the two AD values are non-negative integers summing to DP.
"""

from __future__ import annotations

import io
import json
import os
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from . import _callsets
from .sequence_io import SequenceRecord
from .vcf import (
    VcfFormatError,
    VcfRecord,
    read_vcf,
)

__all__ = [
    "ManifestError",
    "ManifestEntry",
    "VariantSummary",
    "SampleCall",
    "read_manifest",
    "summarize_variants",
    "render_summaries",
]


class ManifestError(ValueError):
    """The manifest or one of its single-sample VCFs is invalid."""


@dataclass(frozen=True)
class ManifestEntry:
    """One manifest row: its 1-based line, a unique ``sample`` and ``vcf``."""

    line_number: int
    sample: str
    vcf: str


@dataclass(frozen=True)
class VariantSummary:
    """One merged variant line.

    ``ref_index`` is the zero-based position of the variant's CHROM in the
    reference input.  ``samples`` holds one :class:`SampleCall` per actual
    caller in manifest order.
    """

    ref_index: int
    chrom: str
    pos: int
    ref: str
    alt: str
    samples: tuple[SampleCall, ...]

    @property
    def sample_count(self) -> int:
        """Number of samples with the variant called."""
        return len(self.samples)

    @property
    def ref_order_key(self) -> tuple:
        """Sorting key: reference order, then POS, REF, ALT."""
        return (self.ref_index, self.pos, self.ref, self.alt)

    @property
    def allele_count(self) -> int:
        """Sum of alt alleles (one for ``0/1``, two for ``1/1``)."""
        return sum(call.alleles for call in self.samples)

    @property
    def depth(self) -> int:
        """Sum of the per-sample DP values."""
        return sum(call.dp for call in self.samples)


@dataclass(frozen=True)
class SampleCall:
    """One sample's call at a merged variant: ``sample``, GT, DP, AD."""

    sample: str
    gt: str
    dp: int
    ad: tuple[int, int]

    @property
    def alleles(self) -> int:
        """Alt allele count for this genotype (``0/1`` -> 1, ``1/1`` -> 2)."""
        return _callsets.GENOTYPE_ALLELES[self.gt]


def read_manifest(
    source: "str | os.PathLike[str] | io.TextIOBase",
) -> tuple[ManifestEntry, ...]:
    """Parse a whole JSON Lines manifest from *source*.

    *source* is either a filesystem path (``str`` or :class:`os.PathLike`)
    or an already-open text stream.  Blank lines are ignored; every other
    line must decode to a JSON object containing exactly the string fields
    ``sample`` (non-empty, tab-free, unique) and ``vcf``.

    Relative ``vcf`` values resolve against the manifest file's directory
    for a path source, or against the current working directory for a
    stream source.

    Raises :class:`ManifestError` carrying the 1-based manifest line
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
        entries: list[ManifestEntry] = []

        for number, raw_line in enumerate(stream, start=1):
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ManifestError(
                    f"{source_name}:{number}: invalid JSON: {exc.msg}"
                ) from None
            if not isinstance(value, dict):
                raise ManifestError(
                    f"{source_name}:{number}: manifest line must be a JSON object"
                )
            if set(value) != {"sample", "vcf"}:
                raise ManifestError(
                    f"{source_name}:{number}: manifest line must contain exactly "
                    "the fields sample and vcf"
                )
            sample = value["sample"]
            vcf = value["vcf"]
            if not isinstance(sample, str) or not sample or "\t" in sample:
                raise ManifestError(
                    f"{source_name}:{number}: sample must be a non-empty string "
                    "without tab characters"
                )
            if not isinstance(vcf, str) or not vcf:
                raise ManifestError(
                    f"{source_name}:{number}: vcf must be a non-empty string path"
                )
            if sample in seen:
                raise ManifestError(
                    f"{source_name}:{number}: duplicate sample {sample!r}"
                )
            seen.add(sample)
            resolved = vcf if os.path.isabs(vcf) else os.path.join(base_dir, vcf)
            entries.append(ManifestEntry(number, sample, resolved))

        return tuple(entries)
    finally:
        if close:
            stream.close()


def _validate_and_extract(
    record: VcfRecord, source: str
) -> tuple[str, int, tuple[int, int]]:
    """Validate one record and return ``(gt, dp, ad)`` from its sample text."""
    _callsets.validate_record(record, source, ManifestError)
    data = _callsets.extract_sample_data(
        record,
        source,
        ManifestError,
        required=("GT", "DP", "AD"),
        missing_message=lambda missing: (
            "FORMAT must contain GT, DP and AD; missing " + ",".join(missing)
        ),
    )
    gt = _callsets.extract_gt(data, record, source, ManifestError)

    try:
        dp = _callsets.non_negative_int(data["DP"])
    except ValueError:
        raise _callsets.record_error(
            source,
            record,
            f"DP must be a non-negative integer, got {data['DP']!r}",
            ManifestError,
        ) from None

    ad_fields = data["AD"].split(",")
    if len(ad_fields) != 2:
        raise _callsets.record_error(
            source,
            record,
            f"AD must contain exactly two values, got {len(ad_fields)}",
            ManifestError,
        )
    ad: list[int] = []
    for field in ad_fields:
        try:
            ad.append(_callsets.non_negative_int(field))
        except ValueError:
            raise _callsets.record_error(
                source,
                record,
                f"AD values must be non-negative integers, got {field!r}",
                ManifestError,
            ) from None
    if ad[0] + ad[1] != dp:
        raise _callsets.record_error(
            source,
            record,
            f"AD values must sum to DP ({dp}), got {ad[0]}+{ad[1]}",
            ManifestError,
        )

    return gt, dp, (ad[0], ad[1])


def summarize_variants(
    manifest: "str | os.PathLike[str] | io.TextIOBase | Sequence[ManifestEntry]",
    reference: "str | os.PathLike[str] | io.TextIOBase | Iterable[SequenceRecord]",
) -> tuple[VariantSummary, ...]:
    """Summarize the manifest's VCFs into merged variant lines.

    *manifest* is a manifest path/stream or an already-parsed sequence of
    :class:`ManifestEntry`; *reference* is a FASTA path/stream or an
    iterable of :class:`~genome_variant.sequence_io.SequenceRecord`.
    Relative VCF paths are resolved as described for :func:`read_manifest`.

    Every VCF is read, structurally and genotypically validated and
    normalized against the reference; equivalent calls (same normalized
    ``CHROM``, ``POS``, ``REF``, ``ALT``) are then merged across samples.
    Two records in one VCF normalizing to the same key are a data error.

    Returns summaries ordered by reference record order, then POS, REF and
    ALT; each summary's ``samples`` follow manifest order.

    Raises :class:`ManifestError` (message carries the manifest line or
    the VCF source and record line), :class:`VcfFormatError` /
    :class:`ReferenceMismatchError` for reference problems, and
    :class:`OSError` for file access failures.
    """
    if isinstance(manifest, (str, os.PathLike, io.IOBase)):
        entries = read_manifest(manifest)
    else:
        entries = tuple(manifest)

    reference_records = _callsets.load_reference(reference)

    # Reference order for sorting.  Every record's CHROM is checked during
    # normalization, so a later lookup only sees present (and unique) CHROMs.
    order = _callsets.reference_order(reference_records)

    # key (chrom, pos, ref, alt) -> list of SampleCall in manifest order
    calls: dict[tuple[str, int, str, str], list[SampleCall]] = {}
    key_order: list[tuple[str, int, str, str]] = []

    for entry in entries:
        source_label = entry.vcf
        try:
            document = read_vcf(entry.vcf)
        except VcfFormatError as exc:
            raise ManifestError(str(exc)) from None

        # Single sample, named exactly as the manifest entry.
        sample = _callsets.single_sample(document, ManifestError)
        if sample != entry.sample:
            raise ManifestError(
                f"{source_label}:1: sample {sample!r} does not match manifest "
                f"sample {entry.sample!r}"
            )

        # Validation first (extract GT/DP/AD per record), then a single
        # reference-aware normalization of the whole document.  Records
        # map one-to-one and keep their order through normalization, so
        # the extracted data pairs back by index.
        extracted: list[tuple[str, int, tuple[int, int]]] = [
            _validate_and_extract(record, source_label) for record in document.records
        ]

        keys = _callsets.normalized_keys(document, reference_records, ManifestError)

        for key, (gt, dp, ad) in zip(keys, extracted):
            if key not in calls:
                calls[key] = []
                key_order.append(key)
            calls[key].append(SampleCall(entry.sample, gt, dp, ad))

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


def render_summaries(summaries: Sequence[VariantSummary]) -> str:
    """Serialize *summaries* to compact UTF-8 JSON Lines text.

    Each line contains, in order, ``chrom``, ``pos``, ``ref``, ``alt``,
    ``sample_count``, ``allele_count``, ``depth`` and ``samples``; each
    sample item contains ``sample``, ``gt``, ``dp`` and ``ad`` (a
    two-element list).  Non-ASCII text is emitted verbatim and lines end
    with ``"\\n"``; an empty result is the empty string and a non-empty
    result ends with exactly one newline.
    """
    lines = []
    for summary in summaries:
        payload = {
            "chrom": summary.chrom,
            "pos": summary.pos,
            "ref": summary.ref,
            "alt": summary.alt,
            "sample_count": summary.sample_count,
            "allele_count": summary.allele_count,
            "depth": summary.depth,
            "samples": [
                {
                    "sample": call.sample,
                    "gt": call.gt,
                    "dp": call.dp,
                    "ad": [call.ad[0], call.ad[1]],
                }
                for call in summary.samples
            ],
        }
        lines.append(
            json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
        )
    return "".join(line + "\n" for line in lines)
