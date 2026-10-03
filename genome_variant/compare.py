"""Batch-level comparison of two single-sample VCF call sets.

Public API:

- :class:`VariantComparisonError` — a VCF declares the wrong number of
  samples or carries records outside the accepted FILTER/ALT/GT subset,
  or the same VCF contains two records normalizing to one variant key.
- :class:`VariantDifference` — one mismatched variant: shared key,
  status and the genotype seen on each side.
- :class:`VariantComparison` — the sample names, the four category
  counts, the concordance figure and the ordered differences.
- :func:`compare_variants` — read, validate and normalize a baseline and
  a candidate VCF against one reference and compare their calls.
- :func:`render_comparison` — serialize a comparison to compact JSON.

Both VCFs must declare exactly one sample and every record must carry a
single ordinary sequence ALT, a ``PASS``/``.`` FILTER and a GT of
``0/1`` or ``1/1`` (extra FORMAT fields are allowed in any order).  The
records are reference-validated, minimized and left-aligned as in
:func:`genome_variant.vcf.normalize_vcf` before comparison; the key is
the normalized ``CHROM``, ``POS``, ``REF`` and ``ALT``.
"""

from __future__ import annotations

import io
import json
import os
from collections.abc import Iterable
from dataclasses import dataclass
from .sequence_io import SequenceRecord, read_sequences
from .vcf import (
    VcfFile,
    VcfFormatError,
    VcfRecord,
    _is_plain_allele,
    normalize_vcf,
    read_vcf,
)

__all__ = [
    "VariantComparisonError",
    "VariantDifference",
    "VariantComparison",
    "compare_variants",
    "render_comparison",
]

#: Genotypes accepted from the FORMAT/GT field.
_GENOTYPES = frozenset(("0/1", "1/1"))

#: Difference status values, in the order the JSON reports them.
STATUS_GENOTYPE_MISMATCH = "genotype_mismatch"
STATUS_BASELINE_ONLY = "baseline_only"
STATUS_CANDIDATE_ONLY = "candidate_only"

_VcfSource = "str | os.PathLike[str] | io.TextIOBase | VcfFile"
_ReferenceSource = (
    "str | os.PathLike[str] | io.TextIOBase | Iterable[SequenceRecord]"
)


class VariantComparisonError(ValueError):
    """A compared VCF violates the single-sample/record subset or lists a
    normalized variant key more than once."""


@dataclass(frozen=True)
class VariantDifference:
    """One variant that is not an exact match.

    Fields are ordered ``chrom``, ``pos``, ``ref``, ``alt``, ``status``,
    ``baseline_gt`` and ``candidate_gt``; *status* is one of
    ``"genotype_mismatch"``, ``"baseline_only"`` or
    ``"candidate_only"`` and the missing side's genotype is ``None``.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    status: str
    baseline_gt: str | None
    candidate_gt: str | None

    @property
    def key(self) -> tuple[str, int, str, str]:
        """The normalized ``(CHROM, POS, REF, ALT)`` variant key."""
        return (self.chrom, self.pos, self.ref, self.alt)


@dataclass(frozen=True)
class VariantComparison:
    """The result of comparing a baseline call set against a candidate.

    ``total`` equals ``matched + genotype_mismatch + baseline_only +
    candidate_only``; ``concordance`` is ``matched / total`` formatted
    with six decimals (``"1.000000"`` when both sides have no records).
    ``differences`` holds only the non-matched categories.
    """

    baseline_sample: str
    candidate_sample: str
    total: int
    matched: int
    genotype_mismatch: int
    baseline_only: int
    candidate_only: int
    concordance: str
    differences: tuple[VariantDifference, ...]


def _record_error(source: str, record: VcfRecord, message: str) -> VariantComparisonError:
    if source and record.line_number:
        return VariantComparisonError(f"{source}:{record.line_number}: {message}")
    if record.line_number:
        return VariantComparisonError(f"line {record.line_number}: {message}")
    return VariantComparisonError(message)


def _load_vcf(source: "str | os.PathLike[str] | io.TextIOBase | VcfFile") -> VcfFile:
    if isinstance(source, VcfFile):
        return source
    if isinstance(source, (str, os.PathLike, io.IOBase)):
        return read_vcf(source)
    raise TypeError("VCF source must be a text path, a text stream or a VcfFile")


def _extract_gt(record: VcfRecord, source: str) -> str:
    """Validate one record's FILTER/ALT subset and return its GT."""
    if record.filter not in ("PASS", "."):
        raise _record_error(
            source, record, f"FILTER must be PASS or '.', got {record.filter!r}"
        )
    if len(record.alt) != 1:
        raise _record_error(
            source,
            record,
            f"expected exactly one ALT allele, got {len(record.alt)}",
        )
    if not _is_plain_allele(record.alt[0]):
        raise _record_error(
            source, record, "ALT must be a single ordinary sequence allele"
        )
    if record.format_text is None:
        raise _record_error(source, record, "missing FORMAT column")

    format_fields = record.format_text.split(":")
    if len(format_fields) != len(set(format_fields)):
        raise _record_error(source, record, "FORMAT contains duplicate field IDs")
    if "GT" not in format_fields:
        raise _record_error(source, record, "FORMAT must contain GT")
    if len(record.sample_text) != 1:
        raise _record_error(source, record, "expected one sample column")

    values = record.sample_text[0].split(":")
    if len(values) != len(format_fields):
        raise _record_error(
            source,
            record,
            f"sample has {len(values)} fields but FORMAT lists {len(format_fields)}",
        )
    gt = dict(zip(format_fields, values))["GT"]
    if gt not in _GENOTYPES:
        raise _record_error(
            source, record, f"GT must be 0/1 or 1/1, got {gt!r}"
        )
    return gt


def _prepare_side(
    document: VcfFile,
    reference_records: tuple[SequenceRecord, ...],
) -> tuple[str, dict[tuple[str, int, str, str], str]]:
    """Validate, normalize and index one VCF as ``key -> GT``."""
    samples = document.header.samples
    if len(samples) != 1:
        raise VariantComparisonError(
            f"{document.source}:1: VCF must declare exactly one sample, "
            f"got {len(samples)}"
        )
    sample = samples[0]
    source = document.source

    # Validate the accepted subset first, then normalize the whole
    # document; records map one-to-one through normalization, so the
    # extracted GTs pair back by index.
    genotypes = [_extract_gt(record, source) for record in document.records]
    normalized = normalize_vcf(document, reference_records)

    calls: dict[tuple[str, int, str, str], str] = {}
    for record, gt in zip(normalized.records, genotypes):
        key = (record.chrom, record.pos, record.ref, record.alt[0])
        if key in calls:
            raise VariantComparisonError(
                f"{source}:{record.line_number}: duplicate normalized variant "
                f"{record.chrom}:{record.pos}:{record.ref}>{record.alt[0]}"
            )
        calls[key] = gt
    return sample, calls


def compare_variants(
    baseline: _VcfSource,
    candidate: _VcfSource,
    reference: _ReferenceSource,
) -> VariantComparison:
    """Compare the calls of *baseline* and *candidate*.

    *baseline* and *candidate* are each a VCF path/text stream or an
    already-parsed :class:`~genome_variant.vcf.VcfFile`; *reference* is
    a FASTA path/text stream or an iterable of
    :class:`~genome_variant.sequence_io.SequenceRecord`.

    Both VCFs must declare exactly one sample; every record must pass
    the FILTER/single-plain-ALT/``0|1``-or-``1|1`` subset described in
    the module docstring.  Records are normalized against *reference*
    and compared on the normalized ``CHROM``/``POS``/``REF``/``ALT``
    key; identical keys with different GT are genotype mismatches and
    keys seen on only one side are baseline-only or candidate-only.

    Differences are ordered by reference record order, then POS, REF and
    ALT.  Raises :class:`VariantComparisonError` for the single-sample /
    record / duplicate-key constraints, :class:`VcfFormatError` or
    :class:`~genome_variant.vcf.ReferenceMismatchError` for structural
    and reference problems, and :class:`OSError` for access failures.
    """
    baseline_document = _load_vcf(baseline)
    candidate_document = _load_vcf(candidate)

    if isinstance(reference, (str, os.PathLike, io.IOBase)):
        reference_records = tuple(read_sequences(reference, format="fasta"))
    else:
        reference_records = tuple(reference)

    baseline_sample, baseline_calls = _prepare_side(
        baseline_document, reference_records
    )
    candidate_sample, candidate_calls = _prepare_side(
        candidate_document, reference_records
    )

    # Every compared CHROM was resolved during normalization, so this
    # order map only contains present identifiers.
    order: dict[str, int] = {}
    for index, record in enumerate(reference_records):
        order.setdefault(record.identifier, index)

    matched = 0
    genotype_mismatch = 0
    baseline_only = 0
    differences: list[VariantDifference] = []

    for key, baseline_gt in baseline_calls.items():
        chrom, pos, ref, alt = key
        candidate_gt = candidate_calls.get(key)
        if candidate_gt is None:
            baseline_only += 1
            differences.append(
                VariantDifference(
                    chrom, pos, ref, alt, STATUS_BASELINE_ONLY, baseline_gt, None
                )
            )
        elif candidate_gt == baseline_gt:
            matched += 1
        else:
            genotype_mismatch += 1
            differences.append(
                VariantDifference(
                    chrom,
                    pos,
                    ref,
                    alt,
                    STATUS_GENOTYPE_MISMATCH,
                    baseline_gt,
                    candidate_gt,
                )
            )

    candidate_only = 0
    for key, candidate_gt in candidate_calls.items():
        if key in baseline_calls:
            continue
        candidate_only += 1
        chrom, pos, ref, alt = key
        differences.append(
            VariantDifference(
                chrom, pos, ref, alt, STATUS_CANDIDATE_ONLY, None, candidate_gt
            )
        )

    total = matched + genotype_mismatch + baseline_only + candidate_only
    concordance = "1.000000" if total == 0 else f"{matched / total:.6f}"

    differences.sort(
        key=lambda difference: (
            order[difference.chrom],
            difference.pos,
            difference.ref,
            difference.alt,
        )
    )

    return VariantComparison(
        baseline_sample=baseline_sample,
        candidate_sample=candidate_sample,
        total=total,
        matched=matched,
        genotype_mismatch=genotype_mismatch,
        baseline_only=baseline_only,
        candidate_only=candidate_only,
        concordance=concordance,
        differences=tuple(differences),
    )


def render_comparison(comparison: VariantComparison) -> str:
    """Serialize *comparison* to one compact UTF-8 JSON object.

    The top-level fields are emitted in the order ``baseline_sample``,
    ``candidate_sample``, ``total``, ``matched``, ``genotype_mismatch``,
    ``baseline_only``, ``candidate_only``, ``concordance`` and
    ``differences``; each difference item carries ``chrom``, ``pos``,
    ``ref``, ``alt``, ``status``, ``baseline_gt`` and ``candidate_gt``
    in that order.  Non-ASCII text is emitted verbatim and the text ends
    with exactly one ``"\\n"``.
    """
    payload = {
        "baseline_sample": comparison.baseline_sample,
        "candidate_sample": comparison.candidate_sample,
        "total": comparison.total,
        "matched": comparison.matched,
        "genotype_mismatch": comparison.genotype_mismatch,
        "baseline_only": comparison.baseline_only,
        "candidate_only": comparison.candidate_only,
        "concordance": comparison.concordance,
        "differences": [
            {
                "chrom": difference.chrom,
                "pos": difference.pos,
                "ref": difference.ref,
                "alt": difference.alt,
                "status": difference.status,
                "baseline_gt": difference.baseline_gt,
                "candidate_gt": difference.candidate_gt,
            }
            for difference in comparison.differences
        ],
    }
    return (
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n"
    )
