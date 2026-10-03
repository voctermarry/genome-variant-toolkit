"""Pairwise comparison of single-sample variant call sets.

Public API:

- :class:`VariantComparisonError` — a VCF fails the single-sample
  call-set constraints (sample count, FILTER, ALT or GT) or two of its
  records normalize to the same key.
- :class:`VariantDifference` — one differing call: its normalized
  ``CHROM``, ``POS``, ``REF``, ``ALT``, a ``status`` and the per-side
  genotypes (``None`` on the side without the call).
- :class:`VariantComparison` — the comparison outcome: both sample
  names, the four category counts and the sorted differences.
- :func:`compare_variants` — read, validate and reference-normalize a
  baseline and a candidate VCF and compare the resulting call sets.
- :func:`render_comparison` — serialize a comparison to one compact
  JSON object followed by a newline.

Each VCF must declare exactly one sample and every record must carry a
``PASS``/``.`` FILTER, a single ordinary sequence ALT and a ``GT`` of
``0/1`` or ``1/1``.  Records are verified against the reference,
minimized and left-aligned (the :func:`~genome_variant.vcf.normalize_vcf`
semantics) and keyed by ``CHROM``, ``POS``, ``REF``, ``ALT``, so input
order and equivalent representations never affect the outcome.
"""

from __future__ import annotations

import io
import json
import os
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from ._callset import (
    reference_order,
    require_single_sample,
    validated_records,
)
from .sequence_io import SequenceRecord, read_sequences
from .vcf import VcfFile, VcfFormatError, read_vcf

__all__ = [
    "VariantComparisonError",
    "VariantDifference",
    "VariantComparison",
    "compare_variants",
    "render_comparison",
]

#: Difference statuses, mirroring the count fields of the same names.
GENOTYPE_MISMATCH = "genotype_mismatch"
BASELINE_ONLY = "baseline_only"
CANDIDATE_ONLY = "candidate_only"


class VariantComparisonError(ValueError):
    """A VCF violates the comparison constraints or has duplicate keys."""


@dataclass(frozen=True)
class VariantDifference:
    """One key present on both sides with different GT, or on one side only.

    ``status`` is ``"genotype_mismatch"``, ``"baseline_only"`` or
    ``"candidate_only"``; ``baseline_gt``/``candidate_gt`` hold the
    per-side genotype and are ``None`` where the call is absent.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    status: str
    baseline_gt: str | None
    candidate_gt: str | None


@dataclass(frozen=True)
class VariantComparison:
    """The outcome of comparing a baseline and a candidate call set.

    ``matched`` counts keys with the same GT on both sides,
    ``genotype_mismatch`` keys shared with different GT and
    ``baseline_only``/``candidate_only`` keys seen on one side only.
    ``differences`` holds one :class:`VariantDifference` per key in the
    last three categories, sorted by reference record order, then POS,
    REF and ALT.
    """

    baseline_sample: str
    candidate_sample: str
    matched: int
    genotype_mismatch: int
    baseline_only: int
    candidate_only: int
    differences: tuple[VariantDifference, ...]

    @property
    def total(self) -> int:
        """Sum of the four category counts."""
        return (
            self.matched
            + self.genotype_mismatch
            + self.baseline_only
            + self.candidate_only
        )

    @property
    def concordance(self) -> str:
        """``matched / total`` with six decimals; ``"1.000000"`` when empty."""
        total = self.total
        if total == 0:
            return "1.000000"
        return f"{self.matched / total:.6f}"


def _extract_calls(
    document: VcfFile,
    reference_records: Sequence[SequenceRecord],
) -> tuple[str, dict[tuple[str, int, str, str], str]]:
    """Validate and normalize *document*; return its sample and key -> GT map.

    Uses the shared single-sample call-set semantics
    (:mod:`genome_variant._callset`) so FILTER, ALT, FORMAT, GT, reference
    normalization, the normalized key and duplicate-key rejection are
    interpreted exactly as in the summary entry; comparison only
    requires a comparable GT and never DP or AD.
    """
    sample = require_single_sample(document, VariantComparisonError)
    calls = validated_records(
        document, reference_records, VariantComparisonError, require_depth=False
    )
    return sample, {call.key: call.gt for call in calls}


def compare_variants(
    baseline: "str | os.PathLike[str] | io.TextIOBase | VcfFile",
    candidate: "str | os.PathLike[str] | io.TextIOBase | VcfFile",
    reference: "str | os.PathLike[str] | io.TextIOBase | Iterable[SequenceRecord]",
) -> VariantComparison:
    """Compare the baseline and candidate single-sample VCF call sets.

    *baseline* and *candidate* are VCF paths, text streams or parsed
    :class:`~genome_variant.vcf.VcfFile` documents; *reference* is a
    FASTA path/stream or an iterable of
    :class:`~genome_variant.sequence_io.SequenceRecord`.

    Both VCFs are read, validated (exactly one sample; every record
    ``PASS``/``.`` FILTER, single ordinary sequence ALT, GT ``0/1`` or
    ``1/1``) and normalized against the reference — REF verified,
    alleles minimized and left-aligned — then keyed by ``CHROM``,
    ``POS``, ``REF``, ``ALT``.  Two records in one VCF normalizing to
    the same key are a data error.  The result is independent of input
    order and of equivalent variant representations.

    Raises :class:`VariantComparisonError` for sample-count, FILTER,
    ALT, GT or duplicate-key violations, :class:`VcfFormatError` for
    malformed VCF or records running past the reference,
    :class:`~genome_variant.vcf.ReferenceMismatchError` for an
    absent/duplicate CHROM or a REF disagreeing with the reference, and
    :class:`OSError` for file access failures.
    """
    if isinstance(reference, (str, os.PathLike, io.IOBase)):
        reference_records = tuple(read_sequences(reference, format="fasta"))
    else:
        reference_records = tuple(reference)

    baseline_document = (
        baseline if isinstance(baseline, VcfFile) else read_vcf(baseline)
    )
    candidate_document = (
        candidate if isinstance(candidate, VcfFile) else read_vcf(candidate)
    )

    baseline_sample, baseline_calls = _extract_calls(
        baseline_document, reference_records
    )
    candidate_sample, candidate_calls = _extract_calls(
        candidate_document, reference_records
    )

    # Reference record order for sorting; every key's CHROM was checked
    # during normalization, so lookups below always succeed.
    order = reference_order(reference_records)

    matched = 0
    genotype_mismatch = 0
    baseline_only = 0
    candidate_only = 0
    differences: list[VariantDifference] = []

    for key, baseline_gt in baseline_calls.items():
        chrom, pos, ref, alt = key
        candidate_gt = candidate_calls.get(key)
        if candidate_gt is None:
            baseline_only += 1
            differences.append(
                VariantDifference(
                    chrom, pos, ref, alt, BASELINE_ONLY, baseline_gt, None
                )
            )
        elif candidate_gt == baseline_gt:
            matched += 1
        else:
            genotype_mismatch += 1
            differences.append(
                VariantDifference(
                    chrom, pos, ref, alt, GENOTYPE_MISMATCH, baseline_gt, candidate_gt
                )
            )
    for key, candidate_gt in candidate_calls.items():
        if key in baseline_calls:
            continue
        chrom, pos, ref, alt = key
        candidate_only += 1
        differences.append(
            VariantDifference(
                chrom, pos, ref, alt, CANDIDATE_ONLY, None, candidate_gt
            )
        )

    differences.sort(key=lambda d: (order[d.chrom], d.pos, d.ref, d.alt))

    return VariantComparison(
        baseline_sample=baseline_sample,
        candidate_sample=candidate_sample,
        matched=matched,
        genotype_mismatch=genotype_mismatch,
        baseline_only=baseline_only,
        candidate_only=candidate_only,
        differences=tuple(differences),
    )


def render_comparison(comparison: VariantComparison) -> str:
    """Serialize *comparison* to one compact JSON object plus a newline.

    The object contains, in order, ``baseline_sample``,
    ``candidate_sample``, ``total``, ``matched``, ``genotype_mismatch``,
    ``baseline_only``, ``candidate_only``, ``concordance`` and
    ``differences``; each difference contains ``chrom``, ``pos``,
    ``ref``, ``alt``, ``status``, ``baseline_gt`` and ``candidate_gt``
    (``null`` on the absent side).  Non-ASCII text is emitted verbatim
    and the result ends with exactly one ``"\\n"``.
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
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n"
