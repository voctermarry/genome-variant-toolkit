"""Shared single-sample call-set semantics for summary and comparison.

Internal to the package — the public entry points remain
:func:`genome_variant.summary.summarize_variants` and
:func:`genome_variant.comparison.compare_variants`.  This module holds
the rules both entries apply identically to one single-sample VCF:
exactly one declared sample, ``PASS``/``.`` FILTER, a single ordinary
sequence ALT, the FORMAT/sample column structure, the accepted GT
values, reference loading and order, reference-aware normalization
keyed by ``CHROM``, ``POS``, ``REF``, ``ALT`` and the rejection of
duplicate normalized keys.

Every helper takes the calling entry's error type so each side keeps
its documented exception boundary (:class:`~genome_variant.summary.ManifestError`
for the summary, :class:`~genome_variant.comparison.VariantComparisonError`
for the comparison) while the checks, messages and key semantics stay
in lockstep.
"""

from __future__ import annotations

import io
import os
from collections.abc import Callable, Iterable, Sequence

from .sequence_io import SequenceRecord, read_sequences
from .vcf import (
    VcfFile,
    VcfRecord,
    _is_plain_allele,
    normalize_vcf,
)

#: Accepted FORMAT/GT values with their alt-allele contribution.
GENOTYPE_ALLELES = {"0/1": 1, "1/1": 2}

#: Normalized variant key: CHROM, POS, REF, ALT.
VariantKey = tuple[str, int, str, str]


def record_error(
    source: str,
    record: VcfRecord,
    message: str,
    error_type: type[ValueError],
) -> ValueError:
    """Build an *error_type* carrying the source name and 1-based line."""
    if source and record.line_number:
        return error_type(f"{source}:{record.line_number}: {message}")
    if record.line_number:
        return error_type(f"line {record.line_number}: {message}")
    return error_type(message)


def non_negative_int(text: str) -> int:
    """Parse *text* as a non-negative integer without a sign or spaces."""
    if not text or any(char < "0" or char > "9" for char in text):
        raise ValueError("not a non-negative integer")
    return int(text)


def load_reference(
    reference: "str | os.PathLike[str] | io.TextIOBase | Iterable[SequenceRecord]",
) -> tuple[SequenceRecord, ...]:
    """Resolve *reference* (FASTA path/stream or records) to a record tuple."""
    if isinstance(reference, (str, os.PathLike, io.IOBase)):
        return tuple(read_sequences(reference, format="fasta"))
    return tuple(reference)


def reference_order(
    reference_records: Sequence[SequenceRecord],
) -> dict[str, int]:
    """Map each CHROM to its first index in the reference input."""
    order: dict[str, int] = {}
    for index, record in enumerate(reference_records):
        order.setdefault(record.identifier, index)
    return order


def single_sample(document: VcfFile, error_type: type[ValueError]) -> str:
    """Return the document's one declared sample, rejecting any other count."""
    samples = document.header.samples
    if len(samples) != 1:
        raise error_type(
            f"{document.source or '<input>'}:1: VCF must declare exactly one "
            f"sample, got {len(samples)}"
        )
    return samples[0]


def validate_record(
    record: VcfRecord, source: str, error_type: type[ValueError]
) -> None:
    """Check the shared record constraints: single plain ALT, PASS/. FILTER."""
    if len(record.alt) != 1:
        raise record_error(
            source,
            record,
            f"expected exactly one ALT allele, got {len(record.alt)}",
            error_type,
        )
    if not _is_plain_allele(record.alt[0]):
        raise record_error(
            source,
            record,
            "ALT must be a single ordinary sequence allele",
            error_type,
        )
    if record.filter not in ("PASS", "."):
        raise record_error(
            source,
            record,
            f"FILTER must be PASS or '.', got {record.filter!r}",
            error_type,
        )


def extract_sample_data(
    record: VcfRecord,
    source: str,
    error_type: type[ValueError],
    required: Sequence[str],
    missing_message: Callable[[list[str]], str],
) -> dict[str, str]:
    """Validate the FORMAT/sample structure and return field -> value.

    The FORMAT column must be present, free of duplicate IDs and cover
    *required*; the record must carry exactly one sample column whose
    value count matches the FORMAT field count.  *missing_message*
    renders the entry-specific text when required fields are absent.
    """
    if record.format_text is None:
        raise record_error(source, record, "missing FORMAT column", error_type)

    format_fields = record.format_text.split(":")
    if len(format_fields) != len(set(format_fields)):
        raise record_error(
            source, record, "FORMAT contains duplicate field IDs", error_type
        )
    missing = [field for field in required if field not in format_fields]
    if missing:
        raise record_error(source, record, missing_message(missing), error_type)

    if len(record.sample_text) != 1:
        raise record_error(source, record, "expected one sample column", error_type)
    values = record.sample_text[0].split(":")
    if len(values) != len(format_fields):
        raise record_error(
            source,
            record,
            f"sample has {len(values)} fields but FORMAT lists {len(format_fields)}",
            error_type,
        )
    return dict(zip(format_fields, values))


def extract_gt(
    data: dict[str, str],
    record: VcfRecord,
    source: str,
    error_type: type[ValueError],
) -> str:
    """Return the GT from parsed sample *data*, restricted to the shared set."""
    gt = data["GT"]
    if gt not in GENOTYPE_ALLELES:
        raise record_error(
            source, record, f"GT must be 0/1 or 1/1, got {gt!r}", error_type
        )
    return gt


def normalized_keys(
    document: VcfFile,
    reference_records: Sequence[SequenceRecord],
    error_type: type[ValueError],
) -> list[VariantKey]:
    """Normalize *document* and return one key per record, in record order.

    Records are verified against the reference, minimized and
    left-aligned (the :func:`~genome_variant.vcf.normalize_vcf`
    semantics) and keyed by ``CHROM``, ``POS``, ``REF``, ``ALT``; two
    records normalizing to the same key are a data error.  Reference
    problems propagate as
    :class:`~genome_variant.vcf.ReferenceMismatchError` /
    :class:`~genome_variant.vcf.VcfFormatError`.
    """
    normalized = normalize_vcf(document, reference_records)
    keys: list[VariantKey] = []
    seen: set[VariantKey] = set()
    for record in normalized.records:
        key = (record.chrom, record.pos, record.ref, record.alt[0])
        if key in seen:
            raise record_error(
                document.source,
                record,
                "duplicate normalized variant "
                f"{record.chrom}:{record.pos}:{record.ref}>{record.alt[0]}",
                error_type,
            )
        seen.add(key)
        keys.append(key)
    return keys
