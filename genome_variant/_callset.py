"""Shared single-sample call-set semantics for summary and comparison.

Internal-only helper module: :mod:`genome_variant.summary` and
:mod:`genome_variant.comparison` both turn a single-sample
:class:`~genome_variant.vcf.VcfFile` into calls keyed by the normalized
``(CHROM, POS, REF, ALT)`` tuple.  The two public entries differ in
*what else* a record must carry (the manifest summary additionally
requires DP/AD with AD summing to DP) and in the exception type their
contract errors take, but the shared rules live here so a given record
is interpreted identically on both paths:

- exactly one declared sample,
- ``PASS``/``.`` FILTER and exactly one ordinary sequence ALT,
- a FORMAT column without duplicate field IDs, exactly one sample
  column with the same number of fields, and a GT of ``0/1`` or ``1/1``,
- reference-aware minimization/left-alignment via
  :func:`~genome_variant.vcf.normalize_vcf`,
- keys of ``(CHROM, POS, REF, ALT)``; two records in one document
  normalizing to the same key are always a contract error.

Nothing in this module is part of the public API.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

from .sequence_io import SequenceRecord
from .vcf import VcfFile, VcfRecord, _is_plain_allele, normalize_vcf

#: Genotypes accepted from the FORMAT/GT field, with allele contribution
#: (one alt allele for ``0/1``, two for ``1/1``).
GENOTYPES: dict[str, int] = {"0/1": 1, "1/1": 2}

#: A normalized variant key: ``(CHROM, POS, REF, ALT)``.
VariantKey = tuple[str, int, str, str]

#: Builds an entry-specific contract error from a finished message.
ErrorType = Callable[[str], Exception]

#: Builds an entry-specific record error qualified by source and line.
RecordErrorBuilder = Callable[[str, VcfRecord, str], Exception]


@dataclass(frozen=True)
class NormalizedCall:
    """One record after shared validation and reference normalization.

    ``key`` is the normalized ``(CHROM, POS, REF, ALT)`` tuple and
    ``record`` the normalized record (retained for its 1-based source
    line when reporting duplicate keys).  ``dp``/``ad`` are ``None`` on
    the comparison path, which only requires a comparable GT.
    """

    key: VariantKey
    record: VcfRecord
    gt: str
    dp: int | None
    ad: tuple[int, int] | None


def record_error(error_type: ErrorType) -> RecordErrorBuilder:
    """Build the shared ``(source, record, message)`` error constructor.

    Messages mirror the VCF reader's shape: ``source:line: message`` when
    both are known, ``line N: message`` for a record without a source and
    the bare message for a record built directly (line 0).
    """

    def build(source: str, record: VcfRecord, message: str) -> Exception:
        if source and record.line_number:
            return error_type(f"{source}:{record.line_number}: {message}")
        if record.line_number:
            return error_type(f"line {record.line_number}: {message}")
        return error_type(message)

    return build


def require_single_sample(document: VcfFile, error_type: ErrorType) -> str:
    """Validate that *document* declares exactly one sample; return its name.

    The sample-count error is qualified by the document source and line
    1 (the line of the ``#CHROM`` header); an unnamed, directly built
    document is qualified as ``<input>``.
    """
    samples = document.header.samples
    if len(samples) != 1:
        raise error_type(
            f"{document.source or '<input>'}:1: VCF must declare exactly one "
            f"sample, got {len(samples)}"
        )
    return samples[0]


def non_negative_int(text: str) -> int:
    """Parse *text* as a non-negative integer without a sign or spaces."""
    if not text or any(char < "0" or char > "9" for char in text):
        raise ValueError("not a non-negative integer")
    return int(text)


def _validate_record(
    record: VcfRecord,
    source: str,
    error: RecordErrorBuilder,
    require_depth: bool,
) -> tuple[str, int | None, tuple[int, int] | None]:
    """Validate one record; return ``(gt, dp, ad)`` (dp/ad ``None`` unless depth)."""
    if record.filter not in ("PASS", "."):
        raise error(
            source, record, f"FILTER must be PASS or '.', got {record.filter!r}"
        )
    if len(record.alt) != 1:
        raise error(
            source,
            record,
            f"expected exactly one ALT allele, got {len(record.alt)}",
        )
    if not _is_plain_allele(record.alt[0]):
        raise error(
            source, record, "ALT must be a single ordinary sequence allele"
        )
    if record.format_text is None:
        raise error(source, record, "missing FORMAT column")

    format_fields = record.format_text.split(":")
    if len(format_fields) != len(set(format_fields)):
        raise error(source, record, "FORMAT contains duplicate field IDs")
    required = ("GT", "DP", "AD") if require_depth else ("GT",)
    missing = [field for field in required if field not in format_fields]
    if missing:
        if required == ("GT",):
            raise error(source, record, "FORMAT must contain GT")
        raise error(
            source,
            record,
            "FORMAT must contain GT, DP and AD; missing " + ",".join(missing),
        )

    if len(record.sample_text) != 1:
        raise error(source, record, "expected one sample column")
    values = record.sample_text[0].split(":")
    if len(values) != len(format_fields):
        raise error(
            source,
            record,
            f"sample has {len(values)} fields but FORMAT lists {len(format_fields)}",
        )
    data = dict(zip(format_fields, values))

    gt = data["GT"]
    if gt not in GENOTYPES:
        raise error(source, record, f"GT must be 0/1 or 1/1, got {gt!r}")

    if not require_depth:
        return gt, None, None

    try:
        dp = non_negative_int(data["DP"])
    except ValueError:
        raise error(
            source, record, f"DP must be a non-negative integer, got {data['DP']!r}"
        ) from None

    ad_fields = data["AD"].split(",")
    if len(ad_fields) != 2:
        raise error(
            source,
            record,
            f"AD must contain exactly two values, got {len(ad_fields)}",
        )
    ad: list[int] = []
    for field in ad_fields:
        try:
            ad.append(non_negative_int(field))
        except ValueError:
            raise error(
                source,
                record,
                f"AD values must be non-negative integers, got {field!r}",
            ) from None
    if ad[0] + ad[1] != dp:
        raise error(
            source,
            record,
            f"AD values must sum to DP ({dp}), got {ad[0]}+{ad[1]}",
        )

    return gt, dp, (ad[0], ad[1])


def normalized_key(record: VcfRecord) -> VariantKey:
    """The shared normalized key of a normalized single-ALT record."""
    return (record.chrom, record.pos, record.ref, record.alt[0])


def validated_records(
    document: VcfFile,
    reference_records: Sequence[SequenceRecord],
    error_type: ErrorType,
    *,
    require_depth: bool,
) -> tuple[NormalizedCall, ...]:
    """Validate, normalize and de-duplicate the records of one VCF document.

    Every record is validated first (FILTER, single ordinary ALT, FORMAT
    shape and GT — plus DP/AD with ``require_depth=True``); the document
    is then normalized against the reference in one pass and two records
    normalizing to the same :func:`normalized_key` are rejected.  Records
    map one-to-one and keep their order through normalization.

    The single-sample requirement is *not* checked here; callers use
    :func:`require_single_sample` first so entry-specific sample checks
    (e.g. the manifest name match) keep their precedence.  All contract
    errors are built with *error_type*, qualified by the document source
    and 1-based record line; reference problems keep their
    :class:`~genome_variant.vcf.ReferenceMismatchError` /
    :class:`~genome_variant.vcf.VcfFormatError` types.
    """
    error = record_error(error_type)

    # Validate every record first, then normalize the whole document in
    # one pass; records map one-to-one and keep their order, so the
    # extracted genotypes pair back by index.
    extracted = [
        _validate_record(record, document.source, error, require_depth)
        for record in document.records
    ]
    normalized = normalize_vcf(document, reference_records)

    calls: list[NormalizedCall] = []
    seen_keys: set[VariantKey] = set()
    for record, (gt, dp, ad) in zip(normalized.records, extracted):
        key = normalized_key(record)
        if key in seen_keys:
            raise error(
                document.source,
                record,
                "duplicate normalized variant "
                f"{record.chrom}:{record.pos}:{record.ref}>{record.alt[0]}",
            )
        seen_keys.add(key)
        calls.append(NormalizedCall(key, record, gt, dp, ad))
    return tuple(calls)


def reference_order(
    reference_records: Sequence[SequenceRecord],
) -> dict[str, int]:
    """Map each reference CHROM to its index, keeping the first occurrence."""
    order: dict[str, int] = {}
    for index, record in enumerate(reference_records):
        order.setdefault(record.identifier, index)
    return order
