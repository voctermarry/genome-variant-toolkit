"""Cross-entry consistency regression tests.

These tests pin the public composition contract: for one shared
multi-record reference, one ordered set of FASTQ samples and one set of
calling thresholds, running ``batch-call-variants`` is byte-for-byte
equivalent to running ``call-variants`` once per sample (in manifest
order) and feeding the resulting single-sample VCFs to
``summarize-variants``.

No production behavior is exercised beyond what the three entry points
already promise; the same contract is checked for the default SNV-only
mode and for ``--call-indels``, from different working directories, for
empty inputs, for a failing later sample, for repeated runs and for
different physical chunking of the input reads.  Invalid thresholds,
duplicate sample names, bad references and file-access failures keep
their existing exit codes and exception types on both paths; the
contract is never achieved by relaxing validation.
"""

from __future__ import annotations

import json
import sys
from io import StringIO
from pathlib import Path

import pytest

from genome_variant.batch import BatchManifestError, batch_call_variants
from genome_variant.calling import VariantCallingError, call_variants
from genome_variant.summary import ManifestError, summarize_variants
from genome_variant.sequence_io import SequenceRecord, write_sequences
from genome_variant.cli import main

# -- Fixture data -----------------------------------------------------------

# Three reference records deliberately ordered so that variants appear
# out of reference order across samples; merged output must be sorted by
# reference order rather than encounter order.
C1 = "ACGTACGTACAA"            # SNV sites at POS 5 (A) and POS 8 (T)
C2 = "TAACGTTTCC"             # reverse-strand evidence site at POS 3 (A)
C3 = "ACGTACAAAATCGAC"        # A-run; an inserted A left-aligns to POS 6
REFERENCE_TEXT = f">c1\n{C1}\n>c2\n{C2}\n>c3\n{C3}\n"
REFERENCE_RECORDS = (
    SequenceRecord("c1", C1),
    SequenceRecord("c2", C2),
    SequenceRecord("c3", C3),
)

HIGH = "I"   # Phred 40, above the default quality floor of 20
LOW = "%"    # Phred 4, below the default quality floor of 20

SAMPLE_ORDER = ("s1", "s2", "s3", "s4", "s5")


def _rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _substitute(sequence: str, index: int, base: str) -> str:
    chars = list(sequence)
    chars[index] = base
    return "".join(chars)


def _record(identifier: str, sequence: str, quality: str | None = None):
    if quality is None:
        quality = HIGH * len(sequence)
    return SequenceRecord(
        identifier,
        sequence,
        quality=tuple(ord(char) - 33 for char in quality),
    )


# c2 POS 3 A>T as a forward read; its reverse complement is how the
# reverse-strand evidence appears in the FASTQ.
_C2_VARIANT_READ = _substitute(C2, 2, "T")

# Insertion of one A into the c3 A-run, sampled at three different raw
# breakpoints/offsets; every raw event normalizes to c3:6 C>CA.
_INS_FULL = C3[:6] + "A" + C3[6:]
_INS_SHORT_B = C3[2:10] + "A" + C3[10:15]
_INS_SHORT_C = C3[1:8] + "A" + C3[8:14]


def _sample_records() -> dict[str, list[SequenceRecord]]:
    """Build the shared sample set (fresh record objects for every call)."""
    return {
        # s1: heterozygous c1:5 A>T plus heterozygous c2:3 A>T.
        # Evidence mixes both strands, a low-quality variant base on each
        # strand (excluded from AC and DP at that site), matching wild
        # reads and an unmapped all-N read.
        "s1": [
            _record("s1f1", _substitute(C1, 4, "T")),
            _record("s1f2", _substitute(C1, 4, "T")),
            _record("s1fw", C1),
            _record("s1flq", _substitute(C1, 4, "T"), HIGH * 4 + LOW + HIGH * 7),
            _record("s1r1", _rc(_C2_VARIANT_READ)),
            _record("s1r2", _rc(_C2_VARIANT_READ)),
            # Low quality at the variant column on the reverse strand:
            # original-read index 7 maps to aligned index 2.
            _record("s1rlq", _rc(_C2_VARIANT_READ), HIGH * 7 + LOW + HIGH * 2),
            _record("s1rw1", C2),
            _record("s1rw2", _rc(C2)),
            _record("s1unmapped", "NNNNNNNN"),
        ],
        # s2: c1:5 A>T stays heterozygous while c1:8 T>A is homozygous;
        # ten reads cover both sites so DP is shared per site.
        "s2": [
            *[_record(f"s2a{i}", _substitute(C1, 4, "T")) for i in range(2)],
            *[_record(f"s2b{i}", _substitute(C1, 7, "A")) for i in range(8)],
        ],
        # s3: no variant anywhere (wild reads, one low-quality leading
        # base and an unmapped read).
        "s3": [
            _record("s3w1", C1),
            _record("s3lq", C2, LOW + HIGH * 9),
            _record("s3unmapped", "NNNNNNNN"),
            _record("s3w2", C1),
        ],
        # s4: heterozygous left-aligned insertion c3:6 C>CA; evidence
        # arrives from different raw breakpoints and both strands.
        "s4": [
            _record("s4a", _INS_FULL),
            _record("s4b", _rc(_INS_SHORT_B)),
            _record("s4c", _rc(_INS_SHORT_C)),
            _record("s4w", C3),
        ],
        # s5: the same left-aligned insertion, homozygous.
        "s5": [
            _record("s5a", _INS_FULL),
            _record("s5b", _INS_SHORT_B),
            _record("s5c", _INS_SHORT_C),
        ],
    }


# Expected per-variant summaries under the default thresholds.
EXPECTED_DEFAULT_SNV = [
    {
        "chrom": "c1",
        "pos": 5,
        "ref": "A",
        "alt": "T",
        "sample_count": 2,
        "allele_count": 2,
        "depth": 13,
        "samples": [
            {"sample": "s1", "gt": "0/1", "dp": 3, "ad": [1, 2]},
            {"sample": "s2", "gt": "0/1", "dp": 10, "ad": [8, 2]},
        ],
    },
    {
        "chrom": "c1",
        "pos": 8,
        "ref": "T",
        "alt": "A",
        "sample_count": 1,
        "allele_count": 2,
        "depth": 10,
        "samples": [
            {"sample": "s2", "gt": "1/1", "dp": 10, "ad": [2, 8]},
        ],
    },
    {
        "chrom": "c2",
        "pos": 3,
        "ref": "A",
        "alt": "T",
        "sample_count": 1,
        "allele_count": 1,
        "depth": 4,
        "samples": [
            {"sample": "s1", "gt": "0/1", "dp": 4, "ad": [2, 2]},
        ],
    },
]

EXPECTED_DEFAULT_INSERTION = {
    "chrom": "c3",
    "pos": 6,
    "ref": "C",
    "alt": "CA",
    "sample_count": 2,
    "allele_count": 3,
    "depth": 7,
    "samples": [
        {"sample": "s4", "gt": "0/1", "dp": 4, "ad": [1, 3]},
        {"sample": "s5", "gt": "1/1", "dp": 3, "ad": [0, 3]},
    ],
}


# -- CLI harness ------------------------------------------------------------

def run_cli(argv, stdin: str = ""):
    old_in, old_out, old_err = sys.stdin, sys.stdout, sys.stderr
    sys.stdin = StringIO(stdin)
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    try:
        try:
            code = main(argv)
        except SystemExit as exc:
            code = exc.code
    finally:
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdin, sys.stdout, sys.stderr = old_in, old_out, old_err
    return code, out, err


@pytest.fixture
def workspace(tmp_path):
    """Materialize the reference, per-sample FASTQs and the reads manifest."""
    (tmp_path / "ref.fa").write_text(REFERENCE_TEXT)

    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (tmp_path / "vcfs").mkdir()

    records = _sample_records()
    for sample in SAMPLE_ORDER:
        with open(reads_dir / f"{sample}.fq", "w", encoding="utf-8", newline="") as fh:
            write_sequences(iter(records[sample]), fh, format="fastq")

    manifest = tmp_path / "reads.jsonl"
    manifest.write_text(
        "".join(
            json.dumps({"sample": sample, "reads": f"reads/{sample}.fq"}) + "\n"
            for sample in SAMPLE_ORDER
        )
    )
    return tmp_path


def _batch_argv(workdir: Path, *extra: str) -> list[str]:
    return [
        "batch-call-variants",
        str(workdir / "reads.jsonl"),
        "--reference",
        str(workdir / "ref.fa"),
        *extra,
    ]


def _call_each_sample_then_summarize(workdir: Path, *extra: str):
    """The decomposed path: call-variants per sample, then summarize."""
    for sample in SAMPLE_ORDER:
        code, out, err = run_cli(
            [
                "call-variants",
                str(workdir / "ref.fa"),
                str(workdir / "reads" / f"{sample}.fq"),
                "--sample-name", sample,
                "--output", str(workdir / "vcfs" / f"{sample}.vcf"),
                *extra,
            ]
        )
        assert code == 0, f"call-variants for {sample} failed: {err}"
        assert out == "" and err == ""

    manifest = workdir / "vcfs.jsonl"
    manifest.write_text(
        "".join(
            json.dumps({"sample": sample, "vcf": f"vcfs/{sample}.vcf"}) + "\n"
            for sample in SAMPLE_ORDER
        )
    )
    return run_cli(
        [
            "summarize-variants",
            str(manifest),
            "--reference",
            str(workdir / "ref.fa"),
        ]
    )


def _assert_default_variant_payload(payload, *, indels: bool) -> None:
    expected = list(EXPECTED_DEFAULT_SNV)
    if indels:
        expected.append(EXPECTED_DEFAULT_INSERTION)
    assert payload == expected


# -- The core byte-for-byte contract ----------------------------------------

@pytest.mark.parametrize(
    "extra,indels",
    [
        ([], False),
        (["--call-indels"], True),
    ],
)
def test_batch_output_byte_identical_to_per_sample_calls_then_summary(
    workspace, extra, indels
) -> None:
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace, *extra))
    assert code_batch == 0 and err_batch == ""

    code_sum, out_sum, err_sum = _call_each_sample_then_summarize(workspace, *extra)
    assert code_sum == 0 and err_sum == ""

    # The contract: successful outputs are byte-for-byte identical.
    assert out_batch == out_sum
    assert out_batch.encode("utf-8") == out_sum.encode("utf-8")

    payload = [json.loads(line) for line in out_batch.splitlines()]
    _assert_default_variant_payload(payload, indels=indels)


def test_output_is_structurally_complete_not_just_record_counts(workspace) -> None:
    code, out, err = run_cli(_batch_argv(workspace, "--call-indels"))
    assert code == 0 and err == ""
    payload = [json.loads(line) for line in out.splitlines()]

    # Normalized variant order: reference order (c1, c2, c3), then POS.
    assert [
        tuple(row[key] for key in ("chrom", "pos", "ref", "alt")) for row in payload
    ] == [
        ("c1", 5, "A", "T"),
        ("c1", 8, "T", "A"),
        ("c2", 3, "A", "T"),
        ("c3", 6, "C", "CA"),
    ]

    first = payload[0]
    # Sample order follows the manifest, not per-VCF iteration quirks.
    assert [sample["sample"] for sample in first["samples"]] == ["s1", "s2"]
    # Genotypes, depths and per-allele depths are carried per sample.
    assert [sample["gt"] for sample in first["samples"]] == ["0/1", "0/1"]
    assert [sample["dp"] for sample in first["samples"]] == [3, 10]
    assert [sample["ad"] for sample in first["samples"]] == [[1, 2], [8, 2]]

    hom = payload[1]
    assert hom["samples"][0]["gt"] == "1/1"
    assert hom["samples"][0]["ad"] == [2, 8]

    insertion = payload[3]
    assert [sample["sample"] for sample in insertion["samples"]] == ["s4", "s5"]
    assert [sample["gt"] for sample in insertion["samples"]] == ["0/1", "1/1"]

    # sample_count, allele_count (1 for 0/1, 2 for 1/1) and total depth.
    assert [row["sample_count"] for row in payload] == [2, 1, 1, 2]
    assert [row["allele_count"] for row in payload] == [2, 2, 1, 3]
    assert [row["depth"] for row in payload] == [13, 10, 4, 7]

    # Only actual callers are listed; the no-variant sample s3 appears
    # nowhere even though it sits between callers in the manifest.
    assert all(
        sample["sample"] != "s3"
        for row in payload
        for sample in row["samples"]
    )

    assert out.endswith("\n") and not out.endswith("\n\n")


def test_same_contract_holds_with_shared_non_default_thresholds(workspace) -> None:
    # min-alt-count 3 drops the heterozygous AC=2 sites; min-base-quality
    # 15 admits the low-quality evidence; the homozygous fraction of 0.7
    # re-genotypes the s4 insertion (3/4) as 1/1.
    extra = [
        "--min-base-quality", "15",
        "--min-alt-count", "3",
        "--min-alt-fraction", "0.3",
        "--homozygous-fraction", "0.7",
        "--call-indels",
    ]
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace, *extra))
    assert code_batch == 0 and err_batch == ""
    code_sum, out_sum, err_sum = _call_each_sample_then_summarize(workspace, *extra)
    assert code_sum == 0 and err_sum == ""
    assert out_batch == out_sum

    payload = [json.loads(line) for line in out_batch.splitlines()]
    assert [
        tuple(row[key] for key in ("chrom", "pos", "ref", "alt")) for row in payload
    ] == [
        ("c1", 8, "T", "A"),
        ("c3", 6, "C", "CA"),
    ]
    assert payload[0]["samples"] == [
        {"sample": "s2", "gt": "1/1", "dp": 10, "ad": [2, 8]}
    ]
    assert payload[1]["sample_count"] == 2
    assert [sample["gt"] for sample in payload[1]["samples"]] == ["1/1", "1/1"]
    assert payload[1]["allele_count"] == 4
    assert payload[1]["depth"] == 7


def test_default_snv_mode_does_not_emit_the_indel(workspace) -> None:
    code, out, err = run_cli(_batch_argv(workspace))
    assert code == 0 and err == ""
    payload = [json.loads(line) for line in out.splitlines()]
    assert all(row["chrom"] != "c3" for row in payload)
    assert len(payload) == 3
    # The decomposed path agrees in SNV-only mode too.
    _, out_sum, err_sum = _call_each_sample_then_summarize(workspace)
    assert err_sum == ""
    assert out == out_sum


# -- Relative paths resolve against the manifest directory on both paths ---

def test_relative_manifest_paths_ignore_the_working_directory(
    workspace, monkeypatch
) -> None:
    monkeypatch.chdir("/")
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace))
    assert code_batch == 0 and err_batch == ""
    # The decomposed manifest stores "vcfs/<sample>.vcf" relative paths;
    # summarize must resolve them against the manifest directory as well.
    code_sum, out_sum, err_sum = _call_each_sample_then_summarize(workspace)
    assert code_sum == 0 and err_sum == ""
    assert out_batch == out_sum


def test_both_paths_agree_from_an_unrelated_working_directory(
    workspace, monkeypatch
) -> None:
    monkeypatch.chdir(workspace / "vcfs")
    code_batch, out_batch, _ = run_cli(_batch_argv(workspace, "--call-indels"))
    code_sum, out_sum, _ = _call_each_sample_then_summarize(
        workspace, "--call-indels"
    )
    assert code_batch == code_sum == 0
    assert out_batch == out_sum


# -- Empty inputs -----------------------------------------------------------

def test_empty_manifest_gives_empty_output_on_both_paths(workspace) -> None:
    (workspace / "reads.jsonl").write_text("\n  \n\n")
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace))
    assert code_batch == 0
    assert out_batch == "" and err_batch == ""

    (workspace / "vcfs.jsonl").write_text("\n  \n\n")
    code_sum, out_sum, err_sum = run_cli(
        [
            "summarize-variants",
            str(workspace / "vcfs.jsonl"),
            "--reference",
            str(workspace / "ref.fa"),
        ]
    )
    assert code_sum == 0
    assert out_sum == "" and err_sum == ""
    assert out_batch == out_sum


def test_all_samples_without_variants_give_empty_output_on_both_paths(
    workspace,
) -> None:
    (workspace / "reads.jsonl").write_text(
        json.dumps({"sample": "s3", "reads": "reads/s3.fq"}) + "\n"
    )
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace))
    assert code_batch == 0
    assert out_batch == "" and err_batch == ""

    # Decomposed path: s3's VCF is header-only, and the summary is empty.
    vcf_path = workspace / "vcfs" / "s3.vcf"
    code, out, err = run_cli(
        [
            "call-variants",
            str(workspace / "ref.fa"),
            str(workspace / "reads" / "s3.fq"),
            "--sample-name", "s3",
            "--output", str(vcf_path),
        ]
    )
    assert code == 0 and out == "" and err == ""
    assert vcf_path.read_text().startswith("##fileformat=VCFv4.2\n")
    (workspace / "vcfs.jsonl").write_text(
        json.dumps({"sample": "s3", "vcf": "vcfs/s3.vcf"}) + "\n"
    )
    code_sum, out_sum, err_sum = run_cli(
        [
            "summarize-variants",
            str(workspace / "vcfs.jsonl"),
            "--reference",
            str(workspace / "ref.fa"),
        ]
    )
    assert code_sum == 0
    assert out_sum == "" and err_sum == ""
    assert out_batch == out_sum


# -- A failing later sample fails the whole batch atomically ----------------

def test_corrupt_later_fastq_fails_batch_without_partial_output(
    workspace,
) -> None:
    # A trailing sample (after five valid ones) has a truncated FASTQ.
    (workspace / "reads" / "broken.fq").write_text("@broken\nACGT\n+broken\nII\n")
    manifest = workspace / "reads.jsonl"
    manifest.write_text(
        manifest.read_text()
        + json.dumps({"sample": "s6", "reads": "reads/broken.fq"}) + "\n"
    )
    target = workspace / "out.jsonl"
    target.write_text("PREEXISTING\n")

    code, out, err = run_cli(_batch_argv(workspace, "--output", str(target)))
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    # The pre-existing destination is untouched and no temp summary stays.
    assert target.read_text() == "PREEXISTING\n"
    assert list(workspace.glob(".batch-call-variants-*")) == []


def test_corrupt_later_fastq_fails_the_per_sample_call_the_same_way(
    workspace,
) -> None:
    # The decomposed path reaches that read through call-variants: the
    # same malformed input must fail there with the same usage exit code,
    # so batch consistency cannot come from relaxed validation.
    (workspace / "reads" / "broken.fq").write_text("@broken\nACGT\n+broken\nII\n")
    code, out, err = run_cli(
        [
            "call-variants",
            str(workspace / "ref.fa"),
            str(workspace / "reads" / "broken.fq"),
            "--sample-name", "s6",
        ]
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


# -- Determinism: repetition and physical read chunking ---------------------

def test_repeated_runs_produce_identical_bytes(workspace) -> None:
    argv = _batch_argv(workspace, "--call-indels")
    first = run_cli(argv)
    second = run_cli(argv)
    assert first == second
    assert first[0] == 0 and first[2] == ""
    # And the decomposed path is repeatable as well.
    decomposed_first = _call_each_sample_then_summarize(workspace, "--call-indels")
    decomposed_second = _call_each_sample_then_summarize(workspace, "--call-indels")
    assert decomposed_first == decomposed_second
    assert first[1] == decomposed_first[1]


def _wrapped_fastq_text(records, chunk: int) -> str:
    """Serialize records with sequence/quality wrapped every *chunk* bases.

    The FASTQ reader accumulates wrapped sequence and quality segments
    into the same logical records, so this varies the physical line
    chunking seen during iteration without changing a single record.
    """
    lines: list[str] = []
    for record in records:
        quality = "".join(chr(score + 33) for score in record.quality)
        lines.append("@" + record.identifier)
        for start in range(0, len(record.sequence), chunk):
            lines.append(record.sequence[start : start + chunk])
        lines.append("+" + record.identifier)
        for start in range(0, len(quality), chunk):
            lines.append(quality[start : start + chunk])
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("chunk", [1, 3, 7])
@pytest.mark.parametrize("extra", [[], ["--call-indels"]])
def test_output_is_independent_of_physical_read_chunking(
    workspace, chunk, extra
) -> None:
    # Rewrite every sample FASTQ with its records wrapped across multiple
    # physical lines; logical records, order and thresholds are unchanged.
    wrapped_dir = workspace / "reads_wrapped"
    wrapped_dir.mkdir()
    for sample, records in _sample_records().items():
        (wrapped_dir / f"{sample}.fq").write_text(
            _wrapped_fastq_text(records, chunk)
        )
    (workspace / "reads_wrapped.jsonl").write_text(
        "".join(
            json.dumps({"sample": sample, "reads": f"reads_wrapped/{sample}.fq"})
            + "\n"
            for sample in SAMPLE_ORDER
        )
    )

    _, out_normal, err_normal = run_cli(_batch_argv(workspace, *extra))
    assert err_normal == ""

    code_wrapped, out_wrapped, err_wrapped = run_cli(
        [
            "batch-call-variants",
            str(workspace / "reads_wrapped.jsonl"),
            "--reference",
            str(workspace / "ref.fa"),
            *extra,
        ]
    )
    assert code_wrapped == 0 and err_wrapped == ""
    assert out_wrapped == out_normal

    # The decomposed path over the wrapped reads produces the same bytes
    # (call-variants parses the wrapped FASTQ identically).
    for sample in SAMPLE_ORDER:
        code, out, err = run_cli(
            [
                "call-variants",
                str(workspace / "ref.fa"),
                str(wrapped_dir / f"{sample}.fq"),
                "--sample-name", sample,
                "--output", str(workspace / "vcfs" / f"{sample}.vcf"),
                *extra,
            ]
        )
        assert code == 0, err
        assert out == "" and err == ""
    manifest = workspace / "vcfs.jsonl"
    manifest.write_text(
        "".join(
            json.dumps({"sample": sample, "vcf": f"vcfs/{sample}.vcf"}) + "\n"
            for sample in SAMPLE_ORDER
        )
    )
    code_sum, out_sum, err_sum = run_cli(
        [
            "summarize-variants",
            str(manifest),
            "--reference",
            str(workspace / "ref.fa"),
        ]
    )
    assert code_sum == 0 and err_sum == ""
    assert out_sum == out_normal


# -- Validation parity: failures stay failures on both paths ----------------

@pytest.mark.parametrize(
    "extra",
    [
        ["--min-base-quality", "94"],
        ["--min-alt-count", "0"],
        ["--min-alt-fraction", "1.5"],
        ["--homozygous-fraction", "nan"],
        ["--min-alt-fraction", "0.8", "--homozygous-fraction", "0.5"],
        ["--max-indel-length", "10"],
    ],
)
def test_invalid_thresholds_exit_2_on_both_paths(workspace, extra) -> None:
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace, *extra))
    assert code_batch == 2
    assert out_batch == ""
    assert err_batch.count("\n") == 1

    # The same threshold combination must fail a single-sample call too.
    code_call, out_call, err_call = run_cli(
        [
            "call-variants",
            str(workspace / "ref.fa"),
            str(workspace / "reads" / "s1.fq"),
            "--sample-name", "s1",
            *extra,
        ]
    )
    assert code_call == 2
    assert out_call == ""
    assert err_call.count("\n") == 1


def test_duplicate_sample_names_rejected_on_both_manifests(workspace) -> None:
    (workspace / "reads.jsonl").write_text(
        json.dumps({"sample": "s1", "reads": "reads/s1.fq"}) + "\n"
        + json.dumps({"sample": "s1", "reads": "reads/s2.fq"}) + "\n"
    )
    code, out, err = run_cli(_batch_argv(workspace))
    assert code == 2 and out == ""
    assert err.count("\n") == 1 and "duplicate" in err

    # The summarize manifest shares the duplicate-sample rule (detected
    # while parsing the manifest, before any VCF is opened).
    (workspace / "vcfs.jsonl").write_text(
        '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
        '{"sample":"s1","vcf":"vcfs/s2.vcf"}\n'
    )
    code, out, err = run_cli(
        [
            "summarize-variants",
            str(workspace / "vcfs.jsonl"),
            "--reference",
            str(workspace / "ref.fa"),
        ]
    )
    assert code == 2 and out == ""
    assert err.count("\n") == 1 and "duplicate" in err


@pytest.mark.parametrize(
    "reference_text",
    [
        "",                        # empty reference
        ">c1\nACGT\n>c1\nTTTT\n",  # duplicate identifier
    ],
)
def test_bad_references_exit_2_on_both_paths(workspace, reference_text) -> None:
    bad_ref = workspace / "bad.fa"
    bad_ref.write_text(reference_text)

    code_batch, out_batch, err_batch = run_cli(
        [
            "batch-call-variants",
            str(workspace / "reads.jsonl"),
            "--reference", str(bad_ref),
        ]
    )
    assert code_batch == 2
    assert out_batch == ""
    assert err_batch.count("\n") == 1

    code_call, out_call, err_call = run_cli(
        [
            "call-variants",
            str(bad_ref),
            str(workspace / "reads" / "s1.fq"),
            "--sample-name", "s1",
        ]
    )
    assert code_call == 2
    assert out_call == ""
    assert err_call.count("\n") == 1


def test_missing_reads_exit_1_and_missing_vcfs_exit_1(workspace) -> None:
    (workspace / "reads.jsonl").write_text(
        json.dumps({"sample": "s1", "reads": "reads/missing.fq"}) + "\n"
    )
    code_batch, out_batch, err_batch = run_cli(_batch_argv(workspace))
    assert code_batch == 1
    assert out_batch == ""
    assert err_batch.count("\n") == 1

    (workspace / "vcfs.jsonl").write_text(
        json.dumps({"sample": "s1", "vcf": "vcfs/missing.vcf"}) + "\n"
    )
    code_sum, out_sum, err_sum = run_cli(
        [
            "summarize-variants",
            str(workspace / "vcfs.jsonl"),
            "--reference",
            str(workspace / "ref.fa"),
        ]
    )
    assert code_sum == 1
    assert out_sum == ""
    assert err_sum.count("\n") == 1


def test_public_exception_types_are_unchanged() -> None:
    """The documented exception types survive on both composition sides."""
    # Invalid thresholds are ValueErrors raised before reads are consumed.
    with pytest.raises(ValueError):
        batch_call_variants((), REFERENCE_RECORDS, min_base_quality=94)
    with pytest.raises(ValueError):
        call_variants(REFERENCE_RECORDS, (), min_base_quality=94)

    # An empty reference is a VariantCallingError on both sides.
    with pytest.raises(VariantCallingError):
        batch_call_variants((), ())
    with pytest.raises(VariantCallingError):
        call_variants((), ())

    # Malformed manifests keep their entry-specific error types.
    with pytest.raises(BatchManifestError):
        batch_call_variants(StringIO("{not json}\n"), REFERENCE_RECORDS)
    with pytest.raises(ManifestError):
        summarize_variants(StringIO("{not json}\n"), REFERENCE_RECORDS)

    # Missing files keep propagating as OSError.
    with pytest.raises(OSError):
        batch_call_variants("/no/such/manifest.jsonl", REFERENCE_RECORDS)
    with pytest.raises(OSError):
        summarize_variants(
            StringIO('{"sample":"s1","vcf":"/no/such/s1.vcf"}\n'),
            REFERENCE_RECORDS,
        )
