"""Cross-entry consistency tests for the batch calling contract.

``batch-call-variants`` promises to be equivalent to running
``call-variants`` per sample followed by ``summarize-variants`` on the
resulting single-sample VCFs.  These tests exercise both paths from the
same multi-record reference, the same manifest-ordered FASTQ samples and
the same calling thresholds, and require the successful outputs to be
byte-identical — in the default SNV-only mode and with ``--call-indels``.

The fixture input covers plus- and minus-strand evidence, low-quality
bases, an unmapped read, a no-variant sample, heterozygous and
homozygous calls, and a homopolymer deletion whose equivalent raw
placements left-align to one normalized indel.  No command, option or
output field is added here, and no validation is relaxed: the existing
error behaviour of every entry point is asserted to stay as published.
"""

from __future__ import annotations

import itertools
import json
import sys
from io import StringIO

import pytest

from genome_variant.batch import batch_call_variants, read_batch_manifest
from genome_variant.calling import call_variants
from genome_variant.cli import main
from genome_variant.sequence_io import read_sequences
from genome_variant.summary import render_summaries
from genome_variant.vcf import render_vcf

# Two-record reference.  c1 carries two SNV sites (POS 10 C>T and
# POS 20 G>A); c2 carries an SNV site (POS 3 C>T) and a homopolymer run
# of seven A's (POS 8-14) in which a one-base deletion left-aligns to
# the normalized event c2:7 CA>C.
C1 = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"
C2 = "GGCTTACAAAAAAATGCGTAC"
REFERENCE = f">c1\n{C1}\n>c2\n{C2}\n"

ALT10 = C1[:9] + "T" + C1[10:]
ALT20 = C1[:19] + "A" + C1[20:]
ALT10_20 = ALT10[:19] + "A" + ALT10[20:]
DEL_C2 = C2[:10] + C2[11:]  # one A removed from the homopolymer run
S4_HAPLOTYPE = DEL_C2[:2] + "T" + DEL_C2[3:]  # C>T at POS 3 plus the deletion

# Manifest order deliberately differs from sorted/sample-file order.
SAMPLE_ORDER = ("s2", "s1", "s4", "s3")

EXPECTED_SNV = [
    {
        "chrom": "c1", "pos": 10, "ref": "C", "alt": "T",
        "sample_count": 2, "allele_count": 3, "depth": 8,
        "samples": [
            {"sample": "s2", "gt": "1/1", "dp": 4, "ad": [0, 4]},
            {"sample": "s1", "gt": "0/1", "dp": 4, "ad": [2, 2]},
        ],
    },
    {
        "chrom": "c1", "pos": 20, "ref": "G", "alt": "A",
        "sample_count": 1, "allele_count": 1, "depth": 4,
        "samples": [
            {"sample": "s2", "gt": "0/1", "dp": 4, "ad": [2, 2]},
        ],
    },
    {
        "chrom": "c2", "pos": 3, "ref": "C", "alt": "T",
        "sample_count": 1, "allele_count": 2, "depth": 3,
        "samples": [
            {"sample": "s4", "gt": "1/1", "dp": 3, "ad": [0, 3]},
        ],
    },
]

EXPECTED_INDEL = EXPECTED_SNV + [
    {
        "chrom": "c2", "pos": 7, "ref": "CA", "alt": "C",
        "sample_count": 2, "allele_count": 3, "depth": 7,
        "samples": [
            {"sample": "s1", "gt": "0/1", "dp": 4, "ad": [1, 3]},
            {"sample": "s4", "gt": "1/1", "dp": 3, "ad": [0, 3]},
        ],
    },
]


def run(argv, stdin: str = ""):
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


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fastq(identifier, sequence, quality=None):
    if quality is None:
        quality = "I" * len(sequence)
    return f"@{identifier}\n{sequence}\n+{identifier}\n{quality}\n"


def _sample_reads():
    # s1: heterozygous c1:10 from one plus-strand and one minus-strand ALT
    # read against two reference reads; a further ALT read whose variant
    # base is low quality (Phred 0) contributes nothing; three deletion
    # reads against one reference read make a heterozygous c2 indel; one
    # all-N read is unmapped.
    low_quality = "I" * 9 + "!" + "I" * (len(C1) - 10)
    s1 = (
        fastq("a1", ALT10)
        + fastq("a2", revcomp(ALT10))
        + fastq("w1", C1)
        + fastq("w2", C1)
        + fastq("lq", ALT10, low_quality)
        + fastq("d1", DEL_C2)
        + fastq("d2", DEL_C2)
        + fastq("d3", DEL_C2)
        + fastq("w3", C2)
        + fastq("un", "N" * 20)
    )
    # s2: homozygous c1:10 (4/4 ALT), heterozygous c1:20 (2/4 ALT).
    s2 = (
        fastq("b1", ALT10)
        + fastq("b2", ALT10)
        + fastq("b3", ALT10_20)
        + fastq("b4", ALT10_20)
    )
    # s3: no variants at all.
    s3 = (
        fastq("e1", C1)
        + fastq("e2", C1)
        + fastq("e3", C1)
        + fastq("e4", C2)
        + fastq("e5", C2)
    )
    # s4: homozygous c2:3 SNV and homozygous c2 deletion on the same
    # reads; the different flanking context still normalizes to the same
    # left-aligned indel event as s1's reads.
    s4 = fastq("f1", S4_HAPLOTYPE) + fastq("f2", S4_HAPLOTYPE) + fastq(
        "f3", S4_HAPLOTYPE
    )
    return {"s1": s1, "s2": s2, "s3": s3, "s4": s4}


@pytest.fixture
def workspace(tmp_path):
    """Reference, per-sample FASTQs and a relative-path reads manifest."""
    (tmp_path / "ref.fa").write_text(REFERENCE)
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    for sample, text in _sample_reads().items():
        (reads_dir / f"{sample}.fq").write_text(text)
    manifest = tmp_path / "m.jsonl"
    manifest.write_text(
        "".join(
            f'{{"sample":"{sample}","reads":"reads/{sample}.fq"}}\n'
            for sample in SAMPLE_ORDER
        )
    )
    return tmp_path


def run_batch(workspace, extra_args=(), manifest="m.jsonl"):
    """The direct path: one ``batch-call-variants`` invocation."""
    return run(
        [
            "batch-call-variants",
            str(workspace / manifest),
            "--reference",
            str(workspace / "ref.fa"),
            *extra_args,
        ]
    )


def run_per_sample(workspace, extra_args=(), samples=SAMPLE_ORDER):
    """The composed path: ``call-variants`` per sample, then summarize."""
    vcf_dir = workspace / "vcfs"
    vcf_dir.mkdir(exist_ok=True)
    for sample in samples:
        code, out, err = run(
            [
                "call-variants",
                str(workspace / "ref.fa"),
                str(workspace / "reads" / f"{sample}.fq"),
                "--sample-name",
                sample,
                *extra_args,
            ]
        )
        assert code == 0, f"call-variants failed for {sample}: {err}"
        (vcf_dir / f"{sample}.vcf").write_text(out)
    manifest = workspace / "vm.jsonl"
    manifest.write_text(
        "".join(
            f'{{"sample":"{sample}","vcf":"vcfs/{sample}.vcf"}}\n'
            for sample in samples
        )
    )
    return run(
        ["summarize-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
    )


def parse_lines(out):
    return [json.loads(line) for line in out.splitlines()]


class TestBatchEqualsCallThenSummarize:
    """The combination contract: both paths, byte-identical output."""

    def test_snv_mode_byte_identical_and_matches_expected(self, workspace):
        batch_code, batch_out, batch_err = run_batch(workspace)
        comp_code, comp_out, comp_err = run_per_sample(workspace)
        assert batch_code == comp_code == 0
        assert batch_err == comp_err == ""
        assert batch_out == comp_out
        # Content acceptance: variant order, sample order, GT/DP/AD and
        # the aggregate counts — not just the number of records.
        assert parse_lines(batch_out) == EXPECTED_SNV
        manifest_index = {sample: i for i, sample in enumerate(SAMPLE_ORDER)}
        for entry in parse_lines(batch_out):
            # Each variant lists only its callers, in manifest order.
            indices = [manifest_index[s["sample"]] for s in entry["samples"]]
            assert indices == sorted(indices)
            assert entry["sample_count"] == len(entry["samples"])
            assert entry["depth"] == sum(s["dp"] for s in entry["samples"])
            assert entry["allele_count"] == sum(
                2 if s["gt"] == "1/1" else 1 for s in entry["samples"]
            )
            for s in entry["samples"]:
                assert s["ad"][0] + s["ad"][1] == s["dp"]

    def test_indel_mode_byte_identical_and_matches_expected(self, workspace):
        batch_code, batch_out, batch_err = run_batch(workspace, ["--call-indels"])
        comp_code, comp_out, comp_err = run_per_sample(workspace, ["--call-indels"])
        assert batch_code == comp_code == 0
        assert batch_err == comp_err == ""
        assert batch_out == comp_out
        entries = parse_lines(batch_out)
        assert entries == EXPECTED_INDEL
        # The homopolymer deletion is left-aligned and merged across
        # samples into a single normalized event.
        indel = entries[-1]
        assert (indel["chrom"], indel["pos"], indel["ref"], indel["alt"]) == (
            "c2",
            7,
            "CA",
            "C",
        )
        assert [s["sample"] for s in indel["samples"]] == ["s1", "s4"]

    def test_thresholds_are_threaded_through_both_paths(self, workspace):
        # A non-default threshold set must shift both paths together.
        args = ["--min-alt-count", "3", "--min-base-quality", "30"]
        batch_code, batch_out, _ = run_batch(workspace, args)
        comp_code, comp_out, _ = run_per_sample(workspace, args)
        assert batch_code == comp_code == 0
        assert batch_out == comp_out
        # c1:20 (AC=2) and the s1 c1:10 support (AC=2) drop out; the
        # low-quality base was already excluded by the default floor.
        keys = [(e["chrom"], e["pos"]) for e in parse_lines(batch_out)]
        assert keys == [("c1", 10), ("c2", 3)]

    def test_file_manifests_resolve_against_manifest_dir_from_any_cwd(
        self, workspace, monkeypatch
    ):
        monkeypatch.chdir(workspace)
        code, batch_local, _ = run(
            ["batch-call-variants", "m.jsonl", "--reference", "ref.fa"]
        )
        assert code == 0
        monkeypatch.chdir("/")
        code, batch_elsewhere, err = run_batch(workspace)
        assert code == 0, err
        assert batch_elsewhere == batch_local

        monkeypatch.chdir(workspace)
        comp_code, comp_local, comp_err = run_per_sample(workspace)
        assert comp_code == 0, comp_err
        monkeypatch.chdir("/")
        comp_code, comp_elsewhere, comp_err = run_per_sample(workspace)
        assert comp_code == 0, comp_err
        assert comp_elsewhere == comp_local == batch_local

    def test_repeated_runs_are_byte_identical(self, workspace):
        argv = [
            "batch-call-variants",
            str(workspace / "m.jsonl"),
            "--reference",
            str(workspace / "ref.fa"),
            "--call-indels",
        ]
        first = workspace / "out1.jsonl"
        second = workspace / "out2.jsonl"
        code1, out1, err1 = run(argv + ["--output", str(first)])
        code2, out2, err2 = run(argv + ["--output", str(second)])
        assert code1 == code2 == 0
        assert out1 == out2 == ""
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        code3, out3, err3 = run(argv)
        assert code3 == 0 and err3 == ""
        assert out3.encode() == first.read_bytes()


class TestEmptyAndDegenerateContracts:
    def test_empty_manifest_produces_empty_output_on_both_paths(self, workspace):
        (workspace / "m.jsonl").write_text("\n\n")
        batch_code, batch_out, batch_err = run_batch(workspace)
        assert (batch_code, batch_out, batch_err) == (0, "", "")

        (workspace / "vm.jsonl").write_text("\n")
        comp_code, comp_out, comp_err = run(
            [
                "summarize-variants",
                str(workspace / "vm.jsonl"),
                "--reference",
                str(workspace / "ref.fa"),
            ]
        )
        assert (comp_code, comp_out, comp_err) == (0, "", "")

    def test_all_samples_without_variants_produce_empty_output(self, workspace):
        (workspace / "m.jsonl").write_text(
            '{"sample":"s3","reads":"reads/s3.fq"}\n'
        )
        batch_code, batch_out, batch_err = run_batch(workspace)
        assert (batch_code, batch_out, batch_err) == (0, "", "")

        comp_code, comp_out, comp_err = run_per_sample(workspace, samples=("s3",))
        assert (comp_code, comp_out, comp_err) == (0, "", "")

    def test_trailing_corrupt_fastq_fails_atomically(self, workspace):
        (workspace / "reads" / "broken.fq").write_text("@x\nACGT\n+x\nII\n")
        (workspace / "m.jsonl").write_text(
            '{"sample":"s2","reads":"reads/s2.fq"}\n'
            '{"sample":"broken","reads":"reads/broken.fq"}\n'
        )
        code, out, err = run_batch(workspace)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

        # With an existing output target the failure preserves it and
        # leaves no partial summary behind.
        target = workspace / "out.jsonl"
        target.write_text("PREVIOUS\n")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(workspace / "ref.fa"),
                "--output",
                str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"
        assert list(workspace.glob(".batch-call-variants-*")) == []

        # The composed path rejects the same corrupt sample as well.
        code, out, err = run(
            [
                "call-variants",
                str(workspace / "ref.fa"),
                str(workspace / "reads" / "broken.fq"),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1


class TestIterationChunking:
    """Byte output must not depend on how iterables are chunked."""

    def test_reads_iterable_chunking_keeps_vcf_bytes(self, workspace):
        reference = list(read_sequences(str(workspace / "ref.fa"), format="fasta"))
        reads = list(read_sequences(str(workspace / "reads" / "s1.fq"), format="fastq"))
        expected = render_vcf(call_variants(reference, reads, sample_name="s1"))
        chunkings = [
            iter(reads),
            (read for read in reads),
            itertools.chain.from_iterable([read] for read in reads),
        ]
        for chunking in chunkings:
            assert (
                render_vcf(call_variants(reference, chunking, sample_name="s1"))
                == expected
            )

    def test_manifest_entry_chunking_keeps_summary_bytes(self, workspace):
        entries = read_batch_manifest(str(workspace / "m.jsonl"))
        reference = str(workspace / "ref.fa")
        expected = render_summaries(
            batch_call_variants(str(workspace / "m.jsonl"), reference)
        )
        for form in (entries, list(entries), iter(entries)):
            assert render_summaries(batch_call_variants(form, reference)) == expected
        # The library-level merge matches the command-line bytes.
        code, out, err = run_batch(workspace)
        assert code == 0 and err == ""
        assert out == expected


class TestErrorContractUnchanged:
    """The consistency tests above must not come from relaxed validation."""

    @pytest.mark.parametrize(
        "option,value",
        [
            ("--min-base-quality", "94"),
            ("--min-base-quality", "-1"),
            ("--min-alt-count", "0"),
            ("--min-alt-fraction", "1.5"),
            ("--homozygous-fraction", "nan"),
        ],
    )
    def test_invalid_thresholds_exit_2_in_both_calling_entries(
        self, workspace, option, value
    ):
        code, out, err = run(
            [
                "call-variants",
                str(workspace / "ref.fa"),
                str(workspace / "reads" / "s1.fq"),
                option,
                value,
            ]
        )
        assert code == 2 and out == "" and err.count("\n") == 1
        code, out, err = run_batch(workspace, [option, value])
        assert code == 2 and out == "" and err.count("\n") == 1

    def test_homozygous_below_min_alt_fraction_exit_2_in_both(self, workspace):
        args = ["--min-alt-fraction", "0.5", "--homozygous-fraction", "0.4"]
        code, out, err = run(
            [
                "call-variants",
                str(workspace / "ref.fa"),
                str(workspace / "reads" / "s1.fq"),
                *args,
            ]
        )
        assert code == 2 and out == "" and err.count("\n") == 1
        code, out, err = run_batch(workspace, args)
        assert code == 2 and out == "" and err.count("\n") == 1

    def test_duplicate_sample_exit_2_in_both_manifest_kinds(self, workspace):
        (workspace / "m.jsonl").write_text(
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
            '{"sample":"s1","reads":"reads/s2.fq"}\n'
        )
        code, out, err = run_batch(workspace)
        assert code == 2 and out == "" and err.count("\n") == 1
        assert "duplicate" in err

        (workspace / "vm.jsonl").write_text(
            '{"sample":"s1","vcf":"a.vcf"}\n{"sample":"s1","vcf":"b.vcf"}\n'
        )
        code, out, err = run(
            [
                "summarize-variants",
                str(workspace / "vm.jsonl"),
                "--reference",
                str(workspace / "ref.fa"),
            ]
        )
        assert code == 2 and out == "" and err.count("\n") == 1
        assert "duplicate" in err

    def test_reference_errors_keep_their_exit_codes(self, workspace):
        empty = workspace / "empty.fa"
        empty.write_text("")
        duplicate = workspace / "dup.fa"
        duplicate.write_text(">c1\nACGT\n>c1\nTTTT\n")
        for reference in (empty, duplicate):
            code, out, err = run(
                [
                    "call-variants",
                    str(reference),
                    str(workspace / "reads" / "s1.fq"),
                ]
            )
            assert code == 2 and out == "" and err.count("\n") == 1
            code, out, err = run(
                [
                    "batch-call-variants",
                    str(workspace / "m.jsonl"),
                    "--reference",
                    str(reference),
                ]
            )
            assert code == 2 and out == "" and err.count("\n") == 1

    def test_file_access_errors_keep_exit_1(self, workspace):
        # Missing manifest.
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "gone.jsonl"),
                "--reference",
                str(workspace / "ref.fa"),
            ]
        )
        assert code == 1 and out == "" and err.count("\n") == 1
        # Missing reads file referenced by the manifest.
        (workspace / "m.jsonl").write_text(
            '{"sample":"s1","reads":"reads/gone.fq"}\n'
        )
        code, out, err = run_batch(workspace)
        assert code == 1 and out == "" and err.count("\n") == 1
        # Missing reference for either calling entry.
        for argv in (
            [
                "call-variants",
                str(workspace / "gone.fa"),
                str(workspace / "reads" / "s1.fq"),
            ],
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(workspace / "gone.fa"),
            ],
        ):
            code, out, err = run(argv)
            assert code == 1 and out == "" and err.count("\n") == 1
        # Missing VCF for the summarize entry.
        (workspace / "vm.jsonl").write_text('{"sample":"s1","vcf":"gone.vcf"}\n')
        code, out, err = run(
            [
                "summarize-variants",
                str(workspace / "vm.jsonl"),
                "--reference",
                str(workspace / "ref.fa"),
            ]
        )
        assert code == 1 and out == "" and err.count("\n") == 1
