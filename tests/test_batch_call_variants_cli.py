"""Tests for the ``batch-call-variants`` command line entry point."""

from __future__ import annotations

import json
import sys
from io import StringIO

import pytest

from genome_variant.cli import main

REFERENCE = ">c1\nACGTACGTACAA\n"
ALT_READ = "ACGTTCGTACAA"
REF_READ = "ACGTACGTACAA"

# A non-repetitive reference; inserting "TT" after the "C" at POS 15 is a
# call-indels event (see the call-variants indel tests).
INDEL_REF_SEQ = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"
INDEL_REFERENCE = f">c1\n{INDEL_REF_SEQ}\n"
INS_READ = INDEL_REF_SEQ[:15] + "TT" + INDEL_REF_SEQ[15:]


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


def fastq(identifier, sequence, quality=None):
    if quality is None:
        quality = "I" * len(sequence)
    return f"@{identifier}\n{sequence}\n+{identifier}\n{quality}\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.fixture
def reference_file(tmp_path):
    return write(tmp_path, "ref.fa", REFERENCE)


@pytest.fixture
def workspace(tmp_path, reference_file):
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "s1.fq").write_text(
        fastq("a", ALT_READ) + fastq("b", ALT_READ) + fastq("c", REF_READ)
    )
    (reads / "s2.fq").write_text(fastq("d", ALT_READ) + fastq("e", ALT_READ))
    manifest = tmp_path / "m.jsonl"
    manifest.write_text(
        '{"sample":"s1","reads":"reads/s1.fq"}\n'
        '{"sample":"s2","reads":"reads/s2.fq"}\n'
    )
    return tmp_path


class TestBatchCallVariantsCli:
    def test_basic_merge_stdout(self, workspace, reference_file) -> None:
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 0
        assert err == ""
        data = [json.loads(line) for line in out.splitlines()]
        assert data == [
            {
                "chrom": "c1",
                "pos": 5,
                "ref": "A",
                "alt": "T",
                "sample_count": 2,
                "allele_count": 3,
                "depth": 5,
                "samples": [
                    {"sample": "s1", "gt": "0/1", "dp": 3, "ad": [1, 2]},
                    {"sample": "s2", "gt": "1/1", "dp": 2, "ad": [0, 2]},
                ],
            }
        ]
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_samples_keep_manifest_order(self, workspace, reference_file) -> None:
        (workspace / "m.jsonl").write_text(
            '{"sample":"s2","reads":"reads/s2.fq"}\n'
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
        )
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 0
        merged = json.loads(out)
        assert [s["sample"] for s in merged["samples"]] == ["s2", "s1"]

    def test_empty_manifest_produces_empty_output(
        self, workspace, reference_file
    ) -> None:
        (workspace / "m.jsonl").write_text("\n\n")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_no_variants_produces_empty_output(
        self, workspace, reference_file
    ) -> None:
        (workspace / "reads" / "w.fq").write_text(fastq("w", REF_READ))
        (workspace / "m.jsonl").write_text('{"sample":"w","reads":"reads/w.fq"}\n')
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_file_manifest_resolves_vs_its_directory(
        self, workspace, reference_file, monkeypatch
    ) -> None:
        monkeypatch.chdir("/")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "c1"

    def test_stdin_manifest_resolves_vs_cwd(
        self, workspace, reference_file, monkeypatch
    ) -> None:
        monkeypatch.chdir(workspace)
        code, out, err = run(
            ["batch-call-variants", "-", "--reference", str(reference_file)],
            '{"sample":"s1","reads":"reads/s1.fq"}\n',
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "c1"

    def test_reference_can_be_stdin_with_file_manifest(
        self, workspace, reference_file
    ) -> None:
        code, out, err = run(
            ["batch-call-variants", str(workspace / "m.jsonl"), "--reference", "-"],
            REFERENCE,
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "c1"

    def test_manifest_and_reference_cannot_both_be_stdin(self) -> None:
        code, out, err = run(
            ["batch-call-variants", "-", "--reference", "-"],
            '{"sample":"s1","reads":"x.fq"}\n',
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_threshold_options_match_call_variants(
        self, workspace, reference_file
    ) -> None:
        (workspace / "m.jsonl").write_text('{"sample":"s1","reads":"reads/s1.fq"}\n')
        argv = [
            "batch-call-variants",
            str(workspace / "m.jsonl"),
            "--reference",
            str(reference_file),
        ]
        # AC=2 of DP=3: dropped by --min-alt-count 3.
        code, out, err = run(argv + ["--min-alt-count", "3"])
        assert code == 0
        assert out == ""
        # Dropped by --min-alt-fraction 0.9.
        code, out, err = run(
            argv + ["--min-alt-fraction", "0.9", "--homozygous-fraction", "0.9"]
        )
        assert code == 0
        assert out == ""
        # Genotyped 1/1 by --homozygous-fraction 0.5.
        code, out, err = run(argv + ["--homozygous-fraction", "0.5"])
        assert code == 0
        assert json.loads(out)["samples"][0]["gt"] == "1/1"
        # Bases below --min-base-quality contribute nothing.
        code, out, err = run(argv + ["--min-base-quality", "41"])
        assert code == 0
        assert out == ""

    def test_call_indels(self, tmp_path) -> None:
        reference = write(tmp_path, "ref.fa", INDEL_REFERENCE)
        inserted = INS_READ
        wild = INDEL_REF_SEQ
        reads = tmp_path / "reads"
        reads.mkdir()
        (reads / "s1.fq").write_text(
            fastq("a", inserted)
            + fastq("b", inserted)
            + fastq("c", inserted)
            + fastq("w", wild)
        )
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"s1","reads":"reads/s1.fq"}\n'
        )
        argv = [
            "batch-call-variants",
            str(manifest),
            "--reference",
            str(reference),
        ]
        code, out, err = run(argv)
        assert code == 0
        assert out == ""  # SNV-only default: the insertion is not called
        code, out, err = run(argv + ["--call-indels"])
        assert code == 0
        merged = json.loads(out)
        assert (merged["chrom"], merged["pos"], merged["ref"], merged["alt"]) == (
            "c1",
            15,
            "C",
            "CTT",
        )
        assert merged["samples"] == [
            {"sample": "s1", "gt": "0/1", "dp": 4, "ad": [1, 3]}
        ]

    def test_max_indel_length_requires_call_indels(
        self, workspace, reference_file
    ) -> None:
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
                "--max-indel-length",
                "10",
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_homozygous_fraction_below_min_alt_fraction(
        self, workspace, reference_file
    ) -> None:
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
                "--min-alt-fraction",
                "0.5",
                "--homozygous-fraction",
                "0.4",
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

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
    def test_bad_thresholds_exit_2(
        self, workspace, reference_file, option, value
    ) -> None:
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
                option,
                value,
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    @pytest.mark.parametrize(
        "line",
        [
            "{not json}",
            '["s1","reads/s1.fq"]',
            '{"sample":"s1"}',
            '{"sample":"s1","reads":"reads/s1.fq","vcf":"x"}',
            '{"sample":"","reads":"reads/s1.fq"}',
            '{"sample":"  ","reads":"reads/s1.fq"}',
            '{"sample":"a\\tb","reads":"reads/s1.fq"}',
            '{"sample":3,"reads":"reads/s1.fq"}',
            '{"sample":"s1","reads":""}',
            '{"sample":"s1","reads":7}',
            '{"sample":"s1","reads":"-"}',
        ],
    )
    def test_bad_manifest_lines_exit_2(
        self, workspace, reference_file, line
    ) -> None:
        (workspace / "m.jsonl").write_text(line + "\n")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":1:" in err

    def test_duplicate_sample_reports_second_line(
        self, workspace, reference_file
    ) -> None:
        (workspace / "m.jsonl").write_text(
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
            '{"sample":"s1","reads":"reads/s2.fq"}\n'
        )
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":2:" in err
        assert "duplicate" in err

    def test_missing_manifest_file_exit_1(self, workspace, reference_file) -> None:
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "gone.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reads_file_exit_1(self, workspace, reference_file) -> None:
        (workspace / "m.jsonl").write_text('{"sample":"s1","reads":"reads/gone.fq"}\n')
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_malformed_fastq_exit_2(self, workspace, reference_file) -> None:
        (workspace / "reads" / "bad.fq").write_text("@x\nACGT\n+x\nII\n")
        (workspace / "m.jsonl").write_text('{"sample":"s1","reads":"reads/bad.fq"}\n')
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2(self, workspace, reference_file) -> None:
        (workspace / "reads" / "bad.fq").write_text(fastq("x", "ACGTZCGTACAA"))
        (workspace / "m.jsonl").write_text('{"sample":"s1","reads":"reads/bad.fq"}\n')
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, workspace) -> None:
        reference = write(workspace, "empty.fa", "")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_duplicate_reference_identifier_exit_2(self, workspace) -> None:
        reference = write(workspace, "dup.fa", ">c1\nACGT\n>c1\nTTTT\n")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_output_file_roundtrip_and_byte_stability(
        self, workspace, reference_file
    ) -> None:
        first = workspace / "out1.jsonl"
        second = workspace / "out2.jsonl"
        argv = [
            "batch-call-variants",
            str(workspace / "m.jsonl"),
            "--reference",
            str(reference_file),
        ]
        code1, out1, err1 = run(argv + ["--output", str(first)])
        code2, out2, err2 = run(argv + ["--output", str(second)])
        assert code1 == code2 == 0
        assert out1 == out2 == ""
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        assert first.read_bytes().endswith(b"\n")
        assert not first.read_bytes().endswith(b"\n\n")

    def test_empty_result_writes_empty_file(self, workspace, reference_file) -> None:
        (workspace / "m.jsonl").write_text("")
        target = workspace / "out.jsonl"
        target.write_text("OLD\n")
        code, _, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
                "--output",
                str(target),
            ]
        )
        assert code == 0
        assert err == ""
        assert target.read_bytes() == b""

    def test_failure_preserves_existing_target_and_leaves_no_temp(
        self, workspace, reference_file
    ) -> None:
        (workspace / "m.jsonl").write_text('{"sample":"s1","reads":"-"}\n')
        target = workspace / "out.jsonl"
        target.write_text("PREVIOUS\n")
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
                "--output",
                str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"
        assert list(workspace.glob(".batch-call-variants-*")) == []

    def test_output_write_failure_exit_1(self, workspace, reference_file) -> None:
        code, out, err = run(
            [
                "batch-call-variants",
                str(workspace / "m.jsonl"),
                "--reference",
                str(reference_file),
                "--output",
                str(workspace / "no-such-dir" / "out.jsonl"),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1
