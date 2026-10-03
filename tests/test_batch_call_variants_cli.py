"""Tests for the ``batch-call-variants`` command line entry point."""

from __future__ import annotations

import json
import sys
from io import StringIO

import pytest

from genome_variant.cli import main

REF = "ACGTACGTACAA"
REF_FASTA = f">c1\n{REF}\n"


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


def alt_read(identifier, base="T", position=4):
    sequence = list(REF)
    sequence[position] = base
    return fastq(identifier, "".join(sequence))


def het_reads():
    return alt_read("a") + alt_read("b") + fastq("c", REF)


def hom_reads():
    return alt_read("a") + alt_read("b") + alt_read("c")


@pytest.fixture
def workspace(tmp_path):
    (tmp_path / "ref.fa").write_text(REF_FASTA)
    (tmp_path / "reads").mkdir()
    return tmp_path


def write_reads(workspace, name, text):
    path = workspace / "reads" / name
    path.write_text(text)
    return path


def write_manifest(workspace, lines):
    manifest = workspace / "m.jsonl"
    manifest.write_text("".join(line + "\n" for line in lines))
    return manifest


class TestBatchCallVariantsCli:
    def test_two_sample_merge(self, workspace) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        write_reads(workspace, "s2.fq", hom_reads())
        manifest = write_manifest(
            workspace,
            [
                '{"sample":"s1","reads":"reads/s1.fq"}',
                '{"sample":"s2","reads":"reads/s2.fq"}',
            ],
        )
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
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
                "depth": 6,
                "samples": [
                    {"sample": "s1", "gt": "0/1", "dp": 3, "ad": [1, 2]},
                    {"sample": "s2", "gt": "1/1", "dp": 3, "ad": [0, 3]},
                ],
            }
        ]
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_samples_keep_manifest_order(self, workspace) -> None:
        write_reads(workspace, "a.fq", het_reads())
        write_reads(workspace, "b.fq", het_reads())
        manifest = write_manifest(
            workspace,
            [
                '{"sample":"b","reads":"reads/b.fq"}',
                '{"sample":"a","reads":"reads/a.fq"}',
            ],
        )
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        merged = json.loads(out)
        assert [s["sample"] for s in merged["samples"]] == ["b", "a"]

    def test_output_sorted_by_reference_then_pos(self, workspace) -> None:
        # Non-repetitive reference so full-length alignments always win.
        c1 = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"
        (workspace / "ref.fa").write_text(f">c2\nTTTTTT\n>c1\n{c1}\n")

        def alt(identifier, position, base):
            sequence = list(c1)
            sequence[position] = base
            return fastq(identifier, "".join(sequence))

        reads = (
            alt("a", 5, "T") + alt("b", 5, "T")      # C->T at POS 6
            + alt("c", 20, "A") + alt("d", 20, "A")  # C->A at POS 21
        )
        write_reads(workspace, "s1.fq", reads)
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        keys = [
            (d["chrom"], d["pos"])
            for d in (json.loads(line) for line in out.splitlines())
        ]
        assert keys == [("c1", 6), ("c1", 21)]

    def test_blank_manifest_lines_ignored(self, workspace) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        manifest = workspace / "m.jsonl"
        manifest.write_text(
            "\n"
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
            "   \n"
        )
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        assert out.count("\n") == 1

    def test_empty_manifest_produces_empty_output(self, workspace) -> None:
        manifest = write_manifest(workspace, [])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_no_variants_produces_empty_output(self, workspace) -> None:
        write_reads(workspace, "s1.fq", fastq("a", REF))
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_stdin_manifest_resolves_vs_cwd(self, workspace, monkeypatch) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        monkeypatch.chdir(workspace)
        code, out, err = run(
            ["batch-call-variants", "-", "--reference", str(workspace / "ref.fa")],
            '{"sample":"s1","reads":"reads/s1.fq"}\n',
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "c1"

    def test_file_manifest_resolves_vs_its_directory(
        self, workspace, monkeypatch
    ) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        # Running from an unrelated cwd must still find the reads.
        monkeypatch.chdir("/")
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "c1"

    def test_reference_can_be_stdin_with_file_manifest(self, workspace) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", "-"],
            REF_FASTA,
        )
        assert code == 0
        assert json.loads(out)["chrom"] == "c1"

    def test_manifest_and_reference_cannot_both_be_stdin(self) -> None:
        code, out, err = run(
            ["batch-call-variants", "-", "--reference", "-"],
            '{"sample":"s1","reads":"x"}\n',
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "standard input" in err

    def test_non_ascii_sample_preserved_verbatim(self, workspace) -> None:
        write_reads(workspace, "u.fq", het_reads())
        manifest = workspace / "m.jsonl"
        manifest.write_text('{"sample":"样","reads":"reads/u.fq"}\n', encoding="utf-8")
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 0
        assert "样" in out
        assert "\\u" not in out

    def test_threshold_options_match_call_variants(self, workspace) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}']
        )
        argv = [
            "batch-call-variants", str(manifest),
            "--reference", str(workspace / "ref.fa"),
        ]
        code, out, _ = run(argv + ["--min-alt-count", "3"])
        assert code == 0
        assert out == ""
        code, out, _ = run(argv + ["--homozygous-fraction", "0.6"])
        assert code == 0
        assert json.loads(out)["samples"][0]["gt"] == "1/1"

    def test_call_indels(self, workspace) -> None:
        indel_ref = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"
        inserted = indel_ref[:15] + "TT" + indel_ref[15:]
        (workspace / "ref.fa").write_text(f">c1\n{indel_ref}\n")
        write_reads(
            workspace,
            "s1.fq",
            fastq("a", inserted) + fastq("b", inserted) + fastq("w", indel_ref),
        )
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        argv = [
            "batch-call-variants", str(manifest),
            "--reference", str(workspace / "ref.fa"),
        ]
        code, out, _ = run(argv)
        assert code == 0
        assert out == ""
        code, out, err = run(argv + ["--call-indels", "--max-indel-length", "10"])
        assert code == 0
        merged = json.loads(out)
        assert (merged["pos"], merged["ref"], merged["alt"]) == (15, "C", "CTT")

    def test_max_indel_length_requires_call_indels(self, workspace) -> None:
        manifest = write_manifest(workspace, [])
        code, out, err = run(
            [
                "batch-call-variants", str(manifest),
                "--reference", str(workspace / "ref.fa"),
                "--max-indel-length", "10",
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_homozygous_fraction_below_min_alt_fraction(self, workspace) -> None:
        manifest = write_manifest(workspace, [])
        code, out, err = run(
            [
                "batch-call-variants", str(manifest),
                "--reference", str(workspace / "ref.fa"),
                "--min-alt-fraction", "0.9",
                "--homozygous-fraction", "0.8",
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_output_file_roundtrip_and_single_newline(self, workspace) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        target = workspace / "out.jsonl"
        code, out, err = run(
            [
                "batch-call-variants", str(manifest),
                "--reference", str(workspace / "ref.fa"),
                "--output", str(target),
            ]
        )
        assert code == 0
        assert out == "" and err == ""
        data = target.read_bytes()
        assert data.endswith(b"\n") and not data.endswith(b"\n\n")

    def test_empty_result_writes_empty_file(self, workspace) -> None:
        manifest = write_manifest(workspace, [])
        target = workspace / "out.jsonl"
        target.write_text("OLD\n")
        code, _, err = run(
            [
                "batch-call-variants", str(manifest),
                "--reference", str(workspace / "ref.fa"),
                "--output", str(target),
            ]
        )
        assert code == 0
        assert err == ""
        assert target.read_bytes() == b""

    def test_byte_stable_across_runs(self, workspace) -> None:
        write_reads(workspace, "s1.fq", het_reads())
        write_reads(workspace, "s2.fq", hom_reads())
        manifest = write_manifest(
            workspace,
            [
                '{"sample":"s1","reads":"reads/s1.fq"}',
                '{"sample":"s2","reads":"reads/s2.fq"}',
            ],
        )
        argv = [
            "batch-call-variants", str(manifest),
            "--reference", str(workspace / "ref.fa"),
        ]
        _, out1, err1 = run(argv)
        _, out2, err2 = run(argv)
        assert err1 == err2 == ""
        assert out1 == out2

    # -- manifest validation ------------------------------------------------

    def test_bad_json_reports_line_number(self, workspace) -> None:
        code, out, err = run(
            ["batch-call-variants", "-", "--reference", str(workspace / "ref.fa")],
            '{"sample":"a","reads":"x"}\n{not json}\n',
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":2:" in err

    def test_line_that_is_not_object(self, workspace) -> None:
        code, _, err = run(
            ["batch-call-variants", "-", "--reference", str(workspace / "ref.fa")],
            '["sample","reads"]\n',
        )
        assert code == 2
        assert ":1:" in err

    @pytest.mark.parametrize(
        "text",
        [
            '{"sample":"a"}\n',                          # missing reads
            '{"reads":"x"}\n',                           # missing sample
            '{"sample":"a","reads":"x","extra":1}\n',    # extra field
            '{"sample":"","reads":"x"}\n',               # empty sample
            '{"sample":"  ","reads":"x"}\n',             # whitespace-only sample
            '{"sample":"a\tb","reads":"x"}\n',           # tab in sample
            '{"sample":3,"reads":"x"}\n',                # non-string sample
            '{"sample":"a","reads":""}\n',               # empty reads
            '{"sample":"a","reads":4}\n',                # non-string reads
            '{"sample":"a","reads":"-"}\n',              # standard input reads
        ],
    )
    def test_bad_manifest_fields_exit_2(self, text, workspace) -> None:
        code, out, err = run(
            ["batch-call-variants", "-", "--reference", str(workspace / "ref.fa")],
            text,
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert ":1:" in err

    def test_duplicate_sample_reports_second_line(self, workspace) -> None:
        code, _, err = run(
            ["batch-call-variants", "-", "--reference", str(workspace / "ref.fa")],
            '{"sample":"a","reads":"x"}\n{"sample":"a","reads":"y"}\n',
        )
        assert code == 2
        assert ":2:" in err
        assert "duplicate" in err

    # -- sequence / quality / reference data errors --------------------------

    def test_fastq_format_error_exit_2(self, workspace) -> None:
        write_reads(workspace, "s1.fq", "@a\nACGT\n+a\n")  # truncated record
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2(self, workspace) -> None:
        write_reads(workspace, "s1.fq", fastq("a", "ACGTACGTACA!"))
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, workspace) -> None:
        (workspace / "ref.fa").write_text("")
        write_reads(workspace, "s1.fq", het_reads())
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_duplicate_reference_identifier_exit_2(self, workspace) -> None:
        (workspace / "ref.fa").write_text(">c1\nACGT\n>c1\nTTTT\n")
        write_reads(workspace, "s1.fq", het_reads())
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/s1.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 2
        assert out == ""
        assert "duplicate" in err

    # -- file handling -------------------------------------------------------

    def test_missing_manifest_file_exit_1(self, workspace) -> None:
        code, out, err = run(
            [
                "batch-call-variants", str(workspace / "nope.jsonl"),
                "--reference", str(workspace / "ref.fa"),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reads_file_exit_1(self, workspace) -> None:
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/gone.fq"}'])
        code, out, err = run(
            ["batch-call-variants", str(manifest), "--reference", str(workspace / "ref.fa")]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reference_file_exit_1(self, workspace) -> None:
        manifest = write_manifest(workspace, ['{"sample":"s1","reads":"reads/x.fq"}'])
        code, out, err = run(
            [
                "batch-call-variants", str(manifest),
                "--reference", str(workspace / "gone.fa"),
            ]
        )
        assert code == 1
        assert out == ""

    def test_failure_preserves_existing_output(self, workspace) -> None:
        write_reads(workspace, "good.fq", het_reads())
        write_reads(workspace, "bad.fq", "@a\nACGT\n+a\n")  # truncated record
        manifest = write_manifest(
            workspace,
            [
                '{"sample":"s1","reads":"reads/good.fq"}',
                '{"sample":"s2","reads":"reads/bad.fq"}',
            ],
        )
        target = workspace / "out.jsonl"
        target.write_text("PREVIOUS\n")
        code, out, err = run(
            [
                "batch-call-variants", str(manifest),
                "--reference", str(workspace / "ref.fa"),
                "--output", str(target),
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS\n"
        # No temporary files left behind.
        assert [p.name for p in workspace.iterdir() if p.name.startswith(".batch-")] == []

    def test_reference_required(self) -> None:
        code, _, err = run(["batch-call-variants", "-"], "")
        assert code == 2
        assert err.count("\n") == 1
