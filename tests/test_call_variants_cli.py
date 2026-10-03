"""Tests for the call-variants subcommand."""

from __future__ import annotations

import sys


def run(argv, stdin: str = ""):
    from io import StringIO

    old_in, old_out, old_err = sys.stdin, sys.stdout, sys.stderr
    sys.stdin = StringIO(stdin)
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    try:
        try:
            from genome_variant.cli import main

            code = main(argv)
        except SystemExit as exc:
            code = exc.code
    finally:
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdin, sys.stdout, sys.stderr = old_in, old_out, old_err
    return code, out, err


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


# Reference and reads with an interior A->T variant shared by two reads.
REFERENCE = ">c1\nACGTACGTACAA\n"
ALT_READ = "ACGTTCGTACAA"
REF_READ = "ACGTACGTACAA"


def fastq(identifier, sequence, quality=None):
    if quality is None:
        quality = "I" * len(sequence)
    return f"@{identifier}\n{sequence}\n+{identifier}\n{quality}\n"


EXPECTED_RECORD = "c1\t5\t.\tA\tT\t.\tPASS\tDP=3;AC=2;AF=0.666667\tGT:DP:AD\t0/1:3:1,2"


class TestCallVariants:
    def test_basic_vcf_output(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(
            tmp_path,
            "q.fq",
            fastq("a", ALT_READ) + fastq("b", ALT_READ) + fastq("c", REF_READ),
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert lines[0] == "##fileformat=VCFv4.2"
        assert lines[-2] == (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        )
        assert lines[-1] == EXPECTED_RECORD
        assert out.endswith("\n")
        assert not out.endswith("\n\n")

    def test_output_is_readable_as_vcf(self, tmp_path) -> None:
        from io import StringIO

        from genome_variant.vcf import read_vcf

        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(
            tmp_path, "q.fq", fastq("a", ALT_READ) + fastq("b", ALT_READ)
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        document = read_vcf(StringIO(out))
        assert document.header.samples == ("SAMPLE",)
        assert len(document.records) == 1
        assert document.records[0].sample_text == ("1/1:2:0,2",)

    def test_no_candidates_still_writes_full_header(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGT\n")
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert lines[0] == "##fileformat=VCFv4.2"
        assert lines[-1].startswith("#CHROM")
        assert out.count("\n") == 8

    def test_multi_reference_ordering(self, tmp_path) -> None:
        reference = write(
            tmp_path,
            "r.fa",
            ">c1\nACGTACGTACAA\n>c2\nTAACGTTTCC\n",
        )
        # c2 read ACGATA maps reverse-strand with A->T at POS 3.
        reads = write(
            tmp_path,
            "q.fq",
            fastq("a", ALT_READ) + fastq("b", "ACGATA", "I" * 6),
        )
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--min-alt-count", "1"]
        )
        assert code == 0
        rows = [line for line in out.splitlines() if not line.startswith("#")]
        assert [(line.split("\t")[0], line.split("\t")[1]) for line in rows] == [
            ("c1", "5"),
            ("c2", "3"),
        ]

    def test_reads_from_stdin(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        code, out, err = run(
            ["call-variants", str(reference), "-"],
            fastq("a", ALT_READ) + fastq("b", ALT_READ),
        )
        assert code == 0
        assert "A\tT" in out

    def test_reference_from_stdin(self, tmp_path) -> None:
        reads = write(
            tmp_path, "q.fq", fastq("a", ALT_READ) + fastq("b", ALT_READ)
        )
        code, out, err = run(
            ["call-variants", "-", str(reads)], REFERENCE
        )
        assert code == 0
        assert "c1\t5" in out

    def test_custom_sample_name(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", fastq("a", ALT_READ) + fastq("b", ALT_READ))
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--sample-name",
                "NA12878",
            ]
        )
        assert code == 0
        assert out.splitlines()[-2].endswith("\tNA12878")

    def test_threshold_options_take_effect(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        # One ALT observation, one reference observation: count below 2.
        reads = write(
            tmp_path, "q.fq", fastq("a", ALT_READ) + fastq("c", REF_READ)
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert not [line for line in out.splitlines() if not line.startswith("#")]

        # Lowering the count admits it; fraction 0.5 calls 0/1.
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--min-alt-count",
                "1",
                "--min-alt-fraction",
                "0.5",
                "--homozygous-fraction",
                "0.9",
            ]
        )
        assert code == 0
        record = [line for line in out.splitlines() if not line.startswith("#")][0]
        assert record.endswith("0/1:2:1,1")

    def test_min_basequality_filters_evidence(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        # Variant base is at sequence index 4; give both ALT reads
        # Phred 19 there, below the default 20.
        quality = ["I"] * len(ALT_READ)
        quality[4] = "4"  # ASCII 52 - 33 = 19
        qtext = "".join(quality)
        reads = write(
            tmp_path,
            "q.fq",
            fastq("a", ALT_READ, qtext) + fastq("b", ALT_READ, qtext),
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert not [line for line in out.splitlines() if not line.startswith("#")]

        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--min-base-quality", "19"]
        )
        assert code == 0
        assert any(line.startswith("c1\t5") for line in out.splitlines())

    def test_fasta_reads_fail_with_quality_error(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fa", f">a\n{ALT_READ}\n")
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--reads-format", "fasta"]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "quality" in err


class TestCallVariantsErrors:
    def test_both_inputs_stdin_rejected(self) -> None:
        code, out, err = run(["call-variants", "-", "-"], REFERENCE)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, tmp_path) -> None:
        empty = write(tmp_path, "empty.fa", "")
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        code, out, err = run(["call-variants", str(empty), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_duplicate_reference_identifier_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">x\nACGT\n>x\nTTTT\n")
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_sample_name_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        for name in ("", "   ", "a\tb"):
            code, out, err = run(
                [
                    "call-variants",
                    str(reference),
                    str(reads),
                    "--sample-name",
                    name,
                ]
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1

    def test_bad_threshold_arguments_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        cases = [
            ["--min-base-quality", "-1"],
            ["--min-base-quality", "94"],
            ["--min-base-quality", "x"],
            ["--min-alt-count", "0"],
            ["--min-alt-count", "1.5"],
            ["--min-alt-fraction", "-0.1"],
            ["--min-alt-fraction", "1.1"],
            ["--min-alt-fraction", "nan"],
            ["--homozygous-fraction", "2"],
            ["--homozygous-fraction", "x"],
        ]
        for options in cases:
            code, out, err = run(
                ["call-variants", str(reference), str(reads), *options]
            )
            assert code == 2, options
            assert out == "", options
            assert err.count("\n") == 1, options

    def test_homozygous_below_alt_fraction_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--min-alt-fraction",
                "0.8",
                "--homozygous-fraction",
                "0.5",
            ]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_quality_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", f"@a\n{ALT_READ}\n+a\nIIIIIII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_malformed_fastq_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", "no header\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", "@a\nAZ\n+a\nII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert "'Z'" in err

    def test_missing_input_file_exit_1(self, tmp_path) -> None:
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        code, out, err = run(
            ["call-variants", str(tmp_path / "missing.fa"), str(reads)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_reads_file_exit_1(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        code, out, err = run(
            ["call-variants", str(reference), str(tmp_path / "missing.fq")]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1


class TestCallVariantsOutput:
    def test_file_output_atomic_and_byte_stable(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(
            tmp_path,
            "q.fq",
            fastq("a", ALT_READ) + fastq("b", ALT_READ) + fastq("c", REF_READ),
        )
        first = tmp_path / "a.vcf"
        second = tmp_path / "b.vcf"
        code1, _, err1 = run(
            ["call-variants", str(reference), str(reads), "--output", str(first)]
        )
        code2, _, err2 = run(
            ["call-variants", str(reference), str(reads), "--output", str(second)]
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        assert first.read_bytes().endswith(b"\n")
        assert not list(tmp_path.glob("*.tmp"))

    def test_failed_run_preserves_existing_output(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", "@a\nAZ\n+a\nII\n")
        target = tmp_path / "out.vcf"
        target.write_text("KEEP\n")
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--output", str(target)]
        )
        assert code == 2
        assert out == ""
        assert target.read_text() == "KEEP\n"
        assert {p.name for p in tmp_path.iterdir()} == {
            "r.fa",
            "q.fq",
            "out.vcf",
        }

    def test_unwritable_output_exit_1(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REFERENCE)
        reads = write(tmp_path, "q.fq", fastq("a", "ACGTACGT"))
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--output",
                str(tmp_path / "no" / "out.vcf"),
            ]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_subcommand_listed_in_help(self) -> None:
        code, out, err = run([])
        assert code == 0
        assert "call-variants" in out


# Non-repetitive reference; INS_READ carries a "TT" insertion after POS 15.
INDEL_REFERENCE = ">c1\nACGATCGTACGGATCCGTAGCTAACCGGTTAC\n"
INDEL_REF_READ = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"
INDEL_INS_READ = "ACGATCGTACGGATCTTCGTAGCTAACCGGTTAC"
EXPECTED_INS_RECORD = (
    "c1\t15\t.\tC\tCTT\t.\tPASS\tDP=4;AC=3;AF=0.750000\tGT:DP:AD\t0/1:4:1,3"
)


def indel_reads(tmp_path):
    reference = write(tmp_path, "r.fa", INDEL_REFERENCE)
    reads = write(
        tmp_path,
        "q.fq",
        fastq("a", INDEL_INS_READ)
        + fastq("b", INDEL_INS_READ)
        + fastq("c", INDEL_INS_READ)
        + fastq("w", INDEL_REF_READ),
    )
    return reference, reads


class TestCallVariantsIndels:
    def test_indels_off_by_default(self, tmp_path) -> None:
        reference, reads = indel_reads(tmp_path)
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        assert [line for line in out.splitlines() if not line.startswith("#")] == []

    def test_call_indels_reports_insertion(self, tmp_path) -> None:
        reference, reads = indel_reads(tmp_path)
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--call-indels"]
        )
        assert code == 0
        assert err == ""
        rows = [line for line in out.splitlines() if not line.startswith("#")]
        assert rows == [EXPECTED_INS_RECORD]

    def test_max_indel_length_limits_events(self, tmp_path) -> None:
        reference, reads = indel_reads(tmp_path)
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--call-indels",
                "--max-indel-length",
                "1",
            ]
        )
        assert code == 0
        assert [line for line in out.splitlines() if not line.startswith("#")] == []
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--call-indels",
                "--max-indel-length",
                "2",
            ]
        )
        assert code == 0
        rows = [line for line in out.splitlines() if not line.startswith("#")]
        assert rows == [EXPECTED_INS_RECORD]

    def test_max_indel_length_requires_call_indels(self, tmp_path) -> None:
        reference, reads = indel_reads(tmp_path)
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--max-indel-length", "30"]
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_max_indel_length_exit_2(self, tmp_path) -> None:
        reference, reads = indel_reads(tmp_path)
        for value in ("0", "-1", "1.5", "x"):
            code, out, err = run(
                [
                    "call-variants",
                    str(reference),
                    str(reads),
                    "--call-indels",
                    "--max-indel-length",
                    value,
                ]
            )
            assert code == 2, value
            assert out == "", value
            assert err.count("\n") == 1, value

    def test_indel_output_round_trips_through_normalize_vcf(self, tmp_path) -> None:
        from genome_variant.sequence_io import read_sequences
        from genome_variant.vcf import normalize_vcf, read_vcf, render_vcf
        from io import StringIO

        reference, reads = indel_reads(tmp_path)
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--call-indels"]
        )
        assert code == 0
        with open(reference) as stream:
            references = list(read_sequences(stream, format="fasta"))
        document = read_vcf(StringIO(out))
        assert len(document.records) == 1
        assert render_vcf(normalize_vcf(document, references)) == out
