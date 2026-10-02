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


HEADER_PREFIX = "##fileformat=VCFv4.2"


def parse(out):
    lines = out.splitlines()
    records = [line for line in lines if line and not line.startswith("#")]
    header = [line for line in lines if line.startswith("#")]
    return header, records


class TestCallVariants:
    def test_basic_vcf_output(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTTCGTAAAA\n+a\n" + chr(40 + 33) * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + chr(40 + 33) * 12 + "\n",
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        header, records = parse(out)
        assert header[0] == HEADER_PREFIX
        assert records == [
            "c1\t5\t.\tA\tT\t.\tPASS\tDP=2;AC=2;AF=1.000000\t"
            "GT:DP:AD\t1/1:2:0,2"
        ]
        assert header[-1] == (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        )

    def test_fasta_reads_rejected_with_quality_error(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(tmp_path, "q.fa", ">a\nACGTTCGTAAAA\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "quality" in err

    def test_reads_from_stdin(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        fastq = (
            "@a\nACGTTCGTAAAA\n+a\n" + chr(73) * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + chr(73) * 12 + "\n"
        )
        code, out, err = run(["call-variants", str(reference), "-"], fastq)
        assert code == 0
        assert err == ""
        assert "A\tT" in out

    def test_reference_from_stdin(self, tmp_path) -> None:
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTTCGTAAAA\n+a\n" + chr(73) * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + chr(73) * 12 + "\n",
        )
        code, out, err = run(["call-variants", "-", str(reads)], ">c1\nACGTACGTAAAA\n")
        assert code == 0
        assert "c1\t5" in out

    def test_no_candidates_writes_full_header(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTACGTAAAA\n+a\n" + chr(73) * 12 + "\n",
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert err == ""
        header, records = parse(out)
        assert records == []
        assert header[0] == HEADER_PREFIX
        assert header[-1].endswith("SAMPLE")

    def test_output_ends_with_exactly_one_newline(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTTCGTAAAA\n+a\n" + chr(73) * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + chr(73) * 12 + "\n",
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert out.endswith("\n")
        assert not out.endswith("\n\n")

    def test_threshold_options_take_effect(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nAAAAAAAAAAAA\n")
        # One alt out of four: count 1 (<2) and fraction 0.25 (>= 0.2).
        body = (
            "@a\nAAAATAAAAAAA\n+a\n" + chr(73) * 12 + "\n"
            "@b\nAAAAAAAAAAAA\n+b\n" + chr(73) * 12 + "\n"
            "@c\nAAAAAAAAAAAA\n+c\n" + chr(73) * 12 + "\n"
            "@d\nAAAAAAAAAAAA\n+d\n" + chr(73) * 12 + "\n"
        )
        reads = write(tmp_path, "q.fq", body)
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert parse(out)[1] == []
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--min-alt-count",
                "1",
            ]
        )
        assert code == 0
        assert parse(out)[1][0].startswith("c1\t5\t.\tA\tT")

    def test_min_base_quality_option(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        # Variant site carries Phred 10 ('+').
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTTCGTAAAA\n+a\n" + "I" * 4 + "+" + "I" * 7 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + "I" * 4 + "+" + "I" * 7 + "\n",
        )
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 0
        assert parse(out)[1] == []
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--min-base-quality", "10"]
        )
        assert code == 0
        assert len(parse(out)[1]) == 1

    def test_sample_name_option(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTACGTAAAA\n+a\n" + "I" * 12 + "\n",
        )
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
        assert out.splitlines()[-1].endswith("NA12878")

    def test_homozygous_fraction_option(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nAAAAAAAAAAAA\n")
        body = (
            "@a\nAAAATAAAAAAA\n+a\n" + "I" * 12 + "\n"
            "@b\nAAAAAAAAAAAA\n+b\n" + "I" * 12 + "\n"
        )
        reads = write(tmp_path, "q.fq", body)
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--min-alt-count",
                "1",
                "--min-alt-fraction",
                "0.4",
                "--homozygous-fraction",
                "0.4",
            ]
        )
        assert code == 0
        assert parse(out)[1][0].endswith("1/1:2:1,1")

    def test_scoring_options_accepted(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTTCGTAAAA\n+a\n" + "I" * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + "I" * 12 + "\n",
        )
        code, out, err = run(
            [
                "call-variants",
                str(reference),
                str(reads),
                "--match",
                "1",
                "--mismatch",
                "1",
                "--gap-open",
                "0",
                "--gap-extend",
                "0",
                "--min-score",
                "1",
            ]
        )
        assert code == 0
        assert err == ""
        assert parse(out)[0][0] == HEADER_PREFIX

    def test_explicit_fastq_format(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n")
        reads = write(
            tmp_path,
            "q.dat",
            "@a\nACGTTCGTAAAA\n+a\n" + "I" * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + "I" * 12 + "\n",
        )
        code, out, err = run(
            ["call-variants", str(reference), str(reads), "--reads-format", "fastq"]
        )
        assert code == 0
        assert len(parse(out)[1]) == 1


class TestCallVariantsErrors:
    def test_both_inputs_stdin_rejected(self) -> None:
        code, out, err = run(["call-variants", "-", "-"], ">c1\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_empty_reference_exit_2(self, tmp_path) -> None:
        empty = write(tmp_path, "empty.fa", "")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        code, out, err = run(["call-variants", str(empty), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_duplicate_reference_identifier_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n>c1\nTTTT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_quality_length_mismatch_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_read_base_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nAZ\n+a\nII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert "'Z'" in err

    def test_malformed_fastq_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_reference_must_be_fasta(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fq", "@c1\nACGT\n+c1\nIIII\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        code, out, err = run(["call-variants", str(reference), str(reads)])
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_bad_min_base_quality_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        for value in ("-1", "94", "x"):
            code, out, err = run(
                [
                    "call-variants",
                    str(reference),
                    str(reads),
                    "--min-base-quality",
                    value,
                ]
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1

    def test_bad_min_alt_count_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        for value in ("0", "-2", "x"):
            code, out, err = run(
                ["call-variants", str(reference), str(reads), "--min-alt-count", value]
            )
            assert code == 2
            assert out == ""
            assert err.count("\n") == 1

    def test_bad_fractions_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        for option, value in (
            ("--min-alt-fraction", "-0.1"),
            ("--min-alt-fraction", "1.1"),
            ("--min-alt-fraction", "nan"),
            ("--min-alt-fraction", "inf"),
            ("--min-alt-fraction", "x"),
            ("--homozygous-fraction", "2"),
        ):
            code, out, err = run(
                ["call-variants", str(reference), str(reads), option, value]
            )
            assert code == 2, (option, value)
            assert out == ""
            assert err.count("\n") == 1

    def test_homozygous_below_alt_fraction_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
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

    def test_bad_sample_name_exit_2(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
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

    def test_missing_input_file_exit_1(self, tmp_path) -> None:
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
        code, out, err = run(
            ["call-variants", str(tmp_path / "missing.fa"), str(reads)]
        )
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1


class TestCallVariantsOutput:
    def test_file_output_is_atomic_and_byte_stable(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGTACGTAAAA\n>c2\nTTTT\n")
        reads = write(
            tmp_path,
            "q.fq",
            "@a\nACGTTCGTAAAA\n+a\n" + "I" * 12 + "\n"
            "@b\nACGTTCGTAAAA\n+b\n" + "I" * 12 + "\n"
            "@c\nACGTACGTAAAA\n+c\n" + "I" * 12 + "\n",
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
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
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
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n")
        reads = write(tmp_path, "q.fq", "@a\nACGT\n+a\nIIII\n")
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

    def test_other_subcommands_still_listed(self) -> None:
        code, out, err = run([])
        assert code == 0
        assert "call-variants" in out
        assert "map-reads" in out
        assert "align-pair" in out
