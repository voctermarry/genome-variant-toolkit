"""Tests for the genome-variant-toolkit command line interface."""

from __future__ import annotations

import json
import sys

import pytest

from genome_variant.cli import main


def run(argv, stdin: str = ""):
    from io import StringIO

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


def test_no_subcommand_prints_help_and_succeeds() -> None:
    code, out, err = run([])
    assert code == 0
    assert "usage:" in out
    assert err == ""


def test_version() -> None:
    code, out, err = run(["version"])
    assert code == 0
    assert out.strip() == "0.1.0"
    assert err == ""


def test_normalize_fasta_stdin_wrapping() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--line-width", "4"],
        ">a gene desc\nacgt\naCGT\n",
    )
    assert code == 0
    assert err == ""
    assert out == ">a gene desc\nACGT\nACGT\n"


def test_normalize_fastq_fixed_layout_ignores_width() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--line-width", "2"],
        "@r1\nAC\nGT\n+r1\nIIII\n@r2\nA\n+\n!\n",
    )
    assert code == 0
    assert out == "@r1\nACGT\n+r1\nIIII\n@r2\nA\n+r2\n!\n"


def test_fasta_to_fastq_is_usage_error() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--output-format", "fastq"],
        ">a\nACGT\n",
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert "no quality" in err


def test_fastq_to_fasta_conversion() -> None:
    code, out, err = run(
        ["normalize-sequences", "-", "--output-format", "fasta"],
        "@x\nacgt\n+x\nIIII\n",
    )
    assert code == 0
    assert out == ">x\nACGT\n"


def test_bad_line_width_exit_2_single_line() -> None:
    for value in ("0", "-3", "abc"):
        code, out, err = run(["normalize-sequences", "-", "--line-width", value], ">a\nA\n")
        assert code == 2
        assert err.count("\n") == 1
        assert out == ""


def test_format_errors_exit_2() -> None:
    code, _, err = run(["normalize-sequences", "-"], "")
    assert code == 2 and err.endswith("\n") and err.count("\n") == 1

    code, _, err = run(["normalize-sequences", "-"], ">a\nACGZ\n")
    assert code == 2 and "position 4" in err

    code, _, err = run(["normalize-sequences", "-"], "@a\nACGT\n+\nIII\n")
    assert code == 2 and "incomplete" in err


def test_input_file_error_exit_1(tmp_path) -> None:
    code, out, err = run(["normalize-sequences", str(tmp_path / "missing.fa")])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1


def test_output_file_error_exit_1(tmp_path) -> None:
    code, _, err = run(
        [
            "normalize-sequences",
            "-",
            "--output",
            str(tmp_path / "no-such-dir" / "out.fa"),
        ],
        ">a\nACGT\n",
    )
    assert code == 1


def test_file_roundtrip_is_byte_stable(tmp_path) -> None:
    source = tmp_path / "in.fq"
    source.write_text("@x\nacgt\n+x\nIIII\n")
    first = tmp_path / "out1.fa"
    second = tmp_path / "out2.fa"
    code1, _, err1 = run(
        ["normalize-sequences", str(source), "--output", str(first)]
    )
    code2, _, err2 = run(
        ["normalize-sequences", str(source), "--output", str(second)]
    )
    assert code1 == code2 == 0
    assert err1 == err2 == ""
    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes() == b"@x\nACGT\n+x\nIIII\n"


def test_unknown_argument_exit_2_single_line() -> None:
    code, out, err = run(["normalize-sequences", "-", "--bogus"], ">a\nA\n")
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1


def q(scores) -> str:
    return "".join(chr(s + 33) for s in scores)


class TestFilterReads:
    def test_basic_filter_stdout(self) -> None:
        kept = "@r1 desc\n" + "A" * 30 + "\n+r1 desc\n" + q([40] * 30) + "\n"
        dropped = "@r2\nACGT\n+r2\n" + q([1] * 4) + "\n"
        code, out, err = run(["filter-reads", "-"], kept + dropped)
        assert code == 0
        assert err == ""
        assert out == kept

    def test_end_trimming_is_reflected_in_output(self) -> None:
        text = "@r\nACGTACGT\n+r\n" + q([19, 20, 30, 40, 40, 30, 20, 19]) + "\n"
        code, out, err = run(
            [
                "filter-reads", "-",
                "--min-end-quality", "20",
                "--min-mean-quality", "0",
                "--min-length", "1",
            ],
            text,
        )
        assert code == 0
        assert err == ""
        assert out == "@r\nCGTACG\n+r\n" + q([20, 30, 40, 40, 30, 20]) + "\n"

    def test_all_filtered_writes_empty_output(self) -> None:
        text = "@r\nACGT\n+r\n" + q([1] * 4) + "\n"
        code, out, err = run(["filter-reads", "-"], text)
        assert code == 0
        assert out == ""
        assert err == ""

    def test_order_preserved(self) -> None:
        a = "@a\n" + "A" * 30 + "\n+a\n" + q([40] * 30) + "\n"
        b = "@b\n" + "C" * 30 + "\n+b\n" + q([41] * 30) + "\n"
        c = "@c\n" + "G" * 30 + "\n+c\n" + q([42] * 30) + "\n"
        code, out, err = run(["filter-reads", "-"], b + a + c)
        assert code == 0
        assert [line for line in out.splitlines() if line[0] == "@"] == ["@b", "@a", "@c"]

    def test_fasta_input_is_rejected_exit_2(self) -> None:
        code, out, err = run(["filter-reads", "-"], ">a\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_malformed_fastq_exit_2(self) -> None:
        code, out, err = run(["filter-reads", "-"], "@x\nACGT\n+\nIII\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2(self) -> None:
        code, out, err = run(["filter-reads", "-"], "@x\nACGZ\n+\n" + q([40] * 4) + "\n")
        assert code == 2
        assert "position 4" in err

    def test_missing_input_file_exit_1(self, tmp_path) -> None:
        code, out, err = run(["filter-reads", str(tmp_path / "missing.fq")])
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    @pytest.mark.parametrize(
        "option,value",
        [
            ("--min-end-quality", "-1"),
            ("--min-end-quality", "94"),
            ("--min-end-quality", "x"),
            ("--min-mean-quality", "100"),
            ("--min-length", "0"),
            ("--min-length", "-2"),
        ],
    )
    def test_bad_thresholds_exit_2_single_line(self, option, value) -> None:
        text = "@r\nACGT\n+r\n" + q([40] * 4) + "\n"
        code, out, err = run(["filter-reads", "-", option, value], text)
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_output_is_byte_stable_and_uses_single_trailing_newline(self, tmp_path) -> None:
        text = "@r\n" + "A" * 30 + "\n+r\n" + q([40] * 30) + "\n"
        first = tmp_path / "o1.fq"
        second = tmp_path / "o2.fq"
        code1, _, err1 = run(["filter-reads", "-", "--output", str(first)], text)
        code2, _, err2 = run(["filter-reads", "-", "--output", str(second)], text)
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        assert first.read_bytes() == second.read_bytes()
        assert first.read_bytes().endswith(b"\n")
        assert not first.read_bytes().endswith(b"\n\n")

    def test_file_output_error_preserves_existing_target(self, tmp_path) -> None:
        target = tmp_path / "out.fq"
        target.write_text("PREVIOUS CONTENT\n")
        # Read error on the second record must not replace the target.
        first = "@r1\n" + "A" * 30 + "\n+r1\n" + q([40] * 30) + "\n"
        bad = "@r2\nACGT\n+r2\n" + q([1, 2]) + "\n"
        code, out, err = run(
            ["filter-reads", "-", "--output", str(target)], first + bad
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS CONTENT\n"

    def test_format_error_preserves_existing_target(self, tmp_path) -> None:
        target = tmp_path / "out.fq"
        target.write_bytes(b"keep me\n")
        code, _, _ = run(
            ["filter-reads", "-", "--output", str(target)],
            "@r\nACGT\n+\nII\n",
        )
        assert code == 2
        assert target.read_bytes() == b"keep me\n"

    def test_all_filtered_leaves_empty_file(self, tmp_path) -> None:
        source = tmp_path / "in.fq"
        source.write_text("@r\nACGT\n+r\n" + q([1] * 4) + "\n")
        target = tmp_path / "out.fq"
        target.write_text("old\n")
        code, _, err = run(
            ["filter-reads", str(source), "--output", str(target)]
        )
        assert code == 0
        assert err == ""
        assert target.read_bytes() == b""

    def test_input_output_file_roundtrip(self, tmp_path) -> None:
        source = tmp_path / "in.fq"
        source.write_text("@r desc\n" + "A" * 30 + "\n+r desc\n" + q([40] * 30) + "\n")
        target = tmp_path / "out.fq"
        code, _, err = run(
            ["filter-reads", str(source), "--output", str(target)]
        )
        assert code == 0
        assert target.read_bytes() == source.read_bytes()


class TestKmerIndex:
    def test_canonical_jsonl_stdout(self) -> None:
        code, out, err = run(["kmer-index", "-", "--k", "2"], ">r\nACGT\n")
        assert code == 0
        assert err == ""
        lines = out.splitlines()
        assert [json.loads(line)["kmer"] for line in lines] == ["AC", "CG"]
        first = json.loads(lines[0])
        assert set(first) == {"kmer", "count", "occurrences"}
        assert first["count"] == 2
        assert first["occurrences"] == [
            {"record": 0, "id": "r", "position": 0, "strand": "+"},
            {"record": 0, "id": "r", "position": 2, "strand": "-"},
        ]
        for occurrence in first["occurrences"]:
            assert set(occurrence) == {"record", "id", "position", "strand"}

    def test_compact_json_has_no_spaces(self) -> None:
        code, out, err = run(["kmer-index", "-", "--k", "2"], ">r\nACGT\n")
        assert code == 0
        assert '"count": 2' not in out
        assert '"count":2' in out
        assert ", " not in out
        assert ": " not in out

    def test_lines_sorted_and_count_matches_occurrences(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--no-canonical"],
            ">r\nTTTTGGGGAAAACCCC\n",
        )
        assert code == 0
        rows = [json.loads(line) for line in out.splitlines()]
        keys = [row["kmer"] for row in rows]
        assert keys == sorted(keys)
        for row in rows:
            assert row["count"] == len(row["occurrences"])

    def test_non_canonical_flag(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--no-canonical"], ">r\nACGT\n"
        )
        assert code == 0
        rows = [json.loads(line) for line in out.splitlines()]
        assert [row["kmer"] for row in rows] == ["AC", "CG", "GT"]
        assert all(o["strand"] == "+" for row in rows for o in row["occurrences"])

    def test_palindrome_is_plus(self) -> None:
        code, out, err = run(["kmer-index", "-", "--k", "4"], ">p\nACGT\n")
        assert code == 0
        row = json.loads(out)
        assert row["kmer"] == "ACGT"
        assert row["occurrences"] == [{"record": 0, "id": "p", "position": 0, "strand": "+"}]

    def test_fastq_input_and_quality_is_ignored(self) -> None:
        # Terrible quality scores must not drop or alter any record.
        low = q([1, 1, 1, 1])
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--input-format", "fastq"],
            f"@r1\nACGT\n+r1\n{low}\n@r2\nacgt\n+r2\n{low}\n",
        )
        assert code == 0
        assert err == ""
        rows = [json.loads(line) for line in out.splitlines()]
        totals = {row["kmer"]: row["count"] for row in rows}
        assert totals == {"AC": 4, "CG": 2}

    def test_auto_detects_fastq(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2"], "@r\nACGT\n+r\n" + q([40] * 4) + "\n"
        )
        assert code == 0
        assert [json.loads(line)["kmer"] for line in out.splitlines()] == ["AC", "CG"]

    def test_duplicate_identifiers_not_merged(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--no-canonical"],
            ">x\nAA\n>x\nAA\n",
        )
        assert code == 0
        row = json.loads(out)
        assert row["kmer"] == "AA"
        assert [o["record"] for o in row["occurrences"]] == [0, 1]
        assert [o["id"] for o in row["occurrences"]] == ["x", "x"]

    def test_non_acgt_windows_skipped_and_short_records(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--no-canonical"],
            ">a\nACGNTAC\n>short\nA\n",
        )
        assert code == 0
        rows = {row["kmer"]: row for row in (json.loads(line) for line in out.splitlines())}
        assert set(rows) == {"AC", "CG", "TA"}
        assert [o["position"] for o in rows["AC"]["occurrences"]] == [0, 5]
        assert all(o["record"] == 0 for o in rows["AC"]["occurrences"])

    def test_empty_index_writes_nothing(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "3"], ">a\nNN\n>b\nAC\n"
        )
        assert code == 0
        assert out == ""
        assert err == ""

    def test_non_ascii_identifier_preserved(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2"], ">gène-甲\nACGT\n"
        )
        assert code == 0
        assert "\\u" not in out
        assert "gène-甲" in out

    def test_single_trailing_newline(self) -> None:
        code, out, err = run(["kmer-index", "-", "--k", "2"], ">r\nACGT\n")
        assert code == 0
        assert out.endswith("\n")
        assert not out.endswith("\n\n")

    def test_k_required(self) -> None:
        code, out, err = run(["kmer-index", "-"], ">r\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    @pytest.mark.parametrize("value", ["0", "-3", "x"])
    def test_bad_k_exit_2(self, value) -> None:
        code, out, err = run(["kmer-index", "-", "--k", value], ">r\nACGT\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_invalid_base_exit_2_single_line(self) -> None:
        code, out, err = run(["kmer-index", "-", "--k", "2"], ">r\nACGZ\n")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert "position 4" in err
        assert "'r'" in err

    def test_format_error_exit_2(self) -> None:
        code, out, err = run(["kmer-index", "-", "--k", "2"], "")
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_wrong_format_flag_exit_2(self) -> None:
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--input-format", "fastq"],
            ">r\nACGT\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1

    def test_missing_input_file_exit_1(self, tmp_path) -> None:
        code, out, err = run(["kmer-index", str(tmp_path / "missing.fa"), "--k", "2"])
        assert code == 1
        assert out == ""
        assert err.count("\n") == 1

    def test_file_output_atomic_and_byte_stable(self, tmp_path) -> None:
        target_a = tmp_path / "a.jsonl"
        target_b = tmp_path / "b.jsonl"
        text = ">r1\nACGTACGT\n>r2 desc\nTTGG\n"
        code1, _, err1 = run(
            ["kmer-index", "-", "--k", "3", "--output", str(target_a)], text
        )
        code2, _, err2 = run(
            ["kmer-index", "-", "--k", "3", "--output", str(target_b)], text
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        assert target_a.read_bytes() == target_b.read_bytes()
        data = target_a.read_bytes()
        assert data.endswith(b"\n")
        assert not data.endswith(b"\n\n")

    def test_file_output_preserves_target_on_sequence_error(self, tmp_path) -> None:
        target = tmp_path / "out.jsonl"
        target.write_text("PREVIOUS CONTENT\n")
        code, out, err = run(
            ["kmer-index", "-", "--k", "2", "--output", str(target)],
            ">r1\nACGT\n>r2\nACGZ\n",
        )
        assert code == 2
        assert out == ""
        assert err.count("\n") == 1
        assert target.read_text() == "PREVIOUS CONTENT\n"

    def test_empty_index_truncates_target_to_empty(self, tmp_path) -> None:
        target = tmp_path / "out.jsonl"
        target.write_text("old\n")
        code, _, err = run(
            ["kmer-index", "-", "--k", "3", "--output", str(target)],
            ">a\nNN\n",
        )
        assert code == 0
        assert err == ""
        assert target.read_bytes() == b""

    def test_deterministic_across_files_and_records(self, tmp_path) -> None:
        source_one = tmp_path / "one.fa"
        source_two = tmp_path / "two.fa"
        source_one.write_text(">a\nACGT\n>b\nGGGT\n")
        source_two.write_text(">c\nACGT\n>d\nGGGT\n")
        out_one = tmp_path / "one.jsonl"
        out_two = tmp_path / "two.jsonl"
        code1, _, err1 = run(
            ["kmer-index", str(source_one), "--k", "2", "--output", str(out_one)]
        )
        code2, _, err2 = run(
            ["kmer-index", str(source_two), "--k", "2", "--output", str(out_two)]
        )
        assert code1 == code2 == 0
        assert err1 == err2 == ""
        # Keys, positions and strands are identical; only record ids differ.
        rows_one = [json.loads(line) for line in out_one.read_text().splitlines()]
        rows_two = [json.loads(line) for line in out_two.read_text().splitlines()]
        assert [r["kmer"] for r in rows_one] == [r["kmer"] for r in rows_two]
        stripped_one = [
            (r["kmer"], [(o["record"], o["position"], o["strand"]) for o in r["occurrences"]])
            for r in rows_one
        ]
        stripped_two = [
            (r["kmer"], [(o["record"], o["position"], o["strand"]) for o in r["occurrences"]])
            for r in rows_two
        ]
        assert stripped_one == stripped_two
