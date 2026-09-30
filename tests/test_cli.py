"""Tests for the genome-variant-toolkit command line interface."""

from __future__ import annotations

import io
import sys

import pytest

from genome_variant.cli import main
from genome_variant import __version__


class CliHarness:
    def __init__(self, tmp_path, monkeypatch, capsys):
        self.tmp_path = tmp_path
        self.monkeypatch = monkeypatch
        self.capsys = capsys

    def run(self, *argv, stdin=""):
        self.monkeypatch.setattr(sys, "stdin", io.StringIO(stdin))
        try:
            code = main(list(argv))
        except SystemExit as exc:
            # argparse exits directly on --help (code 0); anything else is a
            # sign that error handling leaked out of main().
            if exc.code != 0:
                raise AssertionError(f"main() raised SystemExit({exc.code})")
            code = 0
        captured = self.capsys.readouterr()
        return code, captured.out, captured.err


@pytest.fixture
def cli(tmp_path, monkeypatch, capsys):
    return CliHarness(tmp_path, monkeypatch, capsys)


def test_version(cli):
    code, out, err = cli.run("version")
    assert code == 0
    assert out.strip() == __version__


def test_no_subcommand_prints_help_and_succeeds(cli):
    code, out, err = cli.run()
    assert code == 0
    assert "usage" in out.lower()
    assert "normalize-sequences" in out


def test_help_lists_subcommand(cli):
    code, out, _ = cli.run("--help")
    assert code == 0
    assert "normalize-sequences" in out


# ------------------------------------------------------- normalize: fasta

def test_normalize_fasta_defaults_to_fasta(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">s1 desc\n acgt\nACGT\n", encoding="utf-8")
    code, out, err = cli.run("normalize-sequences", str(source))
    assert code == 0
    assert err == ""
    assert out == ">s1 desc\nACGTACGT\n"


def test_normalize_fasta_line_width(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">s\nACGTACGT\n", encoding="utf-8")
    code, out, _ = cli.run("normalize-sequences", str(source), "--line-width", "3")
    assert code == 0
    assert out == ">s\nACG\nTAC\nGT\n"


def test_normalize_fasta_to_output_file(cli, tmp_path):
    source = tmp_path / "in.fa"
    target = tmp_path / "out.fa"
    source.write_text(">s\nacgt\n", encoding="utf-8")
    code, _, err = cli.run(
        "normalize-sequences", str(source), "--output", str(target)
    )
    assert code == 0
    assert err == ""
    assert target.read_bytes() == b">s\nACGT\n"


def test_normalize_stdin_to_stdout(cli):
    code, out, err = cli.run(
        "normalize-sequences", "-", "--input-format", "fasta", stdin=">s\nacgt\n"
    )
    assert code == 0
    assert out == ">s\nACGT\n"
    assert err == ""


def test_normalize_output_dash_means_stdout(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">s\nac\n", encoding="utf-8")
    code, out, _ = cli.run(
        "normalize-sequences", str(source), "--output", "-"
    )
    assert code == 0
    assert out == ">s\nAC\n"


def test_normalize_fastq_keeps_fastq_layout(cli, tmp_path):
    source = tmp_path / "in.fq"
    source.write_text("@r1 d\nacgt\n+r1 d\n!!II\n", encoding="utf-8")
    code, out, _ = cli.run("normalize-sequences", str(source))
    assert code == 0
    assert out == "@r1 d\nACGT\n+\n" + chr(33) * 2 + "II\n"
    assert len(out.splitlines()) == 4


def test_normalize_fastq_line_width_is_ignored(cli, tmp_path):
    source = tmp_path / "in.fq"
    source.write_text("@r\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")
    code, out, _ = cli.run(
        "normalize-sequences", str(source), "--line-width", "3"
    )
    assert code == 0
    assert out.splitlines() == [
        "@r", "ACGTACGT", "+", "IIIIIIII"
    ]


def test_normalize_fastq_input_to_fasta_output(cli, tmp_path):
    source = tmp_path / "in.fq"
    source.write_text("@r\nacgt\n+\nIIII\n", encoding="utf-8")
    code, out, _ = cli.run(
        "normalize-sequences", str(source),
        "--output-format", "fasta", "--line-width", "2",
    )
    assert code == 0
    assert out == ">r\nAC\nGT\n"


def test_normalize_fasta_to_fastq_is_usage_error(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">r\nACGT\n", encoding="utf-8")
    code, out, err = cli.run(
        "normalize-sequences", str(source), "--output-format", "fastq"
    )
    assert code == 2
    assert out == ""
    assert err.count("\n") == 1
    assert "quality" in err.lower()


def test_normalize_bad_line_width_exits_2(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">r\nACGT\n", encoding="utf-8")
    for bad in ("0", "-3", "wide"):
        code, _, err = cli.run(
            "normalize-sequences", str(source), "--line-width", bad
        )
        assert code == 2, bad
        assert err.count("\n") == 1


def test_normalize_unknown_format_exits_2(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">r\nACGT\n", encoding="utf-8")
    code, _, err = cli.run(
        "normalize-sequences", str(source), "--input-format", "genbank"
    )
    assert code == 2
    assert err.count("\n") == 1


def test_normalize_missing_input_file_exits_1(cli, tmp_path):
    code, out, err = cli.run("normalize-sequences", str(tmp_path / "missing.fa"))
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1
    assert "missing.fa" in err


def test_normalize_bad_output_dir_exits_1(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">r\nACGT\n", encoding="utf-8")
    code, _, err = cli.run(
        "normalize-sequences", str(source),
        "--output", str(tmp_path / "no_such_dir" / "out.fa"),
    )
    assert code == 1
    assert err.endswith("\n")


def test_normalize_malformed_input_exits_2(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">r\nACGTX\n", encoding="utf-8")
    code, out, err = cli.run("normalize-sequences", str(source))
    assert code == 2
    assert out == ""
    assert "r" in err and "X" in err and ":2:" in err
    assert err.count("\n") == 1


def test_normalize_empty_input_exits_2(cli):
    code, _, err = cli.run(
        "normalize-sequences", "-", "--input-format", "auto", stdin="\n  \n"
    )
    assert code == 2
    assert err.endswith("\n") and err.count("\n") == 1


def test_normalize_undecidable_input_exits_2(cli):
    code, _, err = cli.run("normalize-sequences", "-", stdin="ACGT\n")
    assert code == 2
    assert "stdin" in err or "<stdin>" in err


def test_normalize_deterministic_byte_output(cli, tmp_path):
    source = tmp_path / "in.fq"
    payload = "@r1\nacgt\n+\nIIII\n@r2 xy\nnn\n+\n##\n"
    source.write_text(payload, encoding="utf-8")
    outs = []
    for _ in range(2):
        code, out, _ = cli.run("normalize-sequences", str(source))
        assert code == 0
        outs.append(out)
    assert outs[0] == outs[1]
    assert outs[0].endswith("\n") and not outs[0].endswith("\n\n")


def test_normalize_trailing_newline_exactly_one(cli, tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">s\nAC\n", encoding="utf-8")
    target = tmp_path / "o.fa"
    code, _, _ = cli.run(
        "normalize-sequences", str(source), "--output", str(target)
    )
    assert code == 0
    data = target.read_bytes()
    assert data.endswith(b"\n") and not data.endswith(b"\n\n")
    assert b"\r" not in data
