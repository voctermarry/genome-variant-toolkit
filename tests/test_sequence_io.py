"""Tests for genome_variant.sequence_io."""

from __future__ import annotations

import io

import pytest

from genome_variant.sequence_io import (
    AUTO,
    FASTA,
    FASTQ,
    SequenceFormatError,
    SequenceRecord,
    SequenceValidationError,
    read_sequences,
    write_sequences,
)


def read_all(text, format=AUTO, *, name="<test>"):
    return list(read_sequences(io.StringIO(text), format, source_name=name))


# ---------------------------------------------------------------- FASTA

def test_fasta_basic_multiline():
    text = ">seq1 first record\nacgt\nACGT\n>seq2\nnnn\n"
    records = read_all(text, FASTA)
    assert records == [
        SequenceRecord("seq1", "ACGTACGT", description="first record"),
        SequenceRecord("seq2", "NNN", description=""),
    ]


def test_fasta_whitespace_inside_sequence_lines_is_stripped():
    text = ">s1\n a C\ttg \n"
    records = read_all(text, FASTA)
    assert records[0].sequence == "ACTG"


def test_fasta_blank_lines_are_ignored():
    text = "\n\n>s1\nACGT\n\n>s2\r\nNN\r\n"
    records = read_all(text, FASTA)
    assert [r.id for r in records] == ["s1", "s2"]
    assert records[1].sequence == "NN"


def test_fasta_all_iupac_symbols_accepted():
    text = ">s\nACGTRYSWKMBDHVN\n"
    records = read_all(text, FASTA)
    assert records[0].sequence == "ACGTRYSWKMBDHVN"


def test_fasta_invalid_base_message():
    text = ">s1\nACGTX\n"
    with pytest.raises(SequenceValidationError) as excinfo:
        read_all(text, FASTA)
    message = str(excinfo.value)
    assert "s1" in message
    assert "'X'" in message
    assert "position 5" in message


def test_fasta_invalid_base_position_spans_multiple_lines():
    text = ">s1\nAC\nG?T\n"
    with pytest.raises(SequenceValidationError) as excinfo:
        read_all(text, FASTA)
    assert "position 4" in str(excinfo.value)


def test_fasta_empty_identifier():
    with pytest.raises(SequenceFormatError, match="empty record identifier"):
        read_all(">   \nACGT\n", FASTA)


def test_fasta_empty_sequence():
    with pytest.raises(SequenceFormatError, match="empty sequence"):
        read_all(">s1\n>s2\nACGT\n", FASTA)


def test_fasta_content_before_header():
    text = "ACGT\n>s1\nACGT\n"
    with pytest.raises(SequenceFormatError) as excinfo:
        read_all(text, FASTA)
    assert ":1:" in str(excinfo.value)
    assert "'>'" in str(excinfo.value)


def test_fasta_empty_input():
    with pytest.raises(SequenceFormatError, match="empty input"):
        read_all("   \n\n", FASTA)


def test_fasta_lazy_iteration():
    stream = io.StringIO(">a\nAC\n>b\n!!\n")
    iterator = read_sequences(stream, FASTA, source_name="lazy")
    first = next(iterator)
    assert first.id == "a"
    with pytest.raises(SequenceValidationError):
        list(iterator)


# ---------------------------------------------------------------- FASTQ

def test_fastq_basic_four_lines():
    text = "@read1 desc\nacgt\n+\nIIII\n"
    records = read_all(text, FASTQ)
    assert len(records) == 1
    record = records[0]
    assert record.id == "read1"
    assert record.description == "desc"
    assert record.sequence == "ACGT"
    assert record.quality == (40, 40, 40, 40)


def test_fastq_repeated_title_must_match_exactly():
    text = "@read1 desc\nACGT\n+read1 desc\nIIII\n"
    records = read_all(text, FASTQ)
    assert records[0].id == "read1"


def test_fastq_repeated_title_mismatch():
    text = "@read1 desc\nACGT\n+read2 desc\nIIII\n"
    with pytest.raises(SequenceFormatError, match="does not match"):
        read_all(text, FASTQ)


def test_fastq_quality_range_boundaries():
    text = "@r1\nAC\n+\n" + chr(33) + chr(126) + "\n"
    records = read_all(text, FASTQ)
    assert records[0].quality == (0, 93)


def test_fastq_quality_out_of_range():
    text = "@r1\nAC\n+\n" + chr(32) + chr(126) + "\n"
    with pytest.raises(SequenceFormatError, match="quality character"):
        read_all(text, FASTQ)


def test_fastq_multiline_sequence_and_quality():
    text = "@r1\nAC\nGT\n+\nII\nII\n"
    records = read_all(text, FASTQ)
    assert records[0].sequence == "ACGT"
    assert records[0].quality == (40,) * 4


def test_fastq_quality_longer_than_sequence():
    text = "@r1\nACGT\n+\nIIIII\n"
    with pytest.raises(SequenceFormatError, match="quality length 5"):
        read_all(text, FASTQ)


def test_fastq_quality_shorter_than_sequence():
    text = "@r1\nACGT\n+\nIII\n"
    with pytest.raises(SequenceFormatError, match="does not match"):
        read_all(text, FASTQ)


def test_fastq_two_records_and_lazy_error():
    text = "@a\nAC\n+\nII\n@b\nX\n+\nI\n"
    iterator = read_sequences(io.StringIO(text), FASTQ, source_name="mq")
    first = next(iterator)
    assert first.id == "a"
    with pytest.raises(SequenceValidationError):
        next(iterator)


def test_fastq_truncated_missing_separator():
    with pytest.raises(SequenceFormatError, match="missing '\\+' separator"):
        read_all("@r1\nACGT\n", FASTQ)


def test_fastq_truncated_missing_quality():
    with pytest.raises(SequenceFormatError, match="incomplete FASTQ record"):
        read_all("@r1\nACGT\n+\nII\n", FASTQ)


def test_fastq_header_does_not_start_with_at():
    text = "r1\nACGT\n+\nIIII\n"
    with pytest.raises(SequenceFormatError) as excinfo:
        read_all(text, FASTQ)
    assert "'@'" in str(excinfo.value)


def test_fastq_empty_input():
    with pytest.raises(SequenceFormatError, match="empty input"):
        read_all("", FASTQ)


def test_fastq_empty_sequence():
    with pytest.raises(SequenceFormatError, match="empty sequence"):
        read_all("@r1\n+\n\n", FASTQ)


def test_fastq_empty_identifier():
    with pytest.raises(SequenceFormatError, match="empty record identifier"):
        read_all("@\nACGT\n+\nIIII\n", FASTQ)


def test_fastq_invalid_base_message():
    with pytest.raises(SequenceValidationError) as excinfo:
        read_all("@r1 x\nACZ\n+\nIII\n", FASTQ)
    message = str(excinfo.value)
    assert "r1" in message and "'Z'" in message and "position 3" in message


# ---------------------------------------------------------------- auto

def test_auto_detects_fasta():
    records = read_all(">s\nACGT\n")
    assert records[0].quality is None


def test_auto_detects_fastq():
    records = read_all("@s\nACGT\n+\nIIII\n")
    assert records[0].quality is not None


def test_auto_skips_leading_blank_lines():
    records = read_all("\n\n>s\nAC\n")
    assert records[0].id == "s"


def test_auto_undecidable():
    with pytest.raises(SequenceFormatError, match="cannot determine"):
        read_all("ACGT\n")


def test_auto_empty():
    with pytest.raises(SequenceFormatError, match="empty input"):
        read_all("\n  \n")


# ---------------------------------------------------------------- paths / sources

def test_read_from_path(tmp_path):
    path = tmp_path / "in.fa"
    path.write_text(">s\nACGT\n", encoding="utf-8")
    records = list(read_sequences(str(path)))
    assert records[0].id == "s"


def test_read_missing_path_raises_oserror(tmp_path):
    missing = tmp_path / "nope.fa"
    with pytest.raises(OSError):
        list(read_sequences(str(missing)))


def test_error_message_contains_source_and_line(tmp_path):
    path = tmp_path / "bad.fa"
    path.write_text(">s\nAC\n!!\n", encoding="utf-8")
    with pytest.raises(SequenceValidationError) as excinfo:
        list(read_sequences(str(path), FASTA))
    assert str(path) in str(excinfo.value)
    assert ":3:" in str(excinfo.value)


def test_unknown_format_argument():
    with pytest.raises(ValueError, match="unknown format"):
        list(read_sequences(io.StringIO(">s\nA\n"), format="genbank"))


# ---------------------------------------------------------------- writing

def test_write_fasta_wraps_at_line_width():
    record = SequenceRecord("s1", "ACGTACGT", description="d")
    stream = io.StringIO()
    write_sequences([record], stream, FASTA, line_width=3)
    assert stream.getvalue() == ">s1 d\nACG\nTAC\nGT\n"


def test_write_fasta_default_width():
    record = SequenceRecord("s1", "A" * 81)
    stream = io.StringIO()
    write_sequences([record], stream, FASTA)
    lines = stream.getvalue().split("\n")
    assert lines[0] == ">s1"
    assert len(lines[1]) == 80 and len(lines[2]) == 1 and lines[3] == ""


def test_write_fastq_roundtrip():
    record = SequenceRecord("r", "ACGT", quality=(0, 10, 40, 93))
    stream = io.StringIO()
    write_sequences([record], stream, FASTQ)
    assert stream.getvalue() == "@r\nACGT\n+\n" + chr(33) + chr(43) + "I" + chr(126) + "\n"


def test_write_fasta_record_without_quality_to_fastq_is_value_error():
    record = SequenceRecord("r", "ACGT")
    with pytest.raises(ValueError, match="no quality scores"):
        write_sequences([record], io.StringIO(), FASTQ)


def test_write_invalid_line_width():
    record = SequenceRecord("r", "ACGT")
    for width in (0, -1):
        with pytest.raises(ValueError, match="positive integer"):
            write_sequences([record], io.StringIO(), FASTA, line_width=width)


def test_write_quality_out_of_range():
    record = SequenceRecord("r", "A", quality=(94,))
    with pytest.raises(ValueError, match="Phred\\+33 range"):
        write_sequences([record], io.StringIO(), FASTQ)


def test_roundtrip_bytes_are_stable(tmp_path):
    original = ">s1 alpha\nacgtacgt\n>s2\nnnn\n"
    path = tmp_path / "x.fa"
    path.write_text(original, encoding="utf-8")
    out1 = tmp_path / "o1.fa"
    out2 = tmp_path / "o2.fa"
    for out in (out1, out2):
        write_sequences(list(read_sequences(str(path))), str(out), FASTA, line_width=5)
    assert out1.read_bytes() == out2.read_bytes()
    assert out1.read_text().endswith("\n") and not out1.read_text().endswith("\n\n")


def test_write_uses_lf_newlines(tmp_path):
    path = tmp_path / "x.fa"
    write_sequences(
        [SequenceRecord("s", "ACGT")], str(path), FASTA, line_width=2
    )
    data = path.read_bytes()
    assert b"\r" not in data
    assert data == b">s\nAC\nGT\n"
