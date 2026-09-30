"""Tests for genome_variant.sequence_io."""

from __future__ import annotations

import io
import os

import pytest

from genome_variant.sequence_io import (
    SequenceFormatError,
    SequenceRecord,
    SequenceValidationError,
    read_sequences,
    write_sequences,
)


def parse(text: str, fmt: str = "auto") -> list[SequenceRecord]:
    return list(read_sequences(io.StringIO(text), format=fmt))


class TestFasta:
    def test_basic_multiline_and_fields(self) -> None:
        records = parse(">seq1 desc here\nacgt\nACGT\r\n  NN  \n\n>s2\nACGT\n")
        assert records == [
            SequenceRecord("seq1", "ACGTACGTNN", "desc here"),
            SequenceRecord("s2", "ACGT", ""),
        ]

    def test_all_iupac_symbols(self) -> None:
        sequence = "ACGTRYSWKMBDHVN"
        records = parse(f">x\n{sequence.lower()}\n")
        assert records[0].sequence == sequence

    def test_leading_blank_lines(self) -> None:
        assert parse("\n\n>x\nACGT\n")[0].sequence == "ACGT"

    def test_empty_input(self) -> None:
        with pytest.raises(SequenceFormatError, match="empty input"):
            parse("")
        with pytest.raises(SequenceFormatError, match="empty input"):
            parse("   \n\n")

    def test_empty_identifier_reports_line(self) -> None:
        with pytest.raises(SequenceFormatError) as info:
            parse(">\nACGT\n", "fasta")
        assert "<stream>:1" in str(info.value)
        assert "empty sequence identifier" in str(info.value)

    def test_empty_sequence_reports_header_line(self) -> None:
        with pytest.raises(SequenceFormatError, match=r"<stream>:1: empty sequence"):
            parse(">x\n\n>y\nACGT\n", "fasta")

    def test_content_before_title(self) -> None:
        with pytest.raises(SequenceFormatError) as info:
            parse("ACGT\n>x\nACGT\n", "fasta")
        assert "<stream>:1" in str(info.value)
        assert "content before first title" in str(info.value)

    def test_invalid_base_validation_error(self) -> None:
        with pytest.raises(SequenceValidationError) as info:
            parse(">x\nacgtZ\n", "fasta")
        message = str(info.value)
        assert "'x'" in message
        assert "'Z'" in message
        assert "position 5" in message

    def test_lazy_iteration(self) -> None:
        records = read_sequences(io.StringIO(">ok\nA\n>bad\nZ\n"), format="fasta")
        assert next(records).identifier == "ok"
        with pytest.raises(SequenceValidationError):
            next(records)


class TestFastq:
    def test_multiline_sequence_and_quality(self) -> None:
        text = "@r1 hello\nAC\nGT\n+\nII\nII\n@r2\nA\n+r2\n!\n"
        records = parse(text)
        assert records[0] == SequenceRecord("r1", "ACGT", "hello", (40, 40, 40, 40))
        assert records[1] == SequenceRecord("r2", "A", "", (0,))

    def test_quality_boundaries(self) -> None:
        assert parse("@x\nAC\n+\n!~\n")[0].quality == (0, 93)

    def test_at_sign_quality_line_is_not_a_header(self) -> None:
        records = parse("@a\nAC\n+\n@I\n@b\nA\n+\nI\n")
        assert [r.identifier for r in records] == ["a", "b"]
        assert records[0].quality == (31, 40)

    def test_mismatched_separator_title(self) -> None:
        with pytest.raises(SequenceFormatError, match="separator title"):
            parse("@x desc\nACGT\n+other\nIIII\n")

    @pytest.mark.parametrize(
        "text",
        [
            "@x\nACGT\n+\nIII\n",       # quality too short, truncated record
            "@x\nAC\n+\nIIIII\n",      # quality too long
            "@x\nACGT\ny\nIIII\n",     # missing separator
            "@x\n\n+\nIIII\n",         # empty sequence
            "@x\n+\nIIII\n",           # separator directly after header
        ],
    )
    def test_structure_errors(self, text: str) -> None:
        with pytest.raises(SequenceFormatError):
            parse(text, "fastq")

    def test_length_mismatch_message_has_line(self) -> None:
        with pytest.raises(SequenceFormatError) as info:
            parse("@x\nAC\n+\nIII\n")
        assert "<stream>:4" in str(info.value)

    def test_quality_character_out_of_range(self) -> None:
        with pytest.raises(SequenceFormatError, match="33-126"):
            parse("@x\nAC\n+\nI I\n")

    def test_invalid_base_is_validation_error(self) -> None:
        with pytest.raises(SequenceValidationError, match="position 3"):
            parse("@x\nACZ\n+\nIII\n")


class TestAutoDetection:
    def test_detect_fasta_and_fastq(self) -> None:
        assert parse(">x\nA\n")[0].quality is None
        assert parse("@x\nA\n+\nI\n")[0].quality is not None

    def test_undecidable(self) -> None:
        with pytest.raises(SequenceFormatError) as info:
            parse("ACGT\n")
        assert "<stream>:1" in str(info.value)
        assert "cannot determine format" in str(info.value)

    def test_explicit_format_mismatch(self) -> None:
        with pytest.raises(SequenceFormatError):
            parse(">x\nA\n", "fastq")
        with pytest.raises(SequenceFormatError):
            parse("@x\nA\n+I\n", "fasta")


class TestPaths:
    def test_roundtrip_and_source_name(self, tmp_path) -> None:
        source = tmp_path / "in.fa"
        source.write_text(">x gene\nacgtacgt\n", newline="")
        records = list(read_sequences(source))
        assert records[0].sequence == "ACGTACGT"

        destination = tmp_path / "out.fa"
        write_sequences(records, destination, format="fasta", line_width=5)
        expected = b">x gene\nACGTA\nCGT\n"
        assert destination.read_bytes() == expected
        # Repeated runs are byte-identical.
        write_sequences(records, destination, format="fasta", line_width=5)
        assert destination.read_bytes() == expected

    def test_missing_path_keeps_oserror(self, tmp_path) -> None:
        with pytest.raises(OSError):
            list(read_sequences(tmp_path / "missing.fa"))

    def test_error_message_uses_path(self, tmp_path) -> None:
        source = tmp_path / "in.fa"
        source.write_text("not a title\n")
        with pytest.raises(SequenceFormatError) as info:
            list(read_sequences(source))
        assert os.fspath(source) in str(info.value)

    def test_iteration_errors_propagate(self) -> None:
        class BrokenStream(io.StringIO):
            def __iter__(self):  # type: ignore[override]
                raise RuntimeError("boom")

        with pytest.raises(RuntimeError, match="boom"):
            list(read_sequences(BrokenStream("@x\nAC\n"), "fastq"))


class TestWriter:
    def test_write_fastq_roundtrip(self) -> None:
        stream = io.StringIO()
        write_sequences([SequenceRecord("r1", "AC", "d", (0, 93))], stream, format="fastq")
        assert stream.getvalue() == "@r1 d\nAC\n+r1 d\n!~\n"

    def test_fasta_without_quality_cannot_become_fastq(self) -> None:
        with pytest.raises(ValueError, match="no quality"):
            write_sequences([SequenceRecord("r1", "AC")], io.StringIO(), format="fastq")

    @pytest.mark.parametrize("width", [0, -1, True, 1.5])
    def test_bad_line_width(self, width) -> None:
        with pytest.raises(ValueError):
            write_sequences([], io.StringIO(), line_width=width)
