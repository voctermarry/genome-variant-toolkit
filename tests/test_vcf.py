"""Tests for the genome_variant.vcf module."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import (
    ReferenceMismatchError,
    VcfFormatError,
    VcfHeader,
    VcfRecord,
    normalize_vcf,
    read_vcf,
    write_vcf,
)

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
FULL_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"


def parse(text: str):
    return read_vcf(StringIO(text))


def norm(line: str, reference, source: str = "test.vcf") -> VcfRecord:
    _, records = parse(HEADER + "\n" + line + "\n")
    return normalize_vcf(records, reference, source=source)[0]


class TestReadVcf:
    def test_meta_lines_preserved_verbatim_and_in_order(self) -> None:
        header, records = parse(
            "##fileformat=VCFv4.2\n"
            "##source=my source, kept as-is\n"
            "##contig=<ID=chr1,length=9>\n"
            + HEADER
            + "\n"
        )
        assert header.meta_lines == (
            "##fileformat=VCFv4.2",
            "##source=my source, kept as-is",
            "##contig=<ID=chr1,length=9>",
        )
        assert header.columns == (
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        )
        assert header.sample_names == ()
        assert records == []

    def test_records_in_input_order_with_samples(self) -> None:
        header, records = parse(
            FULL_HEADER
            + "\n"
            + "chr1\t3\trs1\tA\tG\t50\tPASS\tDP=4\tGT:DP\t0/1:4\t0/0:4\n"
            + "chr1\t1\t.\tAAA\tAA,A\t.\t.\t.\tGT\t0/2\t1/2\n"
        )
        assert header.sample_names == ("s1", "s2")
        first, second = records
        assert first.chrom == "chr1"
        assert first.pos == 3
        assert first.id == "rs1"
        assert first.ref == "A"
        assert first.alts == ("G",)
        assert first.qual == "50"
        assert first.filter == "PASS"
        assert first.info == "DP=4"
        assert first.format == "GT:DP"
        assert first.samples == ("0/1:4", "0/0:4")
        assert first.line_number == 2
        assert second.alts == ("AA", "A")
        assert second.line_number == 3

    def test_sequence_alleles_are_uppercased(self) -> None:
        _, records = parse(HEADER + "\nchr1\t2\tx\taAg\taGt,n\t.\t.\t.\n")
        assert records[0].ref == "AAG"
        assert records[0].alts == ("AGT", "N")

    def test_special_alleles_kept_verbatim(self) -> None:
        _, records = parse(
            HEADER + "\nchr1\t2\tx\tAA\t<DEL>,A]chr2:5],*,.\t.\t.\t.\n"
        )
        assert records[0].alts == ("<DEL>", "A]chr2:5]", "*", ".")

    def test_missing_column_header(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse("##fileformat=VCFv4.2\n")
        assert ":2:" in str(excinfo.value)

    def test_empty_input_reports_line_one(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse("")
        assert ":1:" in str(excinfo.value)

    def test_record_before_column_header(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse("chr1\t1\t.\tA\tG\t.\t.\t.\n")
        assert ":1:" in str(excinfo.value)

    def test_duplicate_column_header(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(HEADER + "\n" + HEADER + "\n")
        assert ":2:" in str(excinfo.value)
        assert "duplicate" in str(excinfo.value)

    def test_invalid_fixed_columns(self) -> None:
        with pytest.raises(VcfFormatError):
            parse("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tWRONG\n")

    def test_ninth_column_must_be_format(self) -> None:
        with pytest.raises(VcfFormatError):
            parse(HEADER + "\tSAMPLE\ts1\n")

    def test_meta_line_after_column_header(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(HEADER + "\n##late=1\n")
        assert ":2:" in str(excinfo.value)

    def test_column_count_mismatch(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\n")
        assert ":2:" in str(excinfo.value)
        assert "7" in str(excinfo.value) and "8" in str(excinfo.value)

    def test_sample_columns_must_match_header(self) -> None:
        with pytest.raises(VcfFormatError):
            parse(FULL_HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")

    @pytest.mark.parametrize("pos", ["0", "-1", "x", "1.5", "+3", ""])
    def test_non_positive_integer_pos(self, pos: str) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(f"{HEADER}\nchr1\t{pos}\t.\tA\tG\t.\t.\t.\n")
        assert ":2:" in str(excinfo.value)

    @pytest.mark.parametrize("ref", ["", ".", "AX", "A*", "<REF>"])
    def test_empty_or_invalid_ref(self, ref: str) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(f"{HEADER}\nchr1\t1\t.\t{ref}\tG\t.\t.\t.\n")
        assert ":2:" in str(excinfo.value)

    @pytest.mark.parametrize("alt", ["", "A,", ",A", "AX", "A<", "A<T"])
    def test_empty_or_invalid_alt(self, alt: str) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(f"{HEADER}\nchr1\t1\t.\tA\t{alt}\t.\t.\t.\n")
        assert ":2:" in str(excinfo.value)

    def test_error_message_contains_source_name(self) -> None:
        stream = StringIO("chr1\t1\t.\tA\tG\t.\t.\t.\n")
        stream.name = "sample.vcf"
        with pytest.raises(VcfFormatError) as excinfo:
            read_vcf(stream)
        assert str(excinfo.value).startswith("sample.vcf:1:")

    def test_read_from_path(self, tmp_path) -> None:
        path = tmp_path / "in.vcf"
        path.write_text(HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\t.\n")
        header, records = read_vcf(path)
        assert len(records) == 1
        with pytest.raises(VcfFormatError) as excinfo:
            bad = tmp_path / "bad.vcf"
            bad.write_text("nope\n")
            read_vcf(bad)
        assert str(excinfo.value).startswith(f"{bad}:1:")


class TestNormalizeVcf:
    def test_snp_unchanged(self) -> None:
        record = norm("chr1\t4\t.\tT\tA\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (4, "T", ("A",))

    def test_common_prefix_trimmed(self) -> None:
        record = norm("chr1\t2\t.\tAA\tAG\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (3, "A", ("G",))

    def test_common_suffix_trimmed(self) -> None:
        record = norm("chr1\t4\t.\tTT\tAT\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (4, "T", ("A",))

    def test_insertion_left_shifted_to_minimal_pos(self) -> None:
        record = norm("chr1\t6\t.\tT\tTT\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (3, "A", ("AT",))

    def test_deletion_left_shifted_to_minimal_pos(self) -> None:
        record = norm("chr1\t5\t.\tTT\tT\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (3, "AT", ("A",))

    def test_left_shift_stops_at_pos_one(self) -> None:
        record = norm("chr1\t3\t.\tA\tAA\t.\t.\t.", {"chr1": "AAAA"})
        assert (record.pos, record.ref, record.alts) == (1, "A", ("AA",))

    def test_multiallelic_trimmed_as_a_group(self) -> None:
        record = norm("chr1\t1\t.\tAAA\tAAA,AA\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (1, "AA", ("AA", "A"))

    def test_alt_order_and_other_columns_preserved(self) -> None:
        record = norm(
            "chr1\t2\trs9\tAA\tAG,AC\t50\tq10\tDP=3",
            {"chr1": "AAATTTAAA"},
        )
        assert record.alts == ("G", "C")
        assert record.id == "rs9"
        assert record.qual == "50"
        assert record.filter == "q10"
        assert record.info == "DP=3"

    @pytest.mark.parametrize("alt", ["<DEL>", "A]chr2:9]", "*", "."])
    def test_special_allele_records_not_trimmed(self, alt: str) -> None:
        record = norm(f"chr1\t1\t.\tAAA\t{alt}\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert (record.pos, record.ref, record.alts) == (1, "AAA", (alt,))

    def test_special_allele_records_still_validate_ref(self) -> None:
        with pytest.raises(ReferenceMismatchError):
            norm("chr1\t1\t.\tACA\t<DEL>\t.\t.\t.", {"chr1": "AAATTTAAA"})

    def test_record_order_preserved(self) -> None:
        _, records = parse(
            HEADER
            + "\nchr1\t6\t.\tT\tTT\t.\t.\t.\nchr1\t4\t.\tT\tA\t.\t.\t.\n"
        )
        normalized = normalize_vcf(records, {"chr1": "AAATTTAAA"})
        assert [r.pos for r in normalized] == [3, 4]
        assert [r.ref for r in normalized] == ["A", "T"]

    def test_ref_mismatch_reports_expected_and_actual(self) -> None:
        with pytest.raises(ReferenceMismatchError) as excinfo:
            norm("chr1\t2\t.\tAC\tA\t.\t.\t.", {"chr1": "AAATTTAAA"})
        exc = excinfo.value
        assert exc.chrom == "chr1"
        assert exc.pos == 2
        assert exc.expected == "AA"
        assert exc.actual == "AC"
        message = str(exc)
        assert "chr1" in message and "2" in message
        assert "AA" in message and "AC" in message

    def test_missing_chrom(self) -> None:
        with pytest.raises(ReferenceMismatchError) as excinfo:
            norm("chr2\t1\t.\tA\tG\t.\t.\t.", {"chr1": "AAATTTAAA"})
        assert excinfo.value.chrom == "chr2"
        assert excinfo.value.pos == 1
        assert "chr2" in str(excinfo.value)

    def test_duplicate_chrom_in_reference(self) -> None:
        reference = [
            SequenceRecord("chr1", "AAATTTAAA"),
            SequenceRecord("chr1", "AAATTTAAA"),
        ]
        with pytest.raises(ReferenceMismatchError) as excinfo:
            norm("chr1\t1\t.\tA\tG\t.\t.\t.", reference)
        assert "chr1" in str(excinfo.value)

    def test_out_of_range_record_is_format_error(self) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            norm("chr1\t9\t.\tAA\tA\t.\t.\t.", {"chr1": "AAATTTAAA"})
        message = str(excinfo.value)
        assert message.startswith("test.vcf:2:")
        assert "chr1" in message and "9" in message

    def test_pos_beyond_reference_is_format_error(self) -> None:
        with pytest.raises(VcfFormatError):
            norm("chr1\t10\t.\tA\tG\t.\t.\t.", {"chr1": "AAATTTAAA"})

    def test_reference_as_sequence_records(self) -> None:
        record = norm(
            "chr1\t6\t.\tT\tTT\t.\t.\t.",
            [SequenceRecord("chr1", "AAATTTAAA")],
        )
        assert (record.pos, record.ref, record.alts) == (3, "A", ("AT",))


class TestWriteVcf:
    def test_round_trip_is_byte_identical(self) -> None:
        text = (
            "##fileformat=VCFv4.2\n"
            "##source=round\n"
            + FULL_HEADER
            + "\n"
            + "chr1\t3\trs1\tA\tG\t50\tPASS\tDP=4\tGT:DP\t0/1:4\t0/0:4\n"
            + "chr1\t1\t.\tAAA\t<DEL>\t.\t.\t.\tGT\t0/1\t0/0\n"
        )
        header, records = parse(text)
        out = StringIO()
        write_vcf(header, records, out)
        assert out.getvalue() == text

    def test_empty_records_still_write_full_header(self) -> None:
        header, records = parse("##a=b\n" + HEADER + "\n")
        out = StringIO()
        write_vcf(header, records, out)
        assert out.getvalue() == "##a=b\n" + HEADER + "\n"

    def test_exactly_one_trailing_newline(self) -> None:
        header, records = parse(HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\t.\n")
        out = StringIO()
        write_vcf(header, records, out)
        assert out.getvalue().endswith("\n")
        assert not out.getvalue().endswith("\n\n")

    def test_column_count_checked_against_header(self) -> None:
        header = VcfHeader((), ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
        record = VcfRecord("chr1", 1, ref="A", alts=("G",), format="GT", samples=("0/1",))
        with pytest.raises(ValueError):
            write_vcf(header, [record], StringIO())

    def test_invalid_header_columns_rejected(self) -> None:
        header = VcfHeader((), ("#CHROM", "POS"))
        with pytest.raises(ValueError):
            write_vcf(header, [], StringIO())

    def test_write_to_path(self, tmp_path) -> None:
        header, records = parse(HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\t.\n")
        path = tmp_path / "out.vcf"
        write_vcf(header, records, path)
        assert path.read_text() == HEADER + "\nchr1\t1\t.\tA\tG\t.\t.\t.\n"
