"""Tests for reference-aware VCF reading and normalization."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import (
    ReferenceMismatchError,
    VcfFile,
    VcfFormatError,
    VcfHeader,
    VcfRecord,
    normalize_record,
    normalize_vcf,
    read_vcf,
    render_vcf,
)

HEADER = "##fileformat=VCFv4.2\n##source=test\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
COLUMNS_SAMPLES = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
)


def parse(text: str) -> VcfFile:
    return read_vcf(StringIO(text))


def normalize(text: str, reference: dict[str, str] | list[SequenceRecord]) -> str:
    document = parse(text)
    if isinstance(reference, dict):
        normalized = normalize_vcf(document, reference)
    else:
        normalized = normalize_vcf(document, reference)
    return render_vcf(normalized)


class TestParsing:
    def test_meta_lines_kept_verbatim_and_in_order(self) -> None:
        document = parse(HEADER + COLUMNS)
        assert document.header.meta_lines == ("##fileformat=VCFv4.2", "##source=test")
        assert document.header.samples == ()
        assert document.records == ()

    def test_records_emitted_in_input_order(self) -> None:
        rows = (
            "chr1\t2\t.\tA\tG\t.\t.\t.\n"
            "chr2\t9\t.\tC\tT\t.\t.\t.\n"
            "chr1\t1\t.\tG\tA\t.\t.\t.\n"
        )
        document = parse(HEADER + COLUMNS + rows)
        assert [record.chrom for record in document.records] == ["chr1", "chr2", "chr1"]
        assert [record.pos for record in document.records] == [2, 9, 1]

    def test_plain_alleles_are_uppercased(self) -> None:
        document = parse(HEADER + COLUMNS + "chr1\t2\t.\taC\tgT,t\t.\t.\t.\n")
        record = document.records[0]
        assert record.ref == "AC"
        assert record.alt == ("GT", "T")

    def test_iupac_symbols_allowed_in_plain_alleles(self) -> None:
        document = parse(HEADER + COLUMNS + "chr1\t2\t.\tAR\tY,N\t.\t.\t.\n")
        record = document.records[0]
        assert record.ref == "AR"
        assert record.alt == ("Y", "N")

    def test_sample_columns_kept_verbatim(self) -> None:
        row = "chr1\t2\t.\tA\tG\t50\tPASS\tNS=1\tGT:DP\t0/1:12\t./.:3\n"
        document = parse(HEADER + COLUMNS_SAMPLES + row)
        record = document.records[0]
        assert record.format_text == "GT:DP"
        assert record.sample_text == ("0/1:12", "./.:3")
        out = render_vcf(document)
        assert out == HEADER + COLUMNS_SAMPLES + row

    @pytest.mark.parametrize(
        "text",
        [
            # No #CHROM at all
            "##fileformat=VCFv4.2\n",
            # Two #CHROM lines
            HEADER + COLUMNS + COLUMNS,
            # Wrong fixed column names
            HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tEXTRA\n",
            # Tail without FORMAT
            HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\ts1\n",
            # FORMAT with no sample
            HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n",
            # Empty sample name
            HEADER
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\t\n",
            # Data before header
            "chr1\t1\t.\tA\tG\t.\t.\t.\n" + COLUMNS,
        ],
    )
    def test_bad_header_raises_format_error(self, text: str) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(text)
        message = str(excinfo.value)
        assert "<stream>" in message
        assert any(f":{n}:" in message for n in range(1, 6))

    @pytest.mark.parametrize(
        "row",
        [
            "chr1\t2\t.\tA\tG\t.\t.\n",  # 7 columns
            "chr1\t2\t.\tA\tG\t.\t.\t.\textra\n",  # 9 columns, no samples
            "chr1\t0\t.\tA\tG\t.\t.\t.\n",  # non-positive POS
            "chr1\t-2\t.\tA\tG\t.\t.\t.\n",
            "chr1\tx\t.\tA\tG\t.\t.\t.\n",
            "chr1\t2\t.\t\tG\t.\t.\t.\n",  # empty REF
            "chr1\t2\t.\tAZ\tG\t.\t.\t.\n",  # illegal REF base
            "chr1\t2\t.\tA\tGZ\t.\t.\t.\n",  # illegal plain ALT base
            "chr1\t2\t.\tA\t\t.\t.\t.\n",  # empty ALT column
            "chr1\t2\t.\tA\tG,\t.\t.\t.\n",  # empty ALT allele
        ],
    )
    def test_bad_rows_raise_format_error_with_source_and_line(self, row: str) -> None:
        with pytest.raises(VcfFormatError) as excinfo:
            parse(HEADER + COLUMNS + row + "\n")
        message = str(excinfo.value)
        assert message.startswith("<stream>:4:")

    def test_sample_column_count_must_match_header(self) -> None:
        row = "chr1\t2\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"
        with pytest.raises(VcfFormatError, match="<stream>:4:"):
            parse(HEADER + COLUMNS_SAMPLES + row)

    def test_blank_line_is_not_accepted(self) -> None:
        with pytest.raises(VcfFormatError):
            parse(HEADER + COLUMNS + "\nchr1\t1\t.\tA\tG\t.\t.\t.\n")

    def test_render_ends_with_single_newline_and_tabs(self) -> None:
        row = "chr1\t2\t.\tA\tG\t.\t.\t.\n"
        out = render_vcf(parse(HEADER + COLUMNS + row))
        assert out.endswith("\n")
        assert not out.endswith("\n\n")
        assert "\r" not in out
        assert out == HEADER + COLUMNS + row

    def test_empty_vcf_still_writes_full_header(self) -> None:
        assert render_vcf(parse(HEADER + COLUMNS)) == HEADER + COLUMNS

    def test_read_from_path(self, tmp_path) -> None:
        path = tmp_path / "in.vcf"
        path.write_text(HEADER + COLUMNS)
        document = read_vcf(path)
        assert document.source == str(path)
        assert document.header.meta_lines[0] == "##fileformat=VCFv4.2"


class TestNormalization:
    REFERENCE = {"chr1": "ACGAACAGG"}

    def test_strip_common_suffix_and_prefix(self) -> None:
        # positions 1-based: 1 A 2 C 3 G 4 A 5 A 6 C 7 A 8 G 9 G
        # chr1:6 C, REF CAG (positions 6-8), ALT CTG (->C...G suffix G,
        # prefix C) so it collapses to A>T at POS 7.
        out = normalize(
            HEADER + COLUMNS + "chr1\t6\t.\tCAG\tCTG\t.\t.\t.\n",
            self.REFERENCE,
        )
        assert "chr1\t7\t.\tA\tT\t.\t.\t.\n" in out

    def test_at_least_one_base_retained_per_allele(self) -> None:
        # REF ACG vs ALT A at POS1 (reference ACG...): shared prefix A
        # cannot be stripped because ALT would vanish.
        out = normalize(
            HEADER + COLUMNS + "chr1\t1\t.\tACG\tA\t.\t.\t.\n",
            {"chr1": "ACGT"},
        )
        assert "chr1\t1\t.\tACG\tA\t.\t.\t.\n" in out

    def test_left_shift_deletion_across_repeat(self) -> None:
        # reference ACGGT: deleting one G is written at POS4 as GT>T and
        # walks left while the last anchor bases agree, stopping at POS2
        # (CG>C) where the preceding A would describe a different
        # haplotype.  Both forms spell the ACGT haplotype.
        out = normalize(
            HEADER + COLUMNS + "chr1\t4\t.\tGT\tT\t.\t.\t.\n",
            {"chr1": "ACGGT"},
        )
        assert "chr1\t2\t.\tCG\tC\t.\t.\t.\n" in out

    def test_left_shift_reaches_chromosome_start(self) -> None:
        # reference AGGGT: the deletion walks all the way to POS1,
        # giving AG>A (haplotype AGGT).
        out = normalize(
            HEADER + COLUMNS + "chr1\t4\t.\tGT\tT\t.\t.\t.\n",
            {"chr1": "AGGGT"},
        )
        assert "chr1\t1\t.\tAG\tA\t.\t.\t.\n" in out

    def test_left_shift_insertion_across_repeat(self) -> None:
        # reference CAA: an extra A written at POS3 (A>AA, haplotype
        # CAAA) walks left across the homopolymer to POS1 (C>CA).
        out = normalize(
            HEADER + COLUMNS + "chr1\t3\t.\tA\tAA\t.\t.\t.\n",
            {"chr1": "CAA"},
        )
        assert "chr1\t1\t.\tC\tCA\t.\t.\t.\n" in out

    def test_no_shift_when_allele_cannot_move_left(self) -> None:
        # CACG: CG>C (haplotype CAC) cannot shift; the extended window
        # AC>A would describe CAG instead, so it is rejected.
        out = normalize(
            HEADER + COLUMNS + "chr1\t3\t.\tCG\tC\t.\t.\t.\n",
            {"chr1": "CACG"},
        )
        assert "chr1\t3\t.\tCG\tC\t.\t.\t.\n" in out

    def test_shift_stops_at_chromosome_start(self) -> None:
        # reference AAA: POS3 A>AA (insertion) should shift to POS2 A>AA,
        # then POS1 A>AA and stop there.
        out = normalize(
            HEADER + COLUMNS + "chr1\t3\t.\tA\tAA\t.\t.\t.\n",
            {"chr1": "AAA"},
        )
        assert "chr1\t1\t.\tA\tAA\t.\t.\t.\n" in out

    def test_equal_length_alleles_are_not_shifted(self) -> None:
        out = normalize(
            HEADER + COLUMNS + "chr1\t2\t.\tC\tT\t.\t.\t.\n",
            {"chr1": "ACG"},
        )
        assert "chr1\t2\t.\tC\tT\t.\t.\t.\n" in out

    @pytest.mark.parametrize(
        "alt",
        [
            "<DEL>",
            "G<DEL>",
            "]chr2:123]G",
            "G[chr2:123[",
            "G.",
            ".G",
            "*",
            "G,*",
            ".",
            "T,.",
        ],
    )
    def test_special_alts_not_trimmed_or_shifted_but_ref_checked(self, alt: str) -> None:
        reference = {"chr1": "ACGTAC"}
        row = f"chr1\t2\t.\tCG\t{alt}\t.\t.\t.\n"
        out = normalize(HEADER + COLUMNS + row, reference)
        assert row in out

    def test_special_alt_with_bad_ref_still_raises(self) -> None:
        with pytest.raises(ReferenceMismatchError) as excinfo:
            normalize(
                HEADER + COLUMNS + "chr1\t2\t.\tCG\t<DEL>\t.\t.\t.\n",
                {"chr1": "TTTTTT"},
            )
        message = str(excinfo.value)
        assert "chr1" in message and "POS 2" in message
        assert "TT" in message and "CG" in message

    def test_alt_order_and_other_columns_preserved(self) -> None:
        row = "chr1\t6\t.\tCAG\tCTG,CCC\t99\tPASS\tNS=2\tGT\t0/1\t0/2\n"
        out = normalize(
            HEADER + COLUMNS_SAMPLES + row, {"chr1": "ACGAACAGG"}
        )
        line = [line for line in out.splitlines() if line.startswith("chr1")][0]
        fields = line.split("\t")
        assert fields[1] == "7"
        assert fields[3:6] == ["AG", "TG,CC", "99"]
        assert fields[6:9] == ["PASS", "NS=2", "GT"]
        assert fields[9:] == ["0/1", "0/2"]

    def test_record_order_preserved(self) -> None:
        rows = (
            "chr1\t4\t.\tGT\tT\t.\t.\t.\n"
            "chr1\t2\t.\tC\tT\t.\t.\t.\n"
        )
        out = normalize(HEADER + COLUMNS + rows, {"chr1": "ACGGT"})
        data = [line for line in out.splitlines() if not line.startswith("#")]
        assert data[0] == "chr1\t2\t.\tCG\tC\t.\t.\t."
        assert data[1] == "chr1\t2\t.\tC\tT\t.\t.\t."

    def test_idempotent(self) -> None:
        reference = {"chr1": "ACGGT"}
        first = normalize(
            HEADER + COLUMNS + "chr1\t4\t.\tGT\tT\t.\t.\t.\n", reference
        )
        second = normalize(first, reference)
        assert first == second


class TestReferenceErrors:
    def test_missing_chrom(self) -> None:
        document = parse(HEADER + COLUMNS + "chr9\t1\t.\tA\tG\t.\t.\t.\n")
        with pytest.raises(ReferenceMismatchError, match="chr9"):
            normalize_vcf(document, {"chr1": "ACGT"})

    def test_duplicate_chrom_in_reference_records(self) -> None:
        document = parse(HEADER + COLUMNS + "chr1\t1\t.\tA\tG\t.\t.\t.\n")
        reference = [
            SequenceRecord("chr1", "ACGT"),
            SequenceRecord("chr1", "TTTT"),
        ]
        with pytest.raises(ReferenceMismatchError, match="chr1"):
            normalize_vcf(document, reference)

    def test_duplicate_chrom_allowed_through_mapping(self) -> None:
        # A plain mapping cannot represent duplicates; it is used as-is.
        document = parse(HEADER + COLUMNS + "chr1\t1\t.\tA\tG\t.\t.\t.\n")
        result = normalize_vcf(document, {"chr1": "ACGT"})
        assert result.records[0].alt == ("G",)

    def test_ref_mismatch_message_has_chrom_pos_expected_actual(self) -> None:
        document = parse(HEADER + COLUMNS + "chr1\t2\t.\tC\tT\t.\t.\t.\n")
        with pytest.raises(ReferenceMismatchError) as excinfo:
            normalize_vcf(document, {"chr1": "TTTT"})
        message = str(excinfo.value)
        assert "chr1" in message
        assert "POS 2" in message
        assert "'T'" in message  # expected
        assert "'C'" in message  # actual

    def test_out_of_bounds_is_format_error_with_source_and_line(self) -> None:
        document = parse(HEADER + COLUMNS + "chr1\t9\t.\tAA\tG\t.\t.\t.\n")
        with pytest.raises(VcfFormatError) as excinfo:
            normalize_vcf(document, {"chr1": "ACGT"})
        message = str(excinfo.value)
        assert "<stream>:4:" in message
        assert "chr1" in message and "POS 9" in message


class TestRecordApi:
    def test_normalize_record_accepts_mapping_and_records(self) -> None:
        record = VcfRecord("chr1", 4, ".", "GT", ("T",), ".", ".", ".")
        shifted = normalize_record(record, {"chr1": "ACGGT"})
        assert (shifted.pos, shifted.ref, shifted.alt) == (2, "CG", ("C",))
        shifted2 = normalize_record(
            record, [SequenceRecord("chr1", "ACGGT")]
        )
        assert shifted == shifted2

    def test_byte_stable_across_runs(self) -> None:
        text = (
            HEADER
            + COLUMNS
            + "chr1\t4\t.\tGT\tT\t.\t.\t.\n"
            + "chr1\t2\t.\tC\tT\t.\t.\t.\n"
        )
        first = render_vcf(normalize_vcf(parse(text), {"chr1": "ACGGT"}))
        second = render_vcf(normalize_vcf(parse(text), {"chr1": "ACGGT"}))
        assert first == second

    def test_header_column_names(self) -> None:
        no_samples = VcfHeader(("##a",), ())
        assert no_samples.column_names == (
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
        )
        with_samples = VcfHeader(("##a",), ("x", "y"))
        assert with_samples.column_names == (
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "x", "y",
        )
