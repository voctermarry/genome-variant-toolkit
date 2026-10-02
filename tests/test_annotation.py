"""Unit tests for :mod:`genome_variant.annotation`."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.annotation import AnnotationFormatError, annotate_vcf
from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import (
    ReferenceMismatchError,
    VcfFile,
    VcfHeader,
    VcfRecord,
    read_vcf,
)

HEADER = "##fileformat=VCFv4.2\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
GVANN_DECLARATION = '##INFO=<ID=GVANN,Number=.,Type=String,'


def gff(text: str) -> StringIO:
    return StringIO(text)


def reference(sequence: str, identifier: str = "chr1") -> list[SequenceRecord]:
    return [SequenceRecord(identifier, sequence)]


def document(text: str) -> VcfFile:
    return read_vcf(StringIO(text))


def record_at(pos: int, ref: str, alt: str, info: str = ".") -> VcfRecord:
    return VcfRecord(
        chrom="chr1",
        pos=pos,
        id=".",
        ref=ref,
        alt=tuple(alt.split(",")),
        qual=".",
        filter=".",
        info=info,
    )


def annotate_text(vcf_text: str, sequence: str, gff_text: str) -> VcfFile:
    return annotate_vcf(document(vcf_text), reference(sequence), gff(gff_text))


def gvann_of(document_: VcfFile, index: int = 0) -> str:
    info = document_.records[index].info
    assert info.startswith("GVANN=")
    return info[len("GVANN="):]


class TestConsequences:
    """ATG AAA TTT GGG CCC ... — five codons on the plus strand."""

    SEQUENCE = "ATGAAATTTGGGCCC"
    CDS = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"

    def test_start_lost(self) -> None:
        text = HEADER + COLUMNS + "chr1\t2\t.\tT\tC\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "C|START_LOST|HIGH|t1|2|ATG>ACG|M>T"

    def test_start_codon_change_is_start_lost(self) -> None:
        # ATG -> GTG (V) under the standard code loses the start.
        text = HEADER + COLUMNS + "chr1\t1\t.\tA\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "G|START_LOST|HIGH|t1|1|ATG>GTG|M>V"

    def test_stop_gained(self) -> None:
        # codon AAA at coding positions 4-6; first base A->T gives TAA.
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "T|STOP_GAINED|HIGH|t1|4|AAA>TAA|K>*"

    def test_synonymous(self) -> None:
        # TTT -> TTC at positions 7-9, both F.
        text = HEADER + COLUMNS + "chr1\t9\t.\tT\tC\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "C|SYNONYMOUS|LOW|t1|9|TTT>TTC|F>F"

    def test_missense(self) -> None:
        # GGG -> AGG at positions 10-12, G -> R.
        text = HEADER + COLUMNS + "chr1\t10\t.\tG\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|MISSENSE|MODERATE|t1|10|GGG>AGG|G>R"

    def test_stop_lost(self) -> None:
        sequence = "ATGAAATAGGGG"  # coding: M K * G
        cds = "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
        # T A G at positions 7-9; T->C gives CAG (Q): STOP_LOST.
        text = HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == "C|STOP_LOST|HIGH|t1|7|TAG>CAG|*>Q"

    def test_stop_to_stop_is_synonymous(self) -> None:
        sequence = "ATGAAATAA"  # TAA at 7-9; third base A->G gives TAG.
        cds = "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
        text = HEADER + COLUMNS + "chr1\t9\t.\tA\tG\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == "G|SYNONYMOUS|LOW|t1|9|TAA>TAG|*>*"


class TestNonCodingAndUnsupported:
    SEQUENCE = "ATGAAATTTGGGCCCAAAAA"
    CDS = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"

    def test_snv_outside_cds(self) -> None:
        text = HEADER + COLUMNS + "chr1\t16\t.\tA\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "G|NON_CODING|MODIFIER|.|.|.|."

    def test_insertion_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tAT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "AT|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_deletion_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tAA\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_ambiguous_alt_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tR\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "R|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_symbolic_alt_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\t<DEL>\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "<DEL>|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_uncalled_missing_alt_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\t.\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == ".|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_same_base_ref_alt_unsupported(self) -> None:
        # A/A is not a substitution even though both symbols are ACGT.
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|UNSUPPORTED|MODIFIER|.|.|.|."


class TestMultiple:
    SEQUENCE = "ATGAAATTTGGGCCC"

    def test_multiple_alts_keep_order(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        text = HEADER + COLUMNS + "chr1\t2\t.\tT\tC,G,A\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        annotations = gvann_of(result).split(",")
        assert [item.split("|", 1)[0] for item in annotations] == ["C", "G", "A"]
        assert annotations[0] == "C|START_LOST|HIGH|t1|2|ATG>ACG|M>T"
        # ATG -> AGG (R) also loses the start.
        assert annotations[1] == "G|START_LOST|HIGH|t1|2|ATG>AGG|M>R"
        assert annotations[2] == "A|START_LOST|HIGH|t1|2|ATG>AAG|M>K"

    def test_multiple_transcripts_lexicographic(self) -> None:
        cds = (
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t_b\n"
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t_a\n"
        )
        text = HEADER + COLUMNS + "chr1\t2\t.\tT\tC\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert [item.split("|")[3] for item in items] == ["t_a", "t_b"]

    def test_mix_of_coding_and_non_coding_alt(self) -> None:
        # First ALT is an in-CDS SNV, second is an indel.
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        text = HEADER + COLUMNS + "chr1\t2\t.\tT\tC,TC\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert items[0] == "C|START_LOST|HIGH|t1|2|ATG>ACG|M>T"
        assert items[1] == "TC|UNSUPPORTED|MODIFIER|.|.|.|."


class TestMinusStrand:
    # Reverse complement of the 9 bases is ATG AA A TAA -> ATGAAATAA.
    SEQUENCE = "TTATTTCAT"
    CDS = "chr1\tx\tCDS\t1\t9\t.\t-\t0\tParent=t1\n"

    def test_codon_is_reverse_complemented_start_lost(self) -> None:
        # Genomic position 9 is coding position 1; T->G makes coding C,
        # turning the start ATG into CTG.
        text = HEADER + COLUMNS + "chr1\t9\t.\tT\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "G|START_LOST|HIGH|t1|1|ATG>CTG|M>L"

    def test_synonymous_stop_on_minus(self) -> None:
        # Genomic position 1 is coding position 9; T->C gives coding G,
        # TAA -> TAG (stop retained).
        text = HEADER + COLUMNS + "chr1\t1\t.\tT\tC\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "C|SYNONYMOUS|LOW|t1|9|TAA>TAG|*>*"

    def test_stop_gained_on_minus(self) -> None:
        # Genomic position 6 is coding position 4; T->A gives coding T,
        # AAA -> TAA.
        text = HEADER + COLUMNS + "chr1\t6\t.\tT\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|STOP_GAINED|HIGH|t1|4|AAA>TAA|K>*"

    def test_split_minus_cds_phase(self) -> None:
        # Two minus-strand fragments: 5-14 (phase 0, ten bases) and
        # 1-3 (phase 2, three bases); coding positions continue across
        # the junction.
        sequence = "TTATTTCATGGGGG"  # 14 bases
        cds = (
            "chr1\tx\tCDS\t5\t14\t.\t-\t0\tParent=t1\n"
            "chr1\tx\tCDS\t1\t3\t.\t-\t2\tParent=t1\n"
        )
        text = HEADER + COLUMNS + "chr1\t2\t.\tT\tC\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == "C|MISSENSE|MODERATE|t1|12|ATA>ATG|I>M"


class TestHeaderAndInfo:
    SEQUENCE = "ATGAAATTTGGGCCC"
    CDS = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"

    def test_unique_gvann_declaration_appended(self) -> None:
        text = HEADER + COLUMNS
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        declarations = [
            line for line in result.header.meta_lines if line.startswith("##INFO=<ID=GVANN")
        ]
        assert len(declarations) == 1
        assert result.header.meta_lines[-1] == declarations[0]
        assert "ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|AA_CHANGE" in (
            declarations[0]
        )

    def test_existing_header_declaration_rejected(self) -> None:
        text = (
            HEADER
            + '##INFO=<ID=GVANN,Number=.,Type=String,Description="x">\n'
            + COLUMNS
        )
        with pytest.raises(AnnotationFormatError):
            annotate_text(text, self.SEQUENCE, self.CDS)

    def test_existing_info_item_rejected(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\tGVANN=x\n"
        with pytest.raises(AnnotationFormatError) as excinfo:
            annotate_text(text, self.SEQUENCE, self.CDS)
        assert "GVANN" in str(excinfo.value)

    def test_info_appended_to_existing_items(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\tDP=5;AF=0.5\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        info = result.records[0].info
        assert info.startswith("DP=5;AF=0.5;GVANN=")

    def test_missing_info_starts_with_gvann(self) -> None:
        text = HEADER + COLUMNS + "chr1\t16\t.\tA\tG\t.\t.\t.\n"
        result = annotate_text(
            text, self.SEQUENCE + "AAA", self.CDS
        )
        assert result.records[0].info.startswith("GVANN=")

    def test_samples_and_fields_preserved(self) -> None:
        header = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        )
        row = "chr1\t4\t.\tA\tT\t99\tPASS\tDP=5\tGT:DP\t0/1:9\t1/1:2\n"
        result = annotate_text(HEADER + header + row, self.SEQUENCE, self.CDS)
        record = result.records[0]
        assert record.format_text == "GT:DP"
        assert record.sample_text == ("0/1:9", "1/1:2")
        assert record.qual == "99"
        assert record.filter == "PASS"
        assert record.alt == ("T",)
        assert record.info.startswith("DP=5;GVANN=")
        assert record.pos == 4 and record.ref == "A"

    def test_records_keep_order(self) -> None:
        rows = (
            "chr1\t4\t.\tA\tT\t.\t.\t.\n"
            "chr1\t7\t.\tT\tA\t.\t.\t.\n"
            "chr1\t10\t.\tG\tA\t.\t.\t.\n"
        )
        result = annotate_text(HEADER + COLUMNS + rows, self.SEQUENCE, self.CDS)
        assert [record.pos for record in result.records] == [4, 7, 10]


class TestGffValidation:
    SEQUENCE = "ATGAAATTTGGGCCC"

    def annotate(self, gff_text: str, vcf_text: str | None = None) -> VcfFile:
        if vcf_text is None:
            vcf_text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        return annotate_text(vcf_text, self.SEQUENCE, gff_text)

    def test_wrong_column_count(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t9\t.\t+\t0\n")
        message = str(excinfo.value)
        assert "<stream>:1:" in message
        assert "9 tab-separated" in message

    def test_bad_coordinates(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\tzero\t9\t.\t+\t0\tParent=t1\n")
        assert ":1:" in str(excinfo.value)

    def test_non_positive_coordinates(self) -> None:
        with pytest.raises(AnnotationFormatError):
            self.annotate("chr1\tx\tCDS\t0\t9\t.\t+\t0\tParent=t1\n")

    def test_start_greater_than_end(self) -> None:
        with pytest.raises(AnnotationFormatError):
            self.annotate("chr1\tx\tCDS\t9\t4\t.\t+\t0\tParent=t1\n")

    def test_bad_strand(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t9\t.\t.\t0\tParent=t1\n")
        assert "strand" in str(excinfo.value)

    def test_bad_phase(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t9\t.\t+\t3\tParent=t1\n")
        assert "phase" in str(excinfo.value)

    def test_missing_parent(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t9\t.\t+\t0\tID=c1\n")
        assert "Parent" in str(excinfo.value) and ":1:" in str(excinfo.value)

    def test_empty_parent_value(self) -> None:
        with pytest.raises(AnnotationFormatError):
            self.annotate("chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=\n")

    def test_fragments_spanning_sequences(self) -> None:
        cds = (
            "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr2\tx\tCDS\t1\t6\t.\t+\t0\tParent=t1\n"
        )
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate(cds)
        message = str(excinfo.value)
        assert ":2:" in message and "multiple sequences" in message

    def test_contradictory_phase(self) -> None:
        cds = (
            "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr1\tx\tCDS\t10\t15\t.\t+\t1\tParent=t1\n"
        )
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate(cds)
        message = str(excinfo.value)
        assert ":2:" in message and "reading frame" in message

    def test_non_cds_features_ignored(self) -> None:
        cds = (
            "chr1\tx\tgene\t1\t15\t.\t+\t.\tID=g1\n"
            "chr1\tx\texon\t1\t15\t.\t+\t.\tParent=t1\n"
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        )
        result = self.annotate(cds)
        assert "t1" in gvann_of(result)

    def test_comments_and_pragmas_ignored(self) -> None:
        cds = (
            "##gff-version 3\n"
            "# a comment\n"
            "\n"
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        )
        result = self.annotate(cds)
        assert "t1" in gvann_of(result)

    def test_coordinates_past_reference_end(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t40\t.\t+\t0\tParent=t1\n")
        message = str(excinfo.value)
        assert ":1:" in message and "past the end" in message


class TestReferenceErrors:
    CDS = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"

    def test_missing_chrom(self) -> None:
        text = HEADER + COLUMNS + "chr9\t4\t.\tA\tT\t.\t.\t.\n"
        with pytest.raises(ReferenceMismatchError):
            annotate_text(text, "ATGAAATTTGGGCCC", self.CDS)

    def test_duplicate_chrom(self) -> None:
        records = [
            SequenceRecord("chr1", "ATGAAATTTGGGCCC"),
            SequenceRecord("chr1", "ATGAAATTTGGGCCC"),
        ]
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        with pytest.raises(ReferenceMismatchError):
            annotate_vcf(document(text), records, gff(self.CDS))

    def test_ref_mismatch(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tC\tT\t.\t.\t.\n"
        with pytest.raises(ReferenceMismatchError):
            annotate_text(text, "ATGAAATTTGGGCCC", self.CDS)

    def test_out_of_bounds_is_vcf_error(self) -> None:
        from genome_variant.vcf import VcfFormatError

        text = HEADER + COLUMNS + "chr1\t99\t.\tA\tT\t.\t.\t.\n"
        with pytest.raises(VcfFormatError):
            annotate_text(text, "ATGAAATTTGGGCCC", self.CDS)


class TestBuiltDirectly:
    def test_record_without_line_number(self) -> None:
        doc = VcfFile(
            VcfHeader(),
            (record_at(4, "A", "T"),),
        )
        result = annotate_vcf(
            doc,
            reference("ATGAAATTTGGGCCC"),
            gff("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"),
        )
        assert result.records[0].info.startswith("GVANN=")
