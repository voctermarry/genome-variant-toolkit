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

    def test_complex_replacement_unsupported(self) -> None:
        # Reference at positions 3-6 is GAAA; GAAA -> GTT both deletes
        # bases and substitutes one, so it is not a pure indel.
        text = HEADER + COLUMNS + "chr1\t3\t.\tGAAA\tGTT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "GTT|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_same_length_mnv_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tAA\tTT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "TT|UNSUPPORTED|MODIFIER|.|.|.|."

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
        # First ALT is an in-CDS SNV, second is an ambiguous allele.
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        text = HEADER + COLUMNS + "chr1\t2\t.\tT\tC,R\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert items[0] == "C|START_LOST|HIGH|t1|2|ATG>ACG|M>T"
        assert items[1] == "R|UNSUPPORTED|MODIFIER|.|.|.|."


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


class TestIndelsPlusStrand:
    """ATG AAA TTT GGG CCC — five codons on chr1, one CDS fragment."""

    SEQUENCE = "ATGAAATTTGGGCCC"
    CDS = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"

    def test_inframe_deletion_three_bases(self) -> None:
        # Delete coding bases 4-6 (AAA): anchor base 3, right base 7.
        text = HEADER + COLUMNS + "chr1\t3\t.\tGAAA\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "G|INFRAME_DELETION|MODERATE|t1|4|.|."
        )

    def test_inframe_deletion_six_bases(self) -> None:
        # Delete coding bases 4-9, anchor base 3 (G), right base 10 (G).
        text = HEADER + COLUMNS + "chr1\t3\t.\tGAAATTT\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "G|INFRAME_DELETION|MODERATE|t1|4|.|."
        )

    def test_frameshift_deletion_one_base(self) -> None:
        # POS 4 REF AA ALT A deletes coding base 5.
        text = HEADER + COLUMNS + "chr1\t4\t.\tAA\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|FRAMESHIFT|HIGH|t1|5|.|."

    def test_frameshift_deletion_two_bases(self) -> None:
        # POS 4 REF AAA ALT A deletes coding bases 5 and 6.
        text = HEADER + COLUMNS + "chr1\t4\t.\tAAA\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|FRAMESHIFT|HIGH|t1|5|.|."

    def test_inframe_insertion_three_bases(self) -> None:
        # Insert AAA after coding base 3; CDS_POS is the right base 4.
        text = HEADER + COLUMNS + "chr1\t3\t.\tG\tGAAA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "GAAA|INFRAME_INSERTION|MODERATE|t1|4|.|."
        )

    def test_frameshift_insertion_one_base(self) -> None:
        # Insert T after coding base 6; the right coding base is 7.
        text = HEADER + COLUMNS + "chr1\t6\t.\tA\tAT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "AT|FRAMESHIFT|HIGH|t1|7|.|."

    def test_frameshift_insertion_four_bases(self) -> None:
        text = HEADER + COLUMNS + "chr1\t3\t.\tG\tGACGT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "GACGT|FRAMESHIFT|HIGH|t1|4|.|."

    def test_insertion_at_cds_end_unsupported(self) -> None:
        # Anchor is the last coding base 15, the right base 16 is outside.
        sequence = self.SEQUENCE + "AAA"
        text = HEADER + COLUMNS + "chr1\t15\t.\tC\tCAAA\t.\t.\t.\n"
        result = annotate_text(text, sequence, self.CDS)
        assert gvann_of(result) == (
            "CAAA|UNSUPPORTED|MODIFIER|t1|15|.|."
        )

    def test_deletion_running_past_cds_end_unsupported(self) -> None:
        # Delete coding bases 14-15 plus non-coding base 16.
        sequence = self.SEQUENCE + "AAA"
        text = HEADER + COLUMNS + "chr1\t13\t.\tCCCA\tC\t.\t.\t.\n"
        result = annotate_text(text, sequence, self.CDS)
        assert gvann_of(result) == "C|UNSUPPORTED|MODIFIER|t1|14|.|."

    def test_insertion_after_cds_is_non_coding(self) -> None:
        sequence = self.SEQUENCE + "AAA"
        text = HEADER + COLUMNS + "chr1\t16\t.\tA\tAT\t.\t.\t.\n"
        result = annotate_text(text, sequence, self.CDS)
        assert gvann_of(result) == "AT|NON_CODING|MODIFIER|.|.|.|."

    def test_deletion_after_cds_is_non_coding(self) -> None:
        sequence = self.SEQUENCE + "AAAAA"
        text = HEADER + COLUMNS + "chr1\t16\t.\tAA\tA\t.\t.\t.\n"
        result = annotate_text(text, sequence, self.CDS)
        assert gvann_of(result) == "A|NON_CODING|MODIFIER|.|.|.|."

    def test_deletion_starting_inside_cds_uses_first_deleted_base(self) -> None:
        # Delete bases 2-4, anchored at coding base 1.
        text = HEADER + COLUMNS + "chr1\t1\t.\tATGA\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A|INFRAME_DELETION|MODERATE|t1|2|.|."


class TestIndelFragments:
    """Two adjacent CDS fragments 1-9 and 10-15 on the plus strand."""

    SEQUENCE = "ATGAAATTTGGGCCC"
    CDS = (
        "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
        "chr1\tx\tCDS\t10\t15\t.\t+\t0\tParent=t1\n"
    )

    def test_deletion_crossing_fragment_junction_unsupported(self) -> None:
        # Delete bases 8-10: anchor 7 in fragment one, right base 11 in two.
        text = HEADER + COLUMNS + "chr1\t7\t.\tTTTG\tT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "T|UNSUPPORTED|MODIFIER|t1|8|.|."

    def test_insertion_crossing_fragment_junction_unsupported(self) -> None:
        # Boundary between coding base 9 (fragment one) and 10 (fragment two).
        text = HEADER + COLUMNS + "chr1\t9\t.\tT\tTAAA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "TAAA|UNSUPPORTED|MODIFIER|t1|10|.|."

    def test_inframe_deletion_within_second_fragment(self) -> None:
        text = HEADER + COLUMNS + "chr1\t10\t.\tGGGC\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "G|INFRAME_DELETION|MODERATE|t1|11|.|."
        )

    def test_inframe_insertion_within_second_fragment(self) -> None:
        text = HEADER + COLUMNS + "chr1\t11\t.\tG\tGAAA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "GAAA|INFRAME_INSERTION|MODERATE|t1|12|.|."
        )

    def test_deletion_spanning_intron_unsupported(self) -> None:
        # Fragments 1-9 and 16-18; delete bases 8-16 (mostly intronic).
        sequence = "ATGAAATTTGGGCCCAAAAA"
        cds = (
            "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr1\tx\tCDS\t16\t18\t.\t+\t0\tParent=t1\n"
        )
        text = HEADER + COLUMNS + "chr1\t7\t.\tTTTGGGCCCA\tT\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == "T|UNSUPPORTED|MODIFIER|t1|8|.|."

    def test_insertion_inside_intron_is_non_coding(self) -> None:
        sequence = "ATGAAATTTGGGCCCAAAAA"
        cds = (
            "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr1\tx\tCDS\t16\t18\t.\t+\t0\tParent=t1\n"
        )
        text = HEADER + COLUMNS + "chr1\t12\t.\tG\tGAAA\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == "GAAA|NON_CODING|MODIFIER|.|.|.|."


class TestIndelsMinusStrand:
    # Reverse complement of the 9 bases is ATG AAA TAA; genomic base 9 is
    # coding base 1 and coding order descends with genomic coordinates.
    SEQUENCE = "TTATTTCAT"
    CDS = "chr1\tx\tCDS\t1\t9\t.\t-\t0\tParent=t1\n"

    def test_inframe_deletion_three_bases(self) -> None:
        # Delete genomic bases 4-6 = coding bases 6,5,4; first affected 4.
        text = HEADER + COLUMNS + "chr1\t3\t.\tATTT\tA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "A|INFRAME_DELETION|MODERATE|t1|4|.|."
        )

    def test_frameshift_deletion_one_base(self) -> None:
        # Delete genomic base 6 (coding base 4), anchored at base 7.
        text = HEADER + COLUMNS + "chr1\t6\t.\tTC\tC\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "C|FRAMESHIFT|HIGH|t1|4|.|."

    def test_inframe_insertion_uses_right_coding_base(self) -> None:
        # Boundary between genomic bases 7 (coding 3) and 8 (coding 2):
        # the downstream side in coding (5'->3') direction is coding
        # position 3, so the insertion reports CDS_POS 3.
        text = HEADER + COLUMNS + "chr1\t7\t.\tC\tCAAA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == (
            "CAAA|INFRAME_INSERTION|MODERATE|t1|3|.|."
        )

    def test_frameshift_insertion_one_base(self) -> None:
        text = HEADER + COLUMNS + "chr1\t7\t.\tC\tCA\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "CA|FRAMESHIFT|HIGH|t1|3|.|."

    def test_split_minus_cds_inframe_deletion(self) -> None:
        # Fragments 5-14 (phase 0) and 1-3 (phase 2); delete genomic
        # bases 10-12 = coding bases 5,4,3, all within fragment 5-14.
        sequence = "TTATTTCATGGGGG"
        cds = (
            "chr1\tx\tCDS\t5\t14\t.\t-\t0\tParent=t1\n"
            "chr1\tx\tCDS\t1\t3\t.\t-\t2\tParent=t1\n"
        )
        text = HEADER + COLUMNS + "chr1\t9\t.\tTGGG\tT\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == (
            "T|INFRAME_DELETION|MODERATE|t1|3|.|."
        )

    def test_split_minus_cds_cross_junction_unsupported(self) -> None:
        # Delete genomic bases 4-5: base 5 is the last coding base of the
        # first fragment, base 4 is intronic and the anchors lie in
        # different fragments.
        sequence = "TTATTTCATGGGGG"
        cds = (
            "chr1\tx\tCDS\t5\t14\t.\t-\t0\tParent=t1\n"
            "chr1\tx\tCDS\t1\t3\t.\t-\t2\tParent=t1\n"
        )
        text = HEADER + COLUMNS + "chr1\t3\t.\tATT\tA\t.\t.\t.\n"
        result = annotate_text(text, sequence, cds)
        assert gvann_of(result) == "A|UNSUPPORTED|MODIFIER|t1|10|.|."


class TestIndelsMultiple:
    SEQUENCE = "ATGAAATTTGGGCCC"

    def test_multiple_indel_alts_keep_order(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        text = HEADER + COLUMNS + "chr1\t3\t.\tGAAA\tG,GAAAACGT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert [item.split("|", 1)[0] for item in items] == [
            "G",
            "GAAAACGT",
        ]
        assert items[0] == "G|INFRAME_DELETION|MODERATE|t1|4|.|."
        assert items[1] == "GAAAACGT|FRAMESHIFT|HIGH|t1|7|.|."

    def test_multiple_transcripts_lexicographic(self) -> None:
        cds = (
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t_b\n"
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t_a\n"
        )
        text = HEADER + COLUMNS + "chr1\t3\t.\tGAAA\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert [item.split("|")[3] for item in items] == ["t_a", "t_b"]
        assert items[0] == "G|INFRAME_DELETION|MODERATE|t_a|4|.|."
        assert items[1] == "G|INFRAME_DELETION|MODERATE|t_b|4|.|."

    def test_one_transcript_coding_one_unaffected(self) -> None:
        # t2's CDS starts at base 10, so the deletion at bases 4-6 only
        # touches t1.
        cds = (
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
            "chr1\tx\tCDS\t10\t15\t.\t+\t0\tParent=t2\n"
        )
        text = HEADER + COLUMNS + "chr1\t3\t.\tGAAA\tG\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert len(items) == 1
        assert items[0] == "G|INFRAME_DELETION|MODERATE|t1|4|.|."

    def test_mixed_snv_and_indel_alts(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT,AT\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, cds)
        items = gvann_of(result).split(",")
        assert items[0] == "T|STOP_GAINED|HIGH|t1|4|AAA>TAA|K>*"
        assert items[1] == "AT|FRAMESHIFT|HIGH|t1|5|.|."


class TestIndelUnsupportedInputs:
    SEQUENCE = "ATGAAATTTGGGCCC"
    CDS = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=t1\n"

    def test_inserted_ambiguous_base_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t3\t.\tG\tGR\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "GR|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_ambiguous_two_base_insertion_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tAN\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "AN|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_symbolic_alt_still_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\t<INS>\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "<INS>|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_spanning_alt_still_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\t*\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "*|UNSUPPORTED|MODIFIER|.|.|.|."

    def test_breakend_alt_still_unsupported(self) -> None:
        text = HEADER + COLUMNS + "chr1\t4\t.\tA\tA]chr2:9]\t.\t.\t.\n"
        result = annotate_text(text, self.SEQUENCE, self.CDS)
        assert gvann_of(result) == "A]chr2:9]|UNSUPPORTED|MODIFIER|.|.|.|."


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


class TestPercentDecoding:
    """GFF3 attribute values: %HH escapes are UTF-8 bytes."""

    SEQUENCE = "ATGAAATTTGGGCCC"

    def annotate(self, gff_text: str, vcf_text: str | None = None) -> VcfFile:
        if vcf_text is None:
            vcf_text = HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        return annotate_text(vcf_text, self.SEQUENCE, gff_text)

    def test_utf8_escape_decodes_to_characters(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%E5%9F%BA%E5%9B%A0\n"
        result = self.annotate(cds)
        assert gvann_of(result) == "T|STOP_GAINED|HIGH|tx基因|4|AAA>TAA|K>*"

    def test_lowercase_hex_is_equivalent(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%e5%9f%ba%e5%9b%a0\n"
        result = self.annotate(cds)
        assert gvann_of(result) == "T|STOP_GAINED|HIGH|tx基因|4|AAA>TAA|K>*"

    def test_escaped_and_plain_forms_group_into_one_transcript(self) -> None:
        # The second fragment's CDS_POS only makes sense if both rows
        # assembled into a single 15-base transcript.
        cds = (
            "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=tx%E5%9F%BA%E5%9B%A0\n"
            "chr1\tx\tCDS\t10\t15\t.\t+\t0\tParent=tx基因\n"
        )
        text = HEADER + COLUMNS + "chr1\t12\t.\tG\tA\t.\t.\t.\n"
        result = self.annotate(cds, text)
        assert gvann_of(result) == "A|SYNONYMOUS|LOW|tx基因|12|GGG>GGA|G>G"

    def test_escaped_comma_is_one_parent(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%2C1\n"
        result = self.annotate(cds)
        # One transcript; the decoded comma is re-encoded in GVANN.
        assert gvann_of(result) == "T|STOP_GAINED|HIGH|tx%2C1|4|AAA>TAA|K>*"

    def test_multi_value_parent_still_splits_on_plain_comma(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=txB,tx%41\n"
        result = self.annotate(cds)
        # Two transcripts, ordered by the decoded Parent (txA < txB).
        assert gvann_of(result) == (
            "T|STOP_GAINED|HIGH|txA|4|AAA>TAA|K>*,"
            "T|STOP_GAINED|HIGH|txB|4|AAA>TAA|K>*"
        )

    def test_transcripts_sorted_by_decoded_parent(self) -> None:
        cds = (
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=txB\n"
            "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%41\n"
        )
        result = self.annotate(cds)
        assert gvann_of(result) == (
            "T|STOP_GAINED|HIGH|txA|4|AAA>TAA|K>*,"
            "T|STOP_GAINED|HIGH|txB|4|AAA>TAA|K>*"
        )

    def test_structural_characters_reencoded_in_transcript(self) -> None:
        # Decoded "a=b;c|d e%" would break the GVANN/VCF structure.
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=a%3Db%3Bc%7Cd%20e%25\n"
        result = self.annotate(cds)
        assert gvann_of(result) == (
            "T|STOP_GAINED|HIGH|a%3Db%3Bc%7Cd%20e%25|4|AAA>TAA|K>*"
        )

    def test_non_ascii_transcript_emitted_directly(self) -> None:
        cds = "chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx基因\n"
        result = self.annotate(cds)
        assert gvann_of(result) == "T|STOP_GAINED|HIGH|tx基因|4|AAA>TAA|K>*"

    def test_truncated_escape_rejected(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%4\n")
        assert "<stream>:1:" in str(excinfo.value)

    def test_lone_percent_rejected(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%\n")
        assert "<stream>:1:" in str(excinfo.value)

    def test_non_hex_escape_rejected(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%GG\n")
        assert "<stream>:1:" in str(excinfo.value)

    def test_non_hex_sign_escape_rejected(self) -> None:
        # "%+1" is not a hexadecimal escape even though int() parses it.
        with pytest.raises(AnnotationFormatError):
            self.annotate("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%+1\n")

    def test_invalid_utf8_escape_bytes_rejected(self) -> None:
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%FF\n")
        assert "<stream>:1:" in str(excinfo.value)

    def test_incomplete_utf8_sequence_rejected(self) -> None:
        # E5 9F starts a three-byte sequence that never completes.
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate("chr1\tx\tCDS\t1\t15\t.\t+\t0\tParent=tx%E5%9F\n")
        assert "<stream>:1:" in str(excinfo.value)

    def test_error_reports_source_line(self) -> None:
        cds = (
            "chr1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
            "chr1\tx\tCDS\t10\t15\t.\t+\t0\tParent=t%FF\n"
        )
        with pytest.raises(AnnotationFormatError) as excinfo:
            self.annotate(cds)
        assert "<stream>:2:" in str(excinfo.value)


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
