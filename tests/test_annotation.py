"""Tests for GFF3-driven VCF consequence annotation."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.annotation import AnnotationFormatError, annotate_vcf
from genome_variant.sequence_io import SequenceRecord
from genome_variant.vcf import (
    ReferenceMismatchError,
    VcfFormatError,
    read_vcf,
    render_vcf,
)

# Positions:
#  1..3   ATG  M
#  4..6   AAA  K
#  7..9   TTT  F
#  10..12 GGG  G
#  13..15 CCC  P
#  16..18 TAA  *
#  19..21 TTA  (forward)
#  22..24 GGG
REFERENCE = "ATGAAATTTGGGCCCTAATTAGGG"
REFERENCE_RECORDS = [SequenceRecord("chr1", REFERENCE)]

HEADER = "##fileformat=VCFv4.2\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
GVANN_META = '##INFO=<ID=GVANN,Number=.,Type=String,Description="Genomic Variant Annotation: ALT|CONSEQUENCE|IMPACT|TRANSCRIPT|CDS_POS|CODON_CHANGE|AA_CHANGE">'


def annotate(
    vcf_text: str,
    gff_text: str,
    reference=REFERENCE_RECORDS,
) -> str:
    document = read_vcf(StringIO(vcf_text))
    annotated = annotate_vcf(document, reference, StringIO(gff_text))
    return render_vcf(annotated)


def record_lines(text: str) -> list[str]:
    return [line for line in text.splitlines() if not line.startswith("#")]


PLUS_GFF = "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx1\n"


class TestConsequences:
    def test_start_lost_first_codon(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", PLUS_GFF
        )
        line = record_lines(out)[0]
        assert line.endswith(
            "GVANN=T|START_LOST|HIGH|tx1|1|ATG>TTG|M1L"
        )

    def test_start_lost_within_first_codon(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t2\t.\tT\tC\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=C|START_LOST|HIGH|tx1|2|ATG>ACG|M1T"
        )

    def test_stop_gained(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=T|STOP_GAINED|HIGH|tx1|4|AAA>TAA|K2*"
        )

    def test_stop_lost(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t16\t.\tT\tC\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=C|STOP_LOST|HIGH|tx1|16|TAA>CAA|*6Q"
        )

    def test_synonymous(self) -> None:
        # AAA (K) -> AAG (K) at CDS position 6.
        out = annotate(
            HEADER + COLUMNS + "chr1\t6\t.\tA\tG\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=G|SYNONYMOUS|LOW|tx1|6|AAA>AAG|K2K"
        )

    def test_missense(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=C|MISSENSE|MODERATE|tx1|7|TTT>CTT|F3L"
        )

    def test_stop_to_stop_is_synonymous(self) -> None:
        # TAA -> TAG at the terminal stop codon.
        out = annotate(
            HEADER + COLUMNS + "chr1\t18\t.\tA\tG\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=G|SYNONYMOUS|LOW|tx1|18|TAA>TAG|*6*"
        )


class TestFallbackConsequences:
    def test_non_coding_snv_outside_cds(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t22\t.\tG\tT\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith("GVANN=T|NON_CODING|MODIFIER||||")

    def test_insertion_unsupported(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\tAT\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith("GVANN=AT|UNSUPPORTED|MODIFIER||||")

    def test_deletion_unsupported(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tAT\tA\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith("GVANN=A|UNSUPPORTED|MODIFIER||||")

    def test_symbolic_unsupported(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\t<DEL>\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith(
            "GVANN=<DEL>|UNSUPPORTED|MODIFIER||||"
        )

    def test_spanning_and_missing_unsupported(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\t*,.\t.\t.\t.\n", PLUS_GFF
        )
        line = record_lines(out)[0]
        assert line.endswith(
            "GVANN=*|UNSUPPORTED|MODIFIER||||,.|UNSUPPORTED|MODIFIER||||"
        )

    def test_ambiguous_alt_unsupported(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\tN\t.\t.\t.\n", PLUS_GFF
        )
        assert record_lines(out)[0].endswith("GVANN=N|UNSUPPORTED|MODIFIER||||")

    def test_ambiguous_base_in_reference_codon_unsupported(self) -> None:
        # Replace the third reference base (G) with N; the start codon
        # cannot be translated, so the SNV is unsupported.
        reference = [SequenceRecord("chr1", "ATN" + REFERENCE[3:])]
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n",
            PLUS_GFF,
            reference=reference,
        )
        assert record_lines(out)[0].endswith("GVANN=T|UNSUPPORTED|MODIFIER||||")

    def test_mixed_alts_keep_order(self) -> None:
        row = "chr1\t1\t.\tA\tT,<DEL>,G\t.\t.\t.\n"
        out = annotate(HEADER + COLUMNS + row, PLUS_GFF)
        line = record_lines(out)[0]
        assert line.endswith(
            "GVANN="
            "T|START_LOST|HIGH|tx1|1|ATG>TTG|M1L,"
            "<DEL>|UNSUPPORTED|MODIFIER||||,"
            "G|START_LOST|HIGH|tx1|1|ATG>GTG|M1V"
        )


class TestStrandAndFrame:
    def test_minus_strand_single_fragment(self) -> None:
        # Forward 19..21 = TTA; reverse complement (coding) = TAA stop.
        gff = "chr1\ttest\tCDS\t19\t21\t.\t-\t0\tID=c2;Parent=tx2\n"
        out = annotate(
            HEADER + COLUMNS + "chr1\t21\t.\tA\tG\t.\t.\t.\n", gff
        )
        assert record_lines(out)[0].endswith(
            "GVANN=G|STOP_LOST|HIGH|tx2|1|TAA>CAA|*1Q"
        )

    def test_minus_strand_cds_position_grows_leftwards(self) -> None:
        gff = "chr1\ttest\tCDS\t19\t21\t.\t-\t0\tID=c2;Parent=tx2\n"
        out = annotate(
            HEADER + COLUMNS + "chr1\t19\t.\tT\tC\t.\t.\t.\n", gff
        )
        # TTA -> CTA on the forward strand; coding TAA -> TAG.
        assert record_lines(out)[0].endswith(
            "GVANN=C|SYNONYMOUS|LOW|tx2|3|TAA>TAG|*1*"
        )

    def test_minus_strand_multi_fragment(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t19\t21\t.\t-\t0\tID=c1;Parent=txB\n"
            "chr1\ttest\tCDS\t16\t18\t.\t-\t0\tID=c2;Parent=txB\n"
        )
        out = annotate(
            HEADER + COLUMNS + "chr1\t17\t.\tA\tC\t.\t.\t.\n", gff
        )
        # Coding = TAA (19..21) + TTA (16..18); genomic 17 is coding
        # position 5, the middle base of the second codon.  A genomic
        # A->C change is a coding T->G change (reverse complement),
        # TTA (L) -> TGA (*).
        assert record_lines(out)[0].endswith(
            "GVANN=C|STOP_GAINED|HIGH|txB|5|TTA>TGA|L2*"
        )

    def test_plus_strand_multi_fragment_phase_zero(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=txA\n"
            "chr1\ttest\tCDS\t10\t12\t.\t+\t0\tID=c2;Parent=txA\n"
        )
        out = annotate(
            HEADER + COLUMNS + "chr1\t10\t.\tG\tA\t.\t.\t.\n", gff
        )
        assert record_lines(out)[0].endswith(
            "GVANN=A|MISSENSE|MODERATE|txA|7|GGG>AGG|G3R"
        )

    def test_internal_phase_must_match_frame(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=txA\n"
            "chr1\ttest\tCDS\t10\t12\t.\t+\t1\tID=c2;Parent=txA\n"
        )
        with pytest.raises(AnnotationFormatError, match="inconsistent with the reading frame"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_internal_phase_two_is_accepted(self) -> None:
        # First fragment 1..4 (4 coding bases): consumed 4, so the next
        # fragment's phase must be (-4) % 3 == 2; second fragment 8..12
        # (5 bases) gives 9 coding bases total.
        gff = (
            "chr1\ttest\tCDS\t1\t4\t.\t+\t0\tID=c1;Parent=txA\n"
            "chr1\ttest\tCDS\t8\t12\t.\t+\t2\tID=c2;Parent=txA\n"
        )
        out = annotate(
            HEADER + COLUMNS + "chr1\t8\t.\tT\tC\t.\t.\t.\n", gff
        )
        # Coding sequence is ATGA + TTTGG: codons ATG, ATT, TGG.
        # Genomic position 8 is CDS position 5, the middle base of the
        # second codon ATT (I); ATT -> ACT (T).
        assert record_lines(out)[0].endswith(
            "GVANN=C|MISSENSE|MODERATE|txA|5|ATT>ACT|I2T"
        )

    def test_phase_two_leader_bases(self) -> None:
        # Standalone CDS 2..11 phase 2: one leader base at the 5' end,
        # coding starts at position 3.
        gff = "chr1\ttest\tCDS\t2\t11\t.\t+\t2\tID=c1;Parent=tx1\n"
        out = annotate(
            HEADER + COLUMNS + "chr1\t3\t.\tG\tC\t.\t.\t.\n", gff
        )
        assert record_lines(out)[0].endswith(
            "GVANN=C|MISSENSE|MODERATE|tx1|1|GAA>CAA|E1Q"
        )

    def test_leader_base_is_non_coding(self) -> None:
        gff = "chr1\ttest\tCDS\t2\t11\t.\t+\t2\tID=c1;Parent=tx1\n"
        out = annotate(
            HEADER + COLUMNS + "chr1\t2\t.\tT\tC\t.\t.\t.\n", gff
        )
        assert record_lines(out)[0].endswith("GVANN=C|NON_CODING|MODIFIER||||")

    def test_phase_one_leader_bases(self) -> None:
        # CDS 2..12 phase 1: two leader bases, coding starts at position 4.
        gff = "chr1\ttest\tCDS\t2\t12\t.\t+\t1\tID=c1;Parent=tx1\n"
        out = annotate(
            HEADER + COLUMNS + "chr1\t4\t.\tA\tT\t.\t.\t.\n", gff
        )
        assert record_lines(out)[0].endswith(
            "GVANN=T|STOP_GAINED|HIGH|tx1|1|AAA>TAA|K1*"
        )


class TestTranscriptOrder:
    def test_multi_transcript_hit_sorted_by_parent(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx_b\n"
            "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c2;Parent=tx_a\n"
        )
        out = annotate(
            HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n", gff
        )
        line = record_lines(out)[0]
        first = line.index("tx_a")
        second = line.index("tx_b")
        assert first < second

    def test_only_transcripts_on_same_sequence_considered(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr9\ttest\tCDS\t1\t18\t.\t+\t0\tID=c2;Parent=tx9\n"
        )
        records = [
            SequenceRecord("chr1", REFERENCE),
            SequenceRecord("chr9", REFERENCE),
        ]
        out = annotate(
            HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n",
            gff,
            reference=records,
        )
        line = record_lines(out)[0]
        assert "tx1" in line
        assert "tx9" not in line


class TestFormatErrors:
    def test_cds_with_eight_columns(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t18\t.\t+\t0\n"
        with pytest.raises(AnnotationFormatError) as info:
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)
        assert ":1:" in str(info.value)

    def test_bad_strand(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t18\t.\t.\t0\tID=c1;Parent=tx1\n"
        with pytest.raises(AnnotationFormatError, match="strand"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    @pytest.mark.parametrize("coords", [("x", "18"), ("1", "z"), ("0", "18"), ("5", "2")])
    def test_bad_coordinates(self, coords) -> None:
        gff = f"chr1\ttest\tCDS\t{coords[0]}\t{coords[1]}\t.\t+\t0\tID=c1;Parent=tx1\n"
        with pytest.raises(AnnotationFormatError, match="coordinates"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    @pytest.mark.parametrize("phase", ["3", ".", "x", ""])
    def test_bad_phase(self, phase) -> None:
        gff = f"chr1\ttest\tCDS\t1\t18\t.\t+\t{phase}\tID=c1;Parent=tx1\n"
        with pytest.raises(AnnotationFormatError, match="phase"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_missing_parent(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1\n"
        with pytest.raises(AnnotationFormatError, match="Parent"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_multiple_parents(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx1,tx2\n"
        with pytest.raises(AnnotationFormatError, match="more than one parent"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_percent_decoded_parent_name(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=some%20tx\n"
        out = annotate(
            HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n", gff
        )
        assert "|some tx|" in record_lines(out)[0]

    def test_fragments_spanning_sequences(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr9\ttest\tCDS\t1\t6\t.\t+\t0\tID=c2;Parent=tx1\n"
        )
        records = [SequenceRecord("chr1", REFERENCE), SequenceRecord("chr9", REFERENCE)]
        with pytest.raises(AnnotationFormatError, match="span multiple sequences"):
            annotate(
                HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n",
                gff,
                reference=records,
            )

    def test_strand_disagreement(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr1\ttest\tCDS\t10\t12\t.\t-\t0\tID=c2;Parent=tx1\n"
        )
        with pytest.raises(AnnotationFormatError, match="strand"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_overlapping_fragments(self) -> None:
        gff = (
            "chr1\ttest\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=tx1\n"
            "chr1\ttest\tCDS\t4\t9\t.\t+\t0\tID=c2;Parent=tx1\n"
        )
        with pytest.raises(AnnotationFormatError, match="overlap"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_frame_length_not_multiple_of_three(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t7\t.\t+\t0\tID=c1;Parent=tx1\n"
        with pytest.raises(AnnotationFormatError, match="in-frame"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_feature_line_number_is_one_based(self) -> None:
        gff = (
            "##gff-version 3\n"
            "# a comment\n"
            "chr1\ttest\tCDS\t1\t18\t.\t+\t9\tID=c1;Parent=tx1\n"
        )
        with pytest.raises(AnnotationFormatError) as info:
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)
        assert "<stream>:3:" in str(info.value)

    def test_cds_past_reference_end(self) -> None:
        gff = "chr1\ttest\tCDS\t1\t99\t.\t+\t0\tID=c1;Parent=tx1\n"
        with pytest.raises(AnnotationFormatError, match="runs past the end"):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_non_cds_features_are_ignored(self) -> None:
        gff = (
            "chr1\ttest\tgene\t1\t18\t.\t+\t.\tID=g1\n"
            "bogus line without enough columns\n"
            "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx1\n"
        )
        out = annotate(
            HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n", gff
        )
        assert "MISSENSE" in record_lines(out)[0]


class TestGvannConflicts:
    def test_gvann_info_declaration_in_header_rejected(self) -> None:
        text = (
            HEADER
            + '##INFO=<ID=GVANN,Number=.,Type=String,Description="x">\n'
            + COLUMNS
            + "chr1\t1\t.\tA\tT\t.\t.\t.\n"
        )
        with pytest.raises(AnnotationFormatError, match="already declares GVANN"):
            annotate(text, PLUS_GFF)

    def test_bare_gvann_meta_rejected(self) -> None:
        text = HEADER + "##GVANN=something\n" + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n"
        with pytest.raises(AnnotationFormatError):
            annotate(text, PLUS_GFF)

    def test_gvann_in_info_rejected(self) -> None:
        text = HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\tGVANN=x\n"
        with pytest.raises(AnnotationFormatError, match="already carries a GVANN"):
            annotate(text, PLUS_GFF)


class TestReferenceAndPreservation:
    def test_missing_chrom_raises_reference_mismatch(self) -> None:
        text = HEADER + COLUMNS + "chr9\t1\t.\tA\tT\t.\t.\t.\n"
        with pytest.raises(ReferenceMismatchError):
            annotate(text, PLUS_GFF)

    def test_ref_mismatch_raises_reference_mismatch(self) -> None:
        text = HEADER + COLUMNS + "chr1\t1\t.\tC\tT\t.\t.\t.\n"
        with pytest.raises(ReferenceMismatchError):
            annotate(text, PLUS_GFF)

    def test_record_out_of_bounds_raises_vcf_format(self) -> None:
        text = HEADER + COLUMNS + "chr1\t99\t.\tA\tT\t.\t.\t.\n"
        with pytest.raises(VcfFormatError):
            annotate(text, PLUS_GFF)

    def test_cds_on_missing_reference_sequence(self) -> None:
        gff = "chr9\ttest\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=tx9\n"
        with pytest.raises(ReferenceMismatchError):
            annotate(HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", gff)

    def test_mapping_reference_accepted(self) -> None:
        document = read_vcf(
            StringIO(HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n")
        )
        annotated = annotate_vcf(document, {"chr1": REFERENCE}, StringIO(PLUS_GFF))
        assert "MISSENSE" in render_vcf(annotated)

    def test_existing_info_fields_kept_and_gvann_appended(self) -> None:
        text = HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\tDP=9;AF=0.5\n"
        out = annotate(text, PLUS_GFF)
        line = record_lines(out)[0]
        assert "\tDP=9;AF=0.5;GVANN=" in line

    def test_missing_info_dot_starts_gvann(self) -> None:
        text = HEADER + COLUMNS + "chr1\t7\t.\tT\tC\t.\t.\t.\n"
        out = annotate(text, PLUS_GFF)
        fields = record_lines(out)[0].split("\t")
        assert fields[7].startswith("GVANN=")

    def test_sample_columns_pass_through(self) -> None:
        columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        row = "chr1\t7\t.\tT\tC\t.\tPASS\tDP=3\tGT:DP\t0/1:9\t1/1:2\n"
        out = annotate(HEADER + columns + row, PLUS_GFF)
        line = record_lines(out)[0]
        assert line.endswith("\tGT:DP\t0/1:9\t1/1:2")
        assert "\tPASS\tDP=3;GVANN=" in line

    def test_other_columns_and_records_unchanged(self) -> None:
        rows = (
            "chr1\t7\t.\tT\tC\t80\tPASS\t.\n"
            "chr1\t22\t.\tG\tT\t.\t.\t.\n"
        )
        out = annotate(HEADER + COLUMNS + rows, PLUS_GFF)
        first, second = record_lines(out)
        assert first.startswith("chr1\t7\t.\tT\tC\t80\tPASS\t")
        assert second.startswith("chr1\t22\t.\tG\tT\t.\t.\t")


class TestHeaderAndRendering:
    def test_unique_gvann_meta_added_last(self) -> None:
        out = annotate(
            HEADER + COLUMNS + "chr1\t1\t.\tA\tT\t.\t.\t.\n", PLUS_GFF
        )
        meta = [line for line in out.splitlines() if line.startswith("##")]
        assert meta == ["##fileformat=VCFv4.2", GVANN_META]

    def test_result_is_deterministic(self) -> None:
        text = HEADER + COLUMNS + "chr1\t1\t.\tA\tG,T\t.\t.\tDP=1\n"
        gff = (
            "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c1;Parent=tx_b\n"
            "chr1\ttest\tCDS\t1\t18\t.\t+\t0\tID=c2;Parent=tx_a\n"
        )
        first = annotate(text, gff)
        second = annotate(text, gff)
        assert first == second
        assert first.endswith("\n")
        assert not first.endswith("\n\n")
