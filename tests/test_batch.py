"""Tests for genome_variant.batch."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.batch import (
    BatchManifestEntry,
    BatchManifestError,
    batch_call_variants,
    read_batch_manifest,
)
from genome_variant.calling import VariantCallingError, call_variants
from genome_variant.sequence_io import (
    SequenceFormatError,
    SequenceRecord,
    read_sequences,
)
from genome_variant.summary import render_summaries, summarize_variants
from genome_variant.vcf import render_vcf

REF = "ACGTACGTACAA"
REF_FASTA = f">c1\n{REF}\n"


def fastq(identifier, sequence, quality=None):
    if quality is None:
        quality = "I" * len(sequence)
    return f"@{identifier}\n{sequence}\n+{identifier}\n{quality}\n"


def alt_read(identifier, base="T", position=4):
    sequence = list(REF)
    sequence[position] = base
    return fastq(identifier, "".join(sequence))


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def het_reads():
    return alt_read("a") + alt_read("b") + fastq("c", REF)


def hom_reads():
    return alt_read("a") + alt_read("b") + alt_read("c")


class TestReadBatchManifest:
    def test_parses_entries_in_order(self, tmp_path) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"a.fq"}\n\n{"sample":"s2","reads":"b.fq"}\n',
        )
        entries = read_batch_manifest(str(manifest))
        assert [entry.sample for entry in entries] == ["s1", "s2"]
        assert [entry.line_number for entry in entries] == [1, 3]
        assert entries[0].reads == str(tmp_path / "a.fq")

    def test_file_manifest_resolves_relative_reads_vs_its_directory(
        self, tmp_path, monkeypatch
    ) -> None:
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s","reads":"r.fq"}\n')
        monkeypatch.chdir("/")
        (entry,) = read_batch_manifest(str(manifest))
        assert entry.reads == str(tmp_path / "r.fq")

    def test_stream_manifest_resolves_relative_reads_vs_cwd(
        self, tmp_path, monkeypatch
    ) -> None:
        monkeypatch.chdir(tmp_path)
        (entry,) = read_batch_manifest(StringIO('{"sample":"s","reads":"r.fq"}\n'))
        assert entry.reads == str(tmp_path / "r.fq")

    def test_absolute_reads_kept_verbatim(self, tmp_path) -> None:
        (entry,) = read_batch_manifest(
            StringIO('{"sample":"s","reads":"/abs/r.fq"}\n')
        )
        assert entry.reads == "/abs/r.fq"

    def test_stream_without_name_uses_placeholder(self) -> None:
        stream = StringIO("{bad json}\n")
        with pytest.raises(BatchManifestError, match=r"<stream>:1:"):
            read_batch_manifest(stream)

    def test_bad_json_reports_source_and_line(self, tmp_path) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"a","reads":"x"}\n{not json}\n',
        )
        with pytest.raises(BatchManifestError) as excinfo:
            read_batch_manifest(str(manifest))
        assert f"{manifest}:2:" in str(excinfo.value)

    def test_line_that_is_not_object(self) -> None:
        with pytest.raises(BatchManifestError, match=r":1: .*JSON object"):
            read_batch_manifest(StringIO('["sample","reads"]\n'))

    @pytest.mark.parametrize(
        "text",
        [
            '{"sample":"a"}\n',                          # missing reads
            '{"reads":"x"}\n',                           # missing sample
            '{"sample":"a","reads":"x","extra":1}\n',    # extra field
            '{"sample":"","reads":"x"}\n',               # empty sample
            '{"sample":"  ","reads":"x"}\n',             # whitespace-only sample
            '{"sample":"a\tb","reads":"x"}\n',           # tab in sample
            '{"sample":3,"reads":"x"}\n',                # non-string sample
            '{"sample":"a","reads":""}\n',               # empty reads
            '{"sample":"a","reads":4}\n',                # non-string reads
            '{"sample":"a","reads":"-"}\n',              # standard input reads
        ],
    )
    def test_bad_fields_raise_batch_manifest_error(self, text) -> None:
        with pytest.raises(BatchManifestError) as excinfo:
            read_batch_manifest(StringIO(text))
        assert ":1:" in str(excinfo.value)

    def test_duplicate_sample_reports_second_line(self) -> None:
        text = '{"sample":"a","reads":"x"}\n{"sample":"a","reads":"y"}\n'
        with pytest.raises(BatchManifestError, match=r":2: duplicate"):
            read_batch_manifest(StringIO(text))

    def test_missing_manifest_file_raises_oserror(self, tmp_path) -> None:
        with pytest.raises(OSError):
            read_batch_manifest(str(tmp_path / "nope.jsonl"))

    def test_bad_source_type_raises_typeerror(self) -> None:
        with pytest.raises(TypeError):
            read_batch_manifest(42)


class TestBatchCallVariants:
    def test_matches_call_variants_plus_summarize(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        reads1 = write(tmp_path, "s1.fq", het_reads())
        reads2 = write(tmp_path, "s2.fq", hom_reads())
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"s1.fq"}\n'
            '{"sample":"s2","reads":"s2.fq"}\n',
        )
        summaries = batch_call_variants(str(manifest), str(reference))
        assert len(summaries) == 1
        summary = summaries[0]
        assert (summary.chrom, summary.pos, summary.ref, summary.alt) == (
            "c1",
            5,
            "A",
            "T",
        )
        assert summary.sample_count == 2
        assert summary.allele_count == 3
        assert summary.depth == 6
        assert [(c.sample, c.gt, c.dp, c.ad) for c in summary.samples] == [
            ("s1", "0/1", 3, (1, 2)),
            ("s2", "1/1", 3, (0, 3)),
        ]

    def test_equivalent_to_summarize_of_per_sample_calls(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        manifest_lines = []
        vcf_lines = []
        for sample, text in (("s1", het_reads()), ("s2", hom_reads())):
            reads_path = write(tmp_path, f"{sample}.fq", text)
            manifest_lines.append(f'{{"sample":"{sample}","reads":"{sample}.fq"}}')
            with open(reads_path) as stream:
                records = list(read_sequences(stream, format="fastq"))
            with open(reference) as stream:
                refs = list(read_sequences(stream, format="fasta"))
            vcf_path = write(
                tmp_path,
                f"{sample}.vcf",
                render_vcf(call_variants(refs, records, sample_name=sample)),
            )
            vcf_lines.append(f'{{"sample":"{sample}","vcf":"{sample}.vcf"}}')
        manifest = write(tmp_path, "m.jsonl", "\n".join(manifest_lines) + "\n")
        vcf_manifest = write(tmp_path, "v.jsonl", "\n".join(vcf_lines) + "\n")

        assert batch_call_variants(str(manifest), str(reference)) == (
            summarize_variants(str(vcf_manifest), str(reference))
        )

    def test_samples_keep_manifest_order(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        write(tmp_path, "a.fq", het_reads())
        write(tmp_path, "b.fq", het_reads())
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"b","reads":"b.fq"}\n{"sample":"a","reads":"a.fq"}\n',
        )
        (summary,) = batch_call_variants(str(manifest), str(reference))
        assert [call.sample for call in summary.samples] == ["b", "a"]

    def test_empty_manifest_gives_empty_result(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        manifest = write(tmp_path, "m.jsonl", "\n")
        assert batch_call_variants(str(manifest), str(reference)) == ()

    def test_no_variants_gives_empty_result(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        write(tmp_path, "s1.fq", fastq("a", REF))
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"s1.fq"}\n')
        assert batch_call_variants(str(manifest), str(reference)) == ()

    def test_accepts_entry_sequence_and_record_iterable(self, tmp_path) -> None:
        write(tmp_path, "s1.fq", het_reads())
        entries = (
            BatchManifestEntry(1, "s1", str(tmp_path / "s1.fq")),
        )
        refs = [SequenceRecord("c1", REF)]
        (summary,) = batch_call_variants(entries, refs)
        assert summary.samples[0].sample == "s1"

    def test_thresholds_match_per_sample_calling(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        write(tmp_path, "s1.fq", het_reads())
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"s1.fq"}\n')
        # AC=2 of DP=3; raising the ALT count threshold above 2 drops the call.
        assert batch_call_variants(
            str(manifest), str(reference), min_alt_count=3
        ) == ()
        (summary,) = batch_call_variants(
            str(manifest), str(reference), homozygous_fraction=0.6
        )
        assert summary.samples[0].gt == "1/1"

    def test_call_indels(self, tmp_path) -> None:
        indel_ref = "ACGATCGTACGGATCCGTAGCTAACCGGTTAC"
        inserted = indel_ref[:15] + "TT" + indel_ref[15:]
        reference = write(tmp_path, "r.fa", f">c1\n{indel_ref}\n")
        write(
            tmp_path,
            "s1.fq",
            fastq("a", inserted) + fastq("b", inserted) + fastq("w", indel_ref),
        )
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"s1.fq"}\n')
        assert batch_call_variants(str(manifest), str(reference)) == ()
        (summary,) = batch_call_variants(
            str(manifest), str(reference), call_indels=True
        )
        assert (summary.pos, summary.ref, summary.alt) == (15, "C", "CTT")
        assert summary.samples[0].ad == (1, 2)

    def test_invalid_thresholds_raise_before_consuming_inputs(
        self, tmp_path
    ) -> None:
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"gone.fq"}\n')
        reference = write(tmp_path, "r.fa", REF_FASTA)
        with pytest.raises(ValueError):
            batch_call_variants(
                str(manifest),
                str(reference),
                min_alt_fraction=0.9,
                homozygous_fraction=0.8,
            )
        with pytest.raises(ValueError):
            batch_call_variants(str(manifest), str(reference), min_base_quality=94)

    def test_empty_reference_raises_calling_error(self, tmp_path) -> None:
        write(tmp_path, "s1.fq", het_reads())
        entries = (BatchManifestEntry(1, "s1", str(tmp_path / "s1.fq")),)
        with pytest.raises(VariantCallingError):
            batch_call_variants(entries, [])

    def test_duplicate_reference_identifier_raises_calling_error(
        self, tmp_path
    ) -> None:
        reference = write(tmp_path, "r.fa", ">c1\nACGT\n>c1\nTTTT\n")
        write(tmp_path, "s1.fq", het_reads())
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"s1.fq"}\n')
        with pytest.raises(VariantCallingError, match="duplicate"):
            batch_call_variants(str(manifest), str(reference))

    def test_manifest_errors_raise_batch_manifest_error(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        manifest = write(tmp_path, "m.jsonl", '{"sample":"a","reads":"-"}\n')
        with pytest.raises(BatchManifestError):
            batch_call_variants(str(manifest), str(reference))

    def test_fastq_format_error_propagates(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        write(tmp_path, "s1.fq", "@a\nACGT\n+a\n")  # truncated record
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"s1.fq"}\n')
        with pytest.raises(SequenceFormatError):
            batch_call_variants(str(manifest), str(reference))

    def test_missing_reads_file_raises_oserror(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"gone.fq"}\n')
        with pytest.raises(OSError):
            batch_call_variants(str(manifest), str(reference))

    def test_byte_stable_across_runs(self, tmp_path) -> None:
        reference = write(tmp_path, "r.fa", REF_FASTA)
        write(tmp_path, "s1.fq", het_reads())
        write(tmp_path, "s2.fq", hom_reads())
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"s1.fq"}\n{"sample":"s2","reads":"s2.fq"}\n',
        )
        first = render_summaries(batch_call_variants(str(manifest), str(reference)))
        second = render_summaries(batch_call_variants(str(manifest), str(reference)))
        assert first == second
