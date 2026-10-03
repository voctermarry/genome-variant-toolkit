"""Tests for the ``genome_variant.batch`` public interface."""

from __future__ import annotations

from io import StringIO

import pytest

from genome_variant.batch import (
    BatchManifestEntry,
    BatchManifestError,
    batch_call_variants,
    read_batch_manifest,
)
from genome_variant.calling import VariantCallingError
from genome_variant.sequence_io import SequenceRecord

REFERENCE = ">c1\nACGTACGTACAA\n"
ALT_READ = "ACGTTCGTACAA"
REF_READ = "ACGTACGTACAA"


def fastq(identifier, sequence):
    return f"@{identifier}\n{sequence}\n+{identifier}\n{'I' * len(sequence)}\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.fixture
def reference_file(tmp_path):
    return write(tmp_path, "ref.fa", REFERENCE)


@pytest.fixture
def reads_dir(tmp_path):
    directory = tmp_path / "reads"
    directory.mkdir()
    (directory / "s1.fq").write_text(
        fastq("a", ALT_READ) + fastq("b", ALT_READ) + fastq("c", REF_READ)
    )
    (directory / "s2.fq").write_text(fastq("d", ALT_READ) + fastq("e", ALT_READ))
    return directory


class TestReadBatchManifest:
    def test_file_manifest_resolves_relative_reads_against_its_directory(
        self, tmp_path, reads_dir, monkeypatch
    ) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
            '{"sample":"s2","reads":"reads/s2.fq"}\n',
        )
        monkeypatch.chdir("/")
        entries = read_batch_manifest(str(manifest))
        assert [entry.sample for entry in entries] == ["s1", "s2"]
        assert entries[0].line_number == 1
        assert entries[1].line_number == 2
        assert entries[0].reads == str(reads_dir / "s1.fq")
        assert entries[1].reads == str(reads_dir / "s2.fq")

    def test_stream_manifest_resolves_relative_reads_against_cwd(
        self, tmp_path, reads_dir, monkeypatch
    ) -> None:
        monkeypatch.chdir(tmp_path)
        entries = read_batch_manifest(StringIO('{"sample":"s1","reads":"reads/s1.fq"}\n'))
        assert entries[0].reads == str(reads_dir / "s1.fq")

    def test_absolute_reads_kept(self, tmp_path, reads_dir) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"%s"}\n' % (reads_dir / "s1.fq"),
        )
        entries = read_batch_manifest(str(manifest))
        assert entries[0].reads == str(reads_dir / "s1.fq")

    def test_blank_lines_ignored(self, tmp_path) -> None:
        manifest = write(
            tmp_path, "m.jsonl", '\n{"sample":"s1","reads":"x.fq"}\n   \n'
        )
        entries = read_batch_manifest(str(manifest))
        assert len(entries) == 1
        assert entries[0].line_number == 2

    def test_empty_manifest(self, tmp_path) -> None:
        manifest = write(tmp_path, "m.jsonl", "\n")
        assert read_batch_manifest(str(manifest)) == ()

    @pytest.mark.parametrize(
        "line",
        [
            "{not json}",
            '["s1","x.fq"]',
            '{"sample":"s1"}',
            '{"reads":"x.fq"}',
            '{"sample":"s1","reads":"x.fq","extra":1}',
            '{"sample":"","reads":"x.fq"}',
            '{"sample":"   ","reads":"x.fq"}',
            '{"sample":"a\\tb","reads":"x.fq"}',
            '{"sample":3,"reads":"x.fq"}',
            '{"sample":"s1","reads":""}',
            '{"sample":"s1","reads":4}',
            '{"sample":"s1","reads":"-"}',
        ],
    )
    def test_invalid_lines_raise_with_source_and_line(self, tmp_path, line) -> None:
        manifest = write(tmp_path, "m.jsonl", line + "\n")
        with pytest.raises(BatchManifestError) as info:
            read_batch_manifest(str(manifest))
        message = str(info.value)
        assert str(manifest) in message
        assert ":1:" in message

    def test_duplicate_sample_reports_second_line(self, tmp_path) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"a.fq"}\n{"sample":"s1","reads":"b.fq"}\n',
        )
        with pytest.raises(BatchManifestError) as info:
            read_batch_manifest(str(manifest))
        assert ":2:" in str(info.value)
        assert "duplicate" in str(info.value)

    def test_missing_manifest_file_raises_oserror(self, tmp_path) -> None:
        with pytest.raises(OSError):
            read_batch_manifest(str(tmp_path / "gone.jsonl"))


class TestBatchCallVariants:
    def test_basic_merge(self, tmp_path, reference_file, reads_dir) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
            '{"sample":"s2","reads":"reads/s2.fq"}\n',
        )
        summaries = batch_call_variants(str(manifest), str(reference_file))
        assert len(summaries) == 1
        summary = summaries[0]
        assert (summary.chrom, summary.pos, summary.ref, summary.alt) == (
            "c1",
            5,
            "A",
            "T",
        )
        assert summary.sample_count == 2
        assert summary.allele_count == 3  # 0/1 (1) + 1/1 (2)
        assert summary.depth == 5
        assert [call.sample for call in summary.samples] == ["s1", "s2"]
        assert summary.samples[0].gt == "0/1"
        assert summary.samples[0].dp == 3
        assert summary.samples[0].ad == (1, 2)
        assert summary.samples[1].gt == "1/1"
        assert summary.samples[1].dp == 2
        assert summary.samples[1].ad == (0, 2)

    def test_empty_manifest_produces_empty_result(
        self, tmp_path, reference_file
    ) -> None:
        manifest = write(tmp_path, "m.jsonl", "")
        assert batch_call_variants(str(manifest), str(reference_file)) == ()

    def test_no_variants_produces_empty_result(
        self, tmp_path, reference_file, reads_dir
    ) -> None:
        (reads_dir / "w.fq").write_text(fastq("w", REF_READ))
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"w","reads":"reads/w.fq"}\n'
        )
        assert batch_call_variants(str(manifest), str(reference_file)) == ()

    def test_thresholds_match_single_sample_calling(
        self, tmp_path, reference_file, reads_dir
    ) -> None:
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"s1","reads":"reads/s1.fq"}\n'
        )
        # AC=2 of DP=3: dropped by a higher count threshold...
        assert (
            batch_call_variants(str(manifest), str(reference_file), min_alt_count=3)
            == ()
        )
        # ...and genotyped 1/1 by a lower homozygous threshold.
        summaries = batch_call_variants(
            str(manifest), str(reference_file), homozygous_fraction=0.5
        )
        assert summaries[0].samples[0].gt == "1/1"

    def test_accepts_entries_and_record_iterables(
        self, tmp_path, reads_dir
    ) -> None:
        entries = (
            BatchManifestEntry(1, "s1", str(reads_dir / "s1.fq")),
            BatchManifestEntry(2, "s2", str(reads_dir / "s2.fq")),
        )
        references = [SequenceRecord("c1", "ACGTACGTACAA")]
        summaries = batch_call_variants(entries, references)
        assert len(summaries) == 1
        assert [call.sample for call in summaries[0].samples] == ["s1", "s2"]

    def test_empty_reference_raises_calling_error(self, tmp_path) -> None:
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"s1","reads":"reads/s1.fq"}\n'
        )
        with pytest.raises(VariantCallingError):
            batch_call_variants(str(manifest), [])

    def test_duplicate_reference_identifier_raises_calling_error(
        self, tmp_path
    ) -> None:
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"s1","reads":"reads/s1.fq"}\n'
        )
        references = [
            SequenceRecord("c1", "ACGT"),
            SequenceRecord("c1", "TTTT"),
        ]
        with pytest.raises(VariantCallingError):
            batch_call_variants(str(manifest), references)

    def test_invalid_threshold_raises_value_error(
        self, tmp_path, reference_file, reads_dir
    ) -> None:
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"s1","reads":"reads/s1.fq"}\n'
        )
        with pytest.raises(ValueError):
            batch_call_variants(
                str(manifest), str(reference_file), min_alt_count=0
            )
        with pytest.raises(ValueError):
            batch_call_variants(
                str(manifest),
                str(reference_file),
                min_alt_fraction=0.9,
                homozygous_fraction=0.5,
            )

    def test_missing_reads_file_raises_oserror(
        self, tmp_path, reference_file
    ) -> None:
        manifest = write(
            tmp_path, "m.jsonl", '{"sample":"s1","reads":"reads/gone.fq"}\n'
        )
        with pytest.raises(OSError):
            batch_call_variants(str(manifest), str(reference_file))

    def test_manifest_error_type(self, tmp_path, reference_file) -> None:
        manifest = write(tmp_path, "m.jsonl", '{"sample":"s1","reads":"-"}\n')
        with pytest.raises(BatchManifestError):
            batch_call_variants(str(manifest), str(reference_file))

    def test_results_are_deterministic(
        self, tmp_path, reference_file, reads_dir
    ) -> None:
        manifest = write(
            tmp_path,
            "m.jsonl",
            '{"sample":"s1","reads":"reads/s1.fq"}\n'
            '{"sample":"s2","reads":"reads/s2.fq"}\n',
        )
        first = batch_call_variants(str(manifest), str(reference_file))
        second = batch_call_variants(str(manifest), str(reference_file))
        assert first == second
