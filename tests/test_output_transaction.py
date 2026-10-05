"""Uniform output-transaction regression tests for every file-output entry.

All subcommands that accept ``--output`` share one observable output
transaction, regardless of whether their payload is text sequences, VCF,
JSON Lines or a single JSON result:

* argument validation, all input reading and the full computation finish
  before the destination directory is touched;
* a regular destination is prepared through a temporary file in that
  directory, replaced exactly once only after the whole UTF-8 payload has
  been written, flushed and closed;
* successful bytes are stable across repeated runs and identical to the
  standard-output bytes (ordering, compact JSON, non-ASCII, newlines,
  empty-result and VCF-header behavior);
* an existing destination is replaced only on success; a failure in a
  later input stage leaves it byte-for-byte untouched and leaves no
  temporary file behind;
* standard output never creates a temporary file;
* an unwritable destination directory and a failing final replacement
  both exit 1 with one stderr line and empty stdout;
* secondary errors while cleaning up a failed transaction never override
  the original return code and error message.

These tests pin the shared contract without relaxing any entry's own
validation or error classification.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass
from io import StringIO

import pytest

from genome_variant import cli
from genome_variant.cli import main

META = "##fileformat=VCFv4.2\n"
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
SAMPLE_COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"


def run(argv, stdin: str = ""):
    old_in, old_out, old_err = sys.stdin, sys.stdout, sys.stderr
    sys.stdin = StringIO(stdin)
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    try:
        try:
            code = main(argv)
        except SystemExit as exc:
            code = exc.code
    finally:
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdin, sys.stdout, sys.stderr = old_in, old_out, old_err
    return code, out, err


def write(tmp, name: str, text: str) -> str:
    path = tmp / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)


def fastq(identifier: str, sequence: str, quality: str | None = None) -> str:
    if quality is None:
        quality = "I" * len(sequence)
    return f"@{identifier}\n{sequence}\n+{identifier}\n{quality}\n"


def sample_vcf(sample: str, rows: list[str]) -> str:
    return META + SAMPLE_COLS.format(s=sample) + "".join(r + "\n" for r in rows)


def sample_row(chrom="chr1", pos=2, ref="C", alt="G", gt="0/1", dp="10", ad="6,4"):
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP:AD\t"
        f"{gt}:{dp}:{ad}"
    )


@dataclass
class EntrySpec:
    command: str
    ext: str
    make_argv: object
    inplace_target: object
    make_late_failure: object | None = None
    make_empty_argv: object | None = None
    make_stdin: object | None = None


# -- Per-entry fixture builders --------------------------------------------

def _normalize(tmp):
    path = write(tmp, "in.fa", ">a géne\nacgt\n>b\nACGT\n")
    return ["normalize-sequences", path]


def _normalize_inplace(tmp):
    return tmp / "in.fa"


def _normalize_late(tmp):
    path = write(tmp, "in.fa", ">a\nACGT\n>b\nACGZ\n")
    return ["normalize-sequences", path], 2


def _normalize_stdin(tmp):
    return ["normalize-sequences", "-"], ">a\nacgt\n"


def _filter(tmp):
    path = write(tmp, "in.fq", fastq("r", "A" * 30))
    return ["filter-reads", path]


def _filter_late(tmp):
    text = fastq("r1", "A" * 30) + "@r2\nACGT\n+r2\nII\n"
    path = write(tmp, "in.fq", text)
    return ["filter-reads", path], 2


def _filter_empty(tmp):
    path = write(tmp, "in.fq", "@r\nACGT\n+r\n" + "!!!!" + "\n")
    return ["filter-reads", path]


def _filter_stdin(tmp):
    return ["filter-reads", "-"], fastq("r", "A" * 30)


def _dedup(tmp):
    path = write(tmp, "in.fa", ">a\nACGT\n>b\nACGT\n")
    return ["deduplicate-reads", path]


def _dedup_late(tmp):
    path = write(tmp, "in.fa", ">a\nACGT\n>b\nACGZ\n")
    return ["deduplicate-reads", path], 2


def _dedup_stdin(tmp):
    return ["deduplicate-reads", "-"], ">a\nACGT\n>b\nACGT\n"


def _kmer(tmp):
    path = write(tmp, "in.fa", ">x\nACGT\n")
    return ["kmer-index", path, "--k", "2"]


def _kmer_late(tmp):
    path = write(tmp, "in.fa", ">x\nACGZ\n")
    return ["kmer-index", path, "--k", "2"], 2


def _kmer_empty(tmp):
    path = write(tmp, "in.fa", ">r1\nNNN\n")
    return ["kmer-index", path, "--k", "3"]


def _kmer_stdin(tmp):
    return ["kmer-index", "-", "--k", "2"], ">x\nACGT\n"


def _align(tmp):
    ref = write(tmp, "ref.fa", ">r\nACGT\n")
    query = write(tmp, "q.fa", ">q\nACGT\n")
    return ["align-pair", ref, query]


def _align_inplace(tmp):
    return tmp / "ref.fa"


def _align_late(tmp):
    ref = write(tmp, "ref.fa", ">r\nAA\n")
    query = write(tmp, "q.fa", ">q\nAZ\n")
    return ["align-pair", ref, query], 2


def _normalize_vcf(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGT\n")
    vcf_path = write(
        tmp, "in.vcf", META + COLS + "chr1\t2\t.\tC\tT\t.\t.\t.\n"
    )
    return ["normalize-vcf", vcf_path, "--reference", ref]


def _vcf_inplace(tmp):
    return tmp / "in.vcf"


def _normalize_vcf_late(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGT\n")
    text = (
        META
        + COLS
        + "chr1\t2\t.\tC\tT\t.\t.\t.\n"
        + "chr1\t2\t.\tX\tY\t.\t.\t.\n"
    )
    vcf_path = write(tmp, "in.vcf", text)
    return ["normalize-vcf", vcf_path, "--reference", ref], 2


def _normalize_vcf_empty(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGT\n")
    vcf_path = write(tmp, "in.vcf", META + COLS)
    return ["normalize-vcf", vcf_path, "--reference", ref]


def _annotate(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nATGAAATTTGGGCCC\n")
    feats = write(
        tmp, "f.gff3", "chr1\ttest\tCDS\t1\t15\t.\t+\t0\tID=c1;Parent=t1\n"
    )
    vcf_path = write(
        tmp, "in.vcf", META + COLS + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
    )
    return ["annotate-vcf", vcf_path, "--reference", ref, "--features", feats]


def _annotate_late(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nATGAAATTTGGGCCC\n")
    feats = write(
        tmp, "f.gff3", "chr1\ttest\tCDS\t1\t15\t.\t+\t0\tID=c1;Parent=t1\n"
    )
    text = (
        META
        + COLS
        + "chr1\t4\t.\tA\tT\t.\t.\t.\n"
        + "chr1\t2\t.\tX\tC\t.\t.\t.\n"
    )
    vcf_path = write(tmp, "in.vcf", text)
    return (
        ["annotate-vcf", vcf_path, "--reference", ref, "--features", feats],
        2,
    )


def _annotate_empty(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nATGAAATTTGGGCCC\n")
    feats = write(
        tmp, "f.gff3", "chr1\ttest\tCDS\t1\t15\t.\t+\t0\tID=c1;Parent=t1\n"
    )
    vcf_path = write(tmp, "in.vcf", META + COLS)
    return ["annotate-vcf", vcf_path, "--reference", ref, "--features", feats]


def _map(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGT\n")
    reads = write(tmp, "reads.fa", ">a\nACGTACGT\n")
    return ["map-reads", ref, reads]


def _ref_inplace(tmp):
    return tmp / "ref.fa"


def _map_late(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGT\n")
    reads = write(tmp, "reads.fa", ">a\nACGTACGT\n>b\nACGZ\n")
    return ["map-reads", ref, reads], 2


def _call(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGT\n")
    reads = write(tmp, "reads.fq", fastq("a", "ACGTTCGT") + fastq("b", "ACGTTCGT"))
    return ["call-variants", ref, reads]


def _call_late(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGT\n")
    reads = write(tmp, "reads.fq", fastq("a", "ACGTACGT") + "@b\nAZ\n+b\nII\n")
    return ["call-variants", ref, reads], 2


def _call_empty(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGT\n")
    reads = write(tmp, "reads.fq", fastq("a", "ACGTACGT"))
    return ["call-variants", ref, reads]


def _coverage(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGT\n")
    reads = write(tmp, "reads.fa", ">a\nACGT\n")
    return ["coverage-report", ref, reads]


def _coverage_late(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGT\n")
    reads = write(tmp, "reads.fa", ">a\nACGT\n>b\nACGZ\n")
    return ["coverage-report", ref, reads], 2


def _batch(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGTACAA\n")
    (tmp / "reads").mkdir(exist_ok=True)
    write(
        tmp,
        "reads/s1.fq",
        fastq("a", "ACGTTCGTACAA") + fastq("b", "ACGTTCGTACAA"),
    )
    manifest = write(
        tmp, "m.jsonl", '{"sample":"s1","reads":"reads/s1.fq"}\n'
    )
    return ["batch-call-variants", manifest, "--reference", ref]


def _manifest_inplace(tmp):
    return tmp / "m.jsonl"


def _batch_late(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGTACAA\n")
    (tmp / "reads").mkdir(exist_ok=True)
    write(tmp, "reads/s1.fq", fastq("a", "ACGTTCGTACAA"))
    write(tmp, "reads/broken.fq", "@broken\nACGT\n+broken\nII\n")
    manifest = write(
        tmp,
        "m.jsonl",
        '{"sample":"s1","reads":"reads/s1.fq"}\n'
        '{"sample":"s2","reads":"reads/broken.fq"}\n',
    )
    return ["batch-call-variants", manifest, "--reference", ref], 2


def _batch_late_oserror(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGTACAA\n")
    (tmp / "reads").mkdir(exist_ok=True)
    write(tmp, "reads/s1.fq", fastq("a", "ACGTTCGTACAA"))
    manifest = write(
        tmp,
        "m.jsonl",
        '{"sample":"s1","reads":"reads/s1.fq"}\n'
        '{"sample":"s2","reads":"reads/missing.fq"}\n',
    )
    return ["batch-call-variants", manifest, "--reference", ref], 1


def _batch_empty(tmp):
    ref = write(tmp, "ref.fa", ">c1\nACGTACGTACAA\n")
    manifest = write(tmp, "blank.jsonl", "\n  \n\n")
    return ["batch-call-variants", manifest, "--reference", ref]


def _summarize(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGT\n")
    (tmp / "vcfs").mkdir(exist_ok=True)
    write(tmp, "vcfs/s1.vcf", sample_vcf("s1", [sample_row()]))
    manifest = write(
        tmp, "m.jsonl", '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
    )
    return ["summarize-variants", manifest, "--reference", ref]


def _summarize_late(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGT\n")
    (tmp / "vcfs").mkdir(exist_ok=True)
    write(tmp, "vcfs/s1.vcf", sample_vcf("s1", [sample_row()]))
    write(tmp, "vcfs/s2.vcf", "not-a-vcf\n")
    manifest = write(
        tmp,
        "m.jsonl",
        '{"sample":"s1","vcf":"vcfs/s1.vcf"}\n'
        '{"sample":"s2","vcf":"vcfs/s2.vcf"}\n',
    )
    return ["summarize-variants", manifest, "--reference", ref], 2


def _summarize_empty(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGT\n")
    (tmp / "vcfs").mkdir(exist_ok=True)
    write(tmp, "vcfs/empty.vcf", META + SAMPLE_COLS.format(s="s1"))
    manifest = write(
        tmp, "m.jsonl", '{"sample":"s1","vcf":"vcfs/empty.vcf"}\n'
    )
    return ["summarize-variants", manifest, "--reference", ref]


def _compare(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGGGGT\n>chr2\nTTTTAAAA\n")
    baseline = write(tmp, "b.vcf", sample_vcf("s1", [sample_row()]))
    candidate = write(
        tmp,
        "c.vcf",
        sample_vcf(
            "s2",
            [sample_row(gt="1/1"),
             sample_row(chrom="chr2", pos=5, ref="A", alt="G", ad="5,5")],
        ),
    )
    return ["compare-variants", baseline, candidate, "--reference", ref]


def _baseline_inplace(tmp):
    return tmp / "b.vcf"


def _compare_late(tmp):
    ref = write(tmp, "ref.fa", ">chr1\nACGGGGGT\n")
    baseline = write(tmp, "b.vcf", sample_vcf("s1", [sample_row()]))
    candidate = write(
        tmp,
        "c.vcf",
        META
        + SAMPLE_COLS.format(s="s2")
        + "chr1\t0\t.\tA\tG\t.\t.\t.\tGT\t0/1\n",
    )
    return ["compare-variants", baseline, candidate, "--reference", ref], 2


ENTRIES: dict[str, EntrySpec] = {
    "normalize-sequences": EntrySpec(
        "normalize-sequences", ".fa", _normalize, _normalize_inplace,
        _normalize_late, None, make_stdin=_normalize_stdin
    ),
    "filter-reads": EntrySpec(
        "filter-reads", ".fq", _filter, lambda t: t / "in.fq",
        _filter_late, _filter_empty, make_stdin=_filter_stdin
    ),
    "deduplicate-reads": EntrySpec(
        "deduplicate-reads", ".fa", _dedup, lambda t: t / "in.fa",
        _dedup_late, None, make_stdin=_dedup_stdin
    ),
    "kmer-index": EntrySpec(
        "kmer-index", ".jsonl", _kmer, lambda t: t / "in.fa",
        _kmer_late, _kmer_empty, make_stdin=_kmer_stdin
    ),
    "align-pair": EntrySpec(
        "align-pair", ".json", _align, _align_inplace, _align_late
    ),
    "normalize-vcf": EntrySpec(
        "normalize-vcf", ".vcf", _normalize_vcf, _vcf_inplace,
        _normalize_vcf_late, _normalize_vcf_empty
    ),
    "annotate-vcf": EntrySpec(
        "annotate-vcf", ".vcf", _annotate, _vcf_inplace,
        _annotate_late, _annotate_empty
    ),
    "map-reads": EntrySpec(
        "map-reads", ".jsonl", _map, _ref_inplace, _map_late
    ),
    "call-variants": EntrySpec(
        "call-variants", ".vcf", _call, _ref_inplace, _call_late, _call_empty
    ),
    "batch-call-variants": EntrySpec(
        "batch-call-variants", ".jsonl", _batch, _manifest_inplace,
        _batch_late, _batch_empty
    ),
    "coverage-report": EntrySpec(
        "coverage-report", ".jsonl", _coverage, _ref_inplace, _coverage_late
    ),
    "summarize-variants": EntrySpec(
        "summarize-variants", ".jsonl", _summarize, _manifest_inplace,
        _summarize_late, _summarize_empty
    ),
    "compare-variants": EntrySpec(
        "compare-variants", ".json", _compare, _baseline_inplace,
        _compare_late
    ),
}

# batch-call-variants also has a late *file-access* failure (exit 1).
LATE_OSERROR = {"batch-call-variants": _batch_late_oserror}

ENTRY_IDS = list(ENTRIES)


# -- Shared successful-transaction contract ---------------------------------

@pytest.mark.parametrize("name", ENTRY_IDS)
def test_file_output_is_byte_deterministic_and_cleans_temp(
    name, tmp_path, monkeypatch
):
    spec = ENTRIES[name]
    first = tmp_path / f"a{spec.ext}"
    second = tmp_path / f"b{spec.ext}"

    real_replace = os.replace
    counts = {"replace": 0}

    def counting_replace(src, dst):
        counts["replace"] += 1
        return real_replace(src, dst)

    monkeypatch.setattr(cli.os, "replace", counting_replace)

    code, out, err = run(spec.make_argv(tmp_path) + ["--output", str(first)])
    assert code == 0, err
    assert out == "" and err == ""
    assert counts["replace"] == 1

    code, _, err = run(spec.make_argv(tmp_path) + ["--output", str(second)])
    assert code == 0, err
    assert first.read_bytes() == second.read_bytes()

    # Repeated runs on the same target stay byte-for-byte identical and
    # the destination is replaced once per run.
    before = first.read_bytes()
    code, _, err = run(spec.make_argv(tmp_path) + ["--output", str(first)])
    assert code == 0 and err == ""
    assert first.read_bytes() == before
    assert counts["replace"] == 3

    # No temporary file survives the transaction.
    leftovers = [
        p.name for p in tmp_path.rglob("*")
        if p.is_file() and (p.name.endswith(".tmp") or p.name.startswith("."))
        and p.name not in {".keep"}
    ]
    assert leftovers == []


@pytest.mark.parametrize("name", ENTRY_IDS)
def test_file_bytes_equal_standard_output_bytes(name, tmp_path, monkeypatch):
    spec = ENTRIES[name]
    target = tmp_path / f"out{spec.ext}"

    real_replace = os.replace
    monkeypatch.setattr(cli.os, "replace", lambda s, d: real_replace(s, d))

    code_file, out_file, err_file = run(
        spec.make_argv(tmp_path) + ["--output", str(target)]
    )
    assert code_file == 0 and out_file == "" and err_file == ""

    # Standard output never goes near the replacement machinery.
    seen = {"replace": 0}

    def no_replace(src, dst):  # pragma: no cover - must never be called
        seen["replace"] += 1
        raise AssertionError("os.replace used on the standard-output path")

    monkeypatch.setattr(cli.os, "replace", no_replace)
    before = {p.name for p in tmp_path.rglob("*")}
    code_std, out_std, err_std = run(spec.make_argv(tmp_path))
    after = {p.name for p in tmp_path.rglob("*")}
    assert code_std == 0 and err_std == ""
    assert seen["replace"] == 0
    assert before == after

    assert target.read_bytes() == out_std.encode("utf-8")


# -- Destination sharing an input path --------------------------------------

@pytest.mark.parametrize("name", ENTRY_IDS)
def test_output_same_path_as_input_succeeds_in_place(name, tmp_path):
    spec = ENTRIES[name]
    argv = spec.make_argv(tmp_path)
    target = spec.inplace_target(tmp_path)
    assert target.exists()

    # The expected bytes are produced independently on the stdout path.
    code_ref, out_ref, err_ref = run(list(argv))
    assert code_ref == 0 and err_ref == ""

    code, out, err = run(list(argv) + ["--output", str(target)])
    assert code == 0, err
    assert out == "" and err == ""
    assert target.read_bytes() == out_ref.encode("utf-8")
    assert not [p for p in tmp_path.rglob("*.tmp")]


# -- Unwritable destination directory ---------------------------------------

@pytest.mark.parametrize("name", ENTRY_IDS)
def test_unwritable_target_directory_is_exit_1_without_target(name, tmp_path):
    spec = ENTRIES[name]
    missing = tmp_path / "no-such-dir" / f"out{spec.ext}"
    code, out, err = run(spec.make_argv(tmp_path) + ["--output", str(missing)])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1 and err.endswith("\n")
    assert not missing.exists()
    assert not (tmp_path / "no-such-dir").exists()
    assert not [p for p in tmp_path.rglob("*.tmp")]


# -- Final replacement failure ----------------------------------------------

@pytest.mark.parametrize("name", ENTRY_IDS)
def test_replacement_failure_is_exit_1_and_cleans_temp(name, tmp_path):
    spec = ENTRIES[name]
    destination = tmp_path / "destination"
    destination.mkdir()  # replacing a directory with a file fails

    code, out, err = run(
        spec.make_argv(tmp_path) + ["--output", str(destination)]
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1 and err.endswith("\n")
    # The directory is untouched and the temp file is removed.
    assert destination.is_dir()
    assert not [p for p in tmp_path.rglob("*.tmp")]


def test_cleanup_errors_never_override_original_failure(tmp_path, monkeypatch):
    source = write(tmp_path, "in.fa", ">a\nACGT\n")
    target = tmp_path / "out.fa"
    target.write_text("KEEP\n")

    def boom_replace(src, dst):
        raise OSError("boom-replace")

    def boom_unlink(path):
        raise OSError("cleanup-boom")

    monkeypatch.setattr(cli.os, "replace", boom_replace)
    monkeypatch.setattr(cli.os, "unlink", boom_unlink)

    code, out, err = run(
        ["normalize-sequences", source, "--output", str(target)]
    )
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1
    assert "boom-replace" in err
    assert "cleanup-boom" not in err
    assert target.read_text() == "KEEP\n"


# -- Failure in a later input stage -----------------------------------------

@pytest.mark.parametrize("name", ENTRY_IDS)
def test_later_input_failure_preserves_existing_target(name, tmp_path):
    spec = ENTRIES[name]
    late_argv, expected_code = spec.make_late_failure(tmp_path)
    target = tmp_path / f"existing{spec.ext}"
    target.write_bytes(b"PREEXISTING-BYTES\n")

    code, out, err = run(late_argv + ["--output", str(target)])
    assert code == expected_code
    assert out == ""
    assert err.count("\n") == 1 and err.endswith("\n")
    assert target.read_bytes() == b"PREEXISTING-BYTES\n"
    assert not [p for p in tmp_path.rglob("*.tmp")]
    assert not [
        p for p in tmp_path.rglob("*")
        if p.is_file() and p.name.startswith(f".{spec.command}-")
    ]


@pytest.mark.parametrize("name", sorted(LATE_OSERROR))
def test_later_input_oserror_is_exit_1_and_preserves_target(name, tmp_path):
    spec = ENTRIES[name]
    late_argv, expected_code = LATE_OSERROR[name](tmp_path)
    assert expected_code == 1
    target = tmp_path / "existing.jsonl"
    target.write_text("KEEP\n")

    code, out, err = run(late_argv + ["--output", str(target)])
    assert code == 1
    assert out == ""
    assert err.count("\n") == 1
    assert target.read_text() == "KEEP\n"
    assert not [p for p in tmp_path.rglob("*.tmp")]


# -- Empty-result representations -------------------------------------------

@pytest.mark.parametrize(
    "name",
    [n for n, s in ENTRIES.items() if s.make_empty_argv is not None],
)
def test_empty_result_is_committed_with_the_canonical_bytes(name, tmp_path):
    spec = ENTRIES[name]
    target = tmp_path / f"empty{spec.ext}"
    target.write_text("SHOULD-BE-REPLACED\n")

    argv = spec.make_empty_argv(tmp_path)
    code_std, out_std, err_std = run(list(argv))
    assert code_std == 0 and err_std == ""

    code, out, err = run(list(argv) + ["--output", str(target)])
    assert code == 0 and out == "" and err == ""
    assert target.read_bytes() == out_std.encode("utf-8")
    # Running it again is byte-stable.
    code, _, err = run(list(argv) + ["--output", str(target)])
    assert code == 0 and err == ""
    assert target.read_bytes() == out_std.encode("utf-8")
    assert not [p for p in tmp_path.rglob("*.tmp")]


# -- Standard input entries never touch the filesystem for output -----------

@pytest.mark.parametrize(
    "name",
    ["normalize-sequences", "filter-reads", "deduplicate-reads", "kmer-index"],
)
def test_stdin_to_stdout_creates_no_output_file(
    name, tmp_path, monkeypatch
):
    spec = ENTRIES[name]
    assert spec.make_stdin is not None
    argv, payload = spec.make_stdin(None)

    seen = {"replace": 0}

    def no_replace(src, dst):  # pragma: no cover - must never be called
        seen["replace"] += 1
        raise AssertionError("os.replace used for standard output")

    monkeypatch.setattr(cli.os, "replace", no_replace)
    monkeypatch.chdir(tmp_path)
    code, out, err = run(argv, payload)
    assert code == 0 and err == ""
    assert out != ""
    assert seen["replace"] == 0
    assert list(tmp_path.iterdir()) == []


# -- Helpers sanity: failure prints a single line ---------------------------

def test_commit_helper_failure_message_format(tmp_path):
    parser = cli._build_parser()
    old_err = sys.stderr
    sys.stderr = StringIO()
    try:
        code = cli._commit_output(
            "x", str(tmp_path / "no" / "o"), parser, ".x-"
        )
        message = sys.stderr.getvalue()
    finally:
        sys.stderr = old_err
    assert code == 1
    assert message.count("\n") == 1
    assert "no" in message
