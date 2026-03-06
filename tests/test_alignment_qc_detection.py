"""Focused tests for alignment QC mutation detection rules."""

import tempfile
from pathlib import Path

from raccoon.utils import alignment_functions as af


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n{sequence}\n")


def _analyze(records: list[tuple[str, str]], **kwargs):
    lengths = {len(sequence) for _, sequence in records}
    assert len(lengths) == 1, "All sequences in an alignment must have identical length"

    with tempfile.TemporaryDirectory() as tmpdir:
        alignment_path = Path(tmpdir) / "test_alignment.fasta"
        _write_fasta(alignment_path, records)
        return af.analyze_alignment(str(alignment_path), **kwargs)


def test_clustered_snps_triggers_when_threshold_met():
    ref = "AAAAAAAAAAAAAAAAAAAA"
    background = ref

    query_list = list(ref)
    for pos in (4, 6, 8):
        query_list[pos] = "T"
    query = "".join(query_list)

    unique, _, _, clustered = _analyze(
        [
            ("ref", ref),
            ("background", background),
            ("query", query),
        ],
        snp_window=5,
        snp_count=3,
    )

    assert unique["query"] == {4, 6, 8}
    assert clustered["query"] == {4, 6, 8}
    assert clustered.get("ref", set()) == set()
    assert clustered.get("background", set()) == set()


def test_clustered_snps_does_not_trigger_when_threshold_not_met():
    ref = "AAAAAAAAAAAAAAAAAAAA"
    background = ref

    query_list = list(ref)
    for pos in (2, 8, 14):
        query_list[pos] = "T"
    query = "".join(query_list)

    unique, _, _, clustered = _analyze(
        [
            ("ref", ref),
            ("background", background),
            ("query", query),
        ],
        snp_window=5,
        snp_count=3,
    )

    assert unique["query"] == {2, 8, 14}
    assert clustered.get("query", set()) == set()


def test_n_adjacent_snps_are_flagged():
    ref = "AAAAAAAAAAAAAAAAAAAA"
    background = ref

    query_list = list(ref)
    query_list[10] = "T"  # unique SNP
    query_list[9] = "N"   # adjacent N within n_window=1
    query = "".join(query_list)

    unique, snps_near_n, _, _ = _analyze(
        [
            ("ref", ref),
            ("background", background),
            ("query", query),
        ],
        n_window=1,
        gap_window=1,
        snp_window=5,
        snp_count=3,
    )

    assert unique["query"] == {10}
    assert snps_near_n["query"] == {10}


def test_gap_adjacent_snps_are_flagged():
    ref = "AAAAAAAAAAAAAAAAAAAA"
    background = ref

    query_list = list(ref)
    query_list[10] = "T"  # unique SNP
    query_list[9] = "-"   # adjacent gap within gap_window=1
    query = "".join(query_list)

    unique, _, snps_near_gap, _ = _analyze(
        [
            ("ref", ref),
            ("background", background),
            ("query", query),
        ],
        n_window=1,
        gap_window=1,
        snp_window=5,
        snp_count=3,
    )

    assert unique["query"] == {10}
    assert snps_near_gap["query"] == {10}
