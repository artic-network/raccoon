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

def test_n_adjacent_snps_flagged_even_when_another_seq_has_n_at_snp_position():
    """Test that SNPs adjacent to N are flagged even if another sequence has N at the SNP position.
    
    This is a regression test for a bug where the N-adjacent check would fail if:
    - Sequence A has a SNP at position 10 with N adjacent at position 9
    - Sequence B has N at position 10 (same position as A's SNP)
    
    The bug was that N doesn't count in col_counter (only ATGC), creating a tie
    in majority calculation, which could make the SNP appear to match the majority,
    thus skipping the N-adjacent check.
    """
    ref = "AAAAAAAAAAAAAAAAAAAA"
    background = ref

    # query has unique SNP at position 10 (T) and N adjacent at position 9
    query_list = list(ref)
    query_list[10] = "T"  # unique SNP
    query_list[9] = "N"   # adjacent N within n_window=1
    query = "".join(query_list)

    # other_seq has N at position 10 (same position as query's SNP)
    # This used to cause the N-adjacent check to fail
    other_seq_list = list(ref)
    other_seq_list[10] = "N"
    other_seq = "".join(other_seq_list)

    unique, snps_near_n, _, _ = _analyze(
        [
            ("ref", ref),
            ("background", background),
            ("query", query),
            ("other_seq", other_seq),
        ],
        n_window=1,
        gap_window=1,
        snp_window=5,
        snp_count=3,
    )

    assert unique["query"] == {10}
    assert snps_near_n["query"] == {10}, "SNP at position 10 should be flagged as N-adjacent"


def test_n_adjacent_detection_handles_lowercase_n_with_uppercase_harmonization():
    ref = "aaaaaaaaaaaaaaaaaaaa"
    background = ref

    # query has lowercase unique SNP at position 10 (g) and lowercase adjacent n at position 9
    query_list = list(ref)
    query_list[10] = "g"
    query_list[9] = "n"
    query = "".join(query_list)

    # another sequence has lowercase n at the SNP position
    # to ensure this remains robust in mixed-case, ambiguous contexts
    other_seq_list = list(ref)
    other_seq_list[10] = "n"
    other_seq = "".join(other_seq_list)

    unique, snps_near_n, _, _ = _analyze(
        [
            ("ref", ref),
            ("background", background),
            ("query", query),
            ("other_seq", other_seq),
        ],
        n_window=1,
        gap_window=1,
        snp_window=5,
        snp_count=3,
    )

    assert unique["query"] == {10}
    assert snps_near_n["query"] == {10}