import csv
from pathlib import Path

from raccoon.utils.reporting import (
    generate_alignment_report,
    generate_combine_report,
    generate_mask_report,
    generate_phylo_report,
)


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for rec_id, seq in records:
            handle.write(f">{rec_id}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def test_generate_combine_report_renders_template(tmp_path: Path) -> None:
    input_a = tmp_path / "a.fasta"
    input_b = tmp_path / "b.fasta"
    input_c = tmp_path / "c.fasta"
    _write_fasta(input_a, [("seqA|Loc1|2024-01-01", "ATGCGTNNNN")])
    _write_fasta(input_b, [("seqB|Loc2|2024-01-02", "ATGCGTATGC")])
    _write_fasta(input_c, [("seqC|Loc1|2024-01-03", "ATGC")])

    output_fasta = tmp_path / "combined.fasta"
    _write_fasta(output_fasta, [("seqA|Loc1|2024-01-01", "ATGCGTNNNN")])

    metadata_csv = tmp_path / "metadata.csv"
    with metadata_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "location", "date"])
        writer.writerow(["seqA", "Loc1", "2024-01-01"])
        writer.writerow(["seqB", "Loc2", "2024-01-02"])

    report_path = generate_combine_report(
        outdir=str(tmp_path),
        output_fasta=str(output_fasta),
        input_fastas=[str(input_a), str(input_b), str(input_c)],
        metadata_paths=[str(metadata_csv)],
        min_length=5,
        max_n_content=0.4,
    )

    html = Path(report_path).read_text()
    assert "Raccoon seq-qc report" in html
    assert "Input files" in html
    assert "Final dataset" in html
    assert "Platform" in html


def test_generate_alignment_report_renders_template(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    _write_fasta(alignment, [
        ("seq1", "ATGNNN"),
        ("seq2", "ATG---"),
    ])

    mask_csv = tmp_path / "mask.csv"
    with mask_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["flagged", "type", "Minimum", "Maximum", "Length", "present_in", "note"])
        writer.writeheader()
        writer.writerow({
            "flagged": "2",
            "type": "site",
            "Minimum": "2",
            "Maximum": "2",
            "Length": "1",
            "present_in": "seq1,seq2",
            "note": "homoplasy",
        })

    report_path = generate_alignment_report(
        outdir=str(tmp_path),
        alignment_path=str(alignment),
        mask_file=str(mask_csv),
    )

    html = Path(report_path).read_text()
    assert "Raccoon aln-qc report" in html
    assert "Alignment N-content" in html
    assert "Flagged sites" in html
    assert "clustered SNPs" in html


def test_generate_tree_report_renders_template(tmp_path: Path) -> None:
    treefile = tmp_path / "tree.nwk"
    treefile.write_text("(A|Loc1|2024-01-01:0.1,B|Loc2|2024-02-01:0.2);")

    flags_csv = tmp_path / "flags.csv"
    with flags_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mutation_type", "site", "mutation"])
        writer.writeheader()
        writer.writerow({"mutation_type": "convergent", "site": "100", "mutation": "A>G"})
        writer.writerow({"mutation_type": "reversion", "site": "200", "mutation": "G>A"})

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=str(flags_csv),
    )

    html = Path(report_path).read_text()
    assert "Raccoon tree-qc report" in html
    assert "Convergent mutations" in html
    assert "Reversions" in html


def test_generate_mask_report_renders_template(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    _write_fasta(alignment, [
        ("seq1", "ATGNNN"),
        ("seq2", "ATG---"),
    ])

    mask_csv = tmp_path / "mask.csv"
    with mask_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["flagged", "type", "Minimum", "Maximum", "Length", "present_in", "note"])
        writer.writeheader()
        writer.writerow({
            "flagged": "2",
            "type": "site",
            "Minimum": "2",
            "Maximum": "2",
            "Length": "1",
            "present_in": "seq1",
            "note": "mask",
        })

    report_path = generate_mask_report(
        outdir=str(tmp_path),
        alignment_path=str(alignment),
        mask_file=str(mask_csv),
        output_alignment=str(tmp_path / "alignment.masked.fasta"),
    )

    html = Path(report_path).read_text()
    assert "Raccoon mask report" in html
    assert "Masked sites by sequence" in html
    assert "Mask file entries" in html