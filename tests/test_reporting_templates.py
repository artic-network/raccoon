import csv
import re
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
    assert "6. Diversity" in html
    assert "smoothed over a 5-base window" in html
    assert "theoretical range is 0 to" in html
    assert "clustered SNPs" in html


def test_generate_tree_report_renders_template(tmp_path: Path) -> None:
    treefile = tmp_path / "tree.nwk"
    treefile.write_text("(A|Loc1|2024-01-01:0.1,B|Loc2|2024-02-01:0.2);")

    flags_csv = tmp_path / "flags.csv"
    with flags_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mutation_type", "site", "mutation", "present_in", "mask_boolean"])
        writer.writeheader()
        writer.writerow({"mutation_type": "convergent", "site": "100.0", "mutation": "A>G", "present_in": "Node1_Node2", "mask_boolean": True})
        writer.writerow({"mutation_type": "reversion", "site": "200.0", "mutation": "G>A", "present_in": "Node2_Node3", "mask_boolean": False})
        writer.writerow({"mutation_type": "adar", "site": "300.0", "mutation": "TC>TT", "present_in": "Node3_Node4", "mask_boolean": True})

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=str(flags_csv),
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    assert "Raccoon tree-qc report" in html
    assert "Convergent mutations" in html
    assert "Reversions" in html
    assert "Signatures of human immune editing" in html
    assert "Primer: the slope is the estimated evolutionary rate" in html
    assert "Slope (rate, subs/site/year):" in html
    assert "tMRCA (x-intercept, decimal year):" in html
    assert "R²:" in html
    assert ">100<" in html
    assert ">200<" in html
    assert ">300<" in html
    assert "100.0" not in html
    assert "200.0" not in html
    assert "300.0" not in html
    assert "mask_boolean" not in html


def test_generate_tree_report_root_to_tip_accepts_mixed_date_precision(tmp_path: Path) -> None:
    treefile = tmp_path / "tree_mixed_dates.nwk"
    treefile.write_text("(A|Loc1|2024:0.1,B|Loc2|2024-02:0.2,C|Loc3|2024-03-15:0.25);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    assert "No root-to-tip distances available." not in html
    assert "Slope (rate, subs/site/year):" in html
    assert "tMRCA (x-intercept, decimal year):" in html


def test_generate_tree_report_supports_midpoint_root_option(tmp_path: Path) -> None:
    treefile = tmp_path / "tree_midpoint.nwk"
    treefile.write_text("((A|Loc1|2024-01-01:0.1,B|Loc2|2024-02-01:0.2):0.1,C|Loc3|2024-03-01:0.3);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        midpoint_root=True,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    assert "Raccoon tree-qc report" in html
    assert "Input tree" in html
    assert "Tree rooting: midpoint rooted" in html


def test_generate_tree_report_displays_rooting_method_outgroup(tmp_path: Path) -> None:
    treefile = tmp_path / "tree.nwk"
    treefile.write_text("(A|Loc1|2024-01-01:0.1,B|Loc2|2024-02-01:0.2);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        midpoint_root=False,
        outgroup_ids=["A"],
    )

    html = Path(report_path).read_text()
    assert "Tree rooting: outgroup rooted" in html


def test_generate_tree_report_displays_rooting_method_unknown(tmp_path: Path) -> None:
    treefile = tmp_path / "tree.nwk"
    treefile.write_text("(A|Loc1|2024-01-01:0.1,B|Loc2|2024-02-01:0.2);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        midpoint_root=False,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    assert "Tree rooting: unknown" in html


def test_generate_tree_report_root_to_tip_renders_plotly_plot(tmp_path: Path) -> None:
    """Test that root-to-tip regression plot renders as Plotly HTML."""
    treefile = tmp_path / "tree.nwk"
    treefile.write_text("(A|Loc1|2024-01-01:0.1,B|Loc2|2024-02-01:0.15,C|Loc3|2024-03-01:0.2);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    # Check for Plotly plot markers
    assert "plotly" in html.lower()
    assert "Root-to-tip regression" in html
    # Check that stats are displayed
    assert "Slope (rate, subs/site/year):" in html
    assert "tMRCA (x-intercept, decimal year):" in html
    assert "R²:" in html
    # Ensure no error message
    assert "No root-to-tip distances available." not in html


def test_generate_tree_report_root_to_tip_preserves_date_precision_in_hover(tmp_path: Path) -> None:
    """Test that date precision is preserved in hover text (year, year-month, full date)."""
    treefile = tmp_path / "tree.nwk"
    # Mix of date precisions: year, year-month, full date
    treefile.write_text("(A|Loc1|2023:0.05,B|Loc2|2024-02:0.1,C|Loc3|2024-03-15:0.15);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    # Check that different date formats are present in hover text
    # Year format should be "2023"
    assert "2023" in html
    # Month format should be "2024-02"
    assert "2024-02" in html
    # Full date format should be "2024-03-15"
    assert "2024-03-15" in html
    # Verify regression rendered
    assert "Slope (rate, subs/site/year):" in html


def test_generate_tree_report_root_to_tip_with_custom_tip_fields(tmp_path: Path) -> None:
    """Test root-to-tip with custom tip_fields where date is not the last field."""
    treefile = tmp_path / "tree.nwk"
    # Custom format: sample|date|location (date in middle)
    treefile.write_text("(A|2024-01-01|LocA:0.1,B|2024-02-01|LocB:0.15,C|2024-03-01|LocC:0.2);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        tip_fields="sample|date|location",
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    # Verify regression rendered with custom template
    assert "Slope (rate, subs/site/year):" in html
    assert "tMRCA (x-intercept, decimal year):" in html
    assert "No root-to-tip distances available." not in html
    # Check dates are extracted correctly
    assert "2024-01-01" in html
    assert "2024-02-01" in html


def test_generate_tree_report_root_to_tip_fallback_to_last_field(tmp_path: Path) -> None:
    """Test that date extraction falls back to last field if specified field is invalid."""
    treefile = tmp_path / "tree.nwk"
    # Template says date is at index 1, but it contains invalid data
    # Last field contains valid date - should fallback
    treefile.write_text("(A|NotADate|LocA|2024-01:0.1,B|NotADate|LocB|2024-02:0.15);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        tip_fields="sample|date|location|actualdate",
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    # Should fallback to last field and render successfully
    assert "Slope (rate, subs/site/year):" in html
    assert "No root-to-tip distances available." not in html


def test_generate_tree_report_root_to_tip_with_insufficient_data(tmp_path: Path) -> None:
    """Test that root-to-tip gracefully handles insufficient data (< 2 tips with dates)."""
    treefile = tmp_path / "tree.nwk"
    # Only one tip with valid date
    treefile.write_text("(A|Loc1|2024-01-01:0.1,B|Loc2|NotADate:0.2);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    # Should show error message when insufficient data
    assert "No root-to-tip distances available." in html


def test_generate_tree_report_root_to_tip_calculates_stats_correctly(tmp_path: Path) -> None:
    """Test that root-to-tip regression stats are calculated and displayed."""
    treefile = tmp_path / "tree.nwk"
    # Create tips with temporal signal
    treefile.write_text("(A|Loc1|2024-01-01:0.001,B|Loc2|2024-06-01:0.002,C|Loc3|2024-12-01:0.003);")

    report_path = generate_phylo_report(
        outdir=str(tmp_path),
        treefile=str(treefile),
        flags_csv=None,
        outgroup_ids=None,
    )

    html = Path(report_path).read_text()
    # Check that all stats labels are present
    assert "Slope (rate, subs/site/year):" in html
    assert "tMRCA (x-intercept, decimal year):" in html
    assert "R²:" in html
    # Check that regression plot is present
    assert "Root-to-tip regression" in html
    assert "plotly" in html.lower()
    # Should not show error message
    assert "No root-to-tip distances available." not in html


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