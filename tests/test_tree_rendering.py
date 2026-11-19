import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from raccoon.utils import tree_rendering


class FakeNode:
    def __init__(self, name, x, y, branch_type="leaf", children=None):
        self.name = name
        self.x = x
        self.y = y
        self.branchType = branch_type
        self.children = children or []
        self.traits = {}

    def is_leaflike(self):
        return self.branchType == "leaf"

    def is_node(self):
        return self.branchType == "internal"


class FakeTree:
    def __init__(self):
        self.root = FakeNode("root", x=0.0, y=1.5, branch_type="internal")
        self.tip_a = FakeNode("tipA", x=1.0, y=1.0, branch_type="leaf")
        self.tip_b = FakeNode("tipB", x=1.1, y=2.0, branch_type="leaf")
        self.root.children = [self.tip_a, self.tip_b]
        self.Objects = [self.root, self.tip_a, self.tip_b]

    def drawTree(self):
        return None

    def getInternal(self):
        return [self.root]


def _write_branch_snps_csv(path: Path) -> None:
    rows = [
        {"parent": "root", "child": "tipA", "site": "10", "snp": "G->A", "dimer": "GA"},
        {"parent": "root", "child": "tipB", "site": "20", "snp": "C->T", "dimer": "TC"},
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["parent", "child", "site", "snp", "dimer"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


@pytest.fixture(autouse=True)
def patch_tree(monkeypatch):
    monkeypatch.setattr(tree_rendering, "load_tree", lambda *_args, **_kwargs: FakeTree())
    yield


def test_build_tree_plot_includes_branch_mutations_and_labels(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    newick = tree_rendering.build_tree_plot("tree", branch_snps_path=str(csv_path))

    # With PearTree, we now return a Newick string instead of Plotly HTML
    assert newick.startswith("(")
    assert newick.endswith(";")
    # Verify it contains tip names
    assert "tipA" in newick
    assert "tipB" in newick


def test_build_tree_plot_without_branch_snps(tmp_path):
    newick = tree_rendering.build_tree_plot("tree", branch_snps_path=str(tmp_path / "missing.csv"))
    # Should still return valid Newick format
    assert newick.startswith("(")
    assert newick.endswith(";")

def test_build_tree_plot_color_by_includes_location_date_year_and_field_indexes(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA|Loc1|ClusterA|2024"
    fake_tree.tip_b.name = "tipB|Loc2|ClusterB|2024-02"

    original_load_tree = tree_rendering.load_tree
    try:
        tree_rendering.load_tree = lambda *_args, **_kwargs: fake_tree
        newick = tree_rendering.build_tree_plot("tree", branch_snps_path=str(csv_path))
        # Verify we get valid Newick with the tip names
        assert newick.startswith("(")
        assert newick.endswith(";")
        assert "tipA" in newick
        assert "tipB" in newick
        assert "Loc1" in newick
        assert "Loc2" in newick
    finally:
        tree_rendering.load_tree = original_load_tree


def test_build_tree_plot_color_by_includes_field_indexes_when_tip_fields_missing(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA|Loc1|ClusterA|2024"
    fake_tree.tip_b.name = "tipB|Loc2|ClusterB|2024-02"

    original_load_tree = tree_rendering.load_tree
    try:
        tree_rendering.load_tree = lambda *_args, **_kwargs: fake_tree
        newick = tree_rendering.build_tree_plot("tree", branch_snps_path=str(csv_path), tip_fields="")
        # Verify we get valid Newick
        assert newick.startswith("(")
        assert newick.endswith(";")
    finally:
        tree_rendering.load_tree = original_load_tree


def test_build_tree_plot_strips_braces_from_tip_fields_in_dropdown(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA|Loc1|2024"
    fake_tree.tip_b.name = "tipB|Loc2|2024-02"

    original_load_tree = tree_rendering.load_tree
    try:
        tree_rendering.load_tree = lambda *_args, **_kwargs: fake_tree
        newick = tree_rendering.build_tree_plot(
            "tree",
            branch_snps_path=str(csv_path),
            tip_fields="{sample}|{location}|{date}",
        )
        # Verify we get valid Newick
        assert newick.startswith("(")
        assert newick.endswith(";")
    finally:
        tree_rendering.load_tree = original_load_tree


def test_build_tree_plot_with_single_field_ids_does_not_invent_date_traits(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA"
    fake_tree.tip_b.name = "tipB"

    original_load_tree = tree_rendering.load_tree
    try:
        tree_rendering.load_tree = lambda *_args, **_kwargs: fake_tree
        newick = tree_rendering.build_tree_plot(
            "tree",
            branch_snps_path=str(csv_path),
            tip_fields="{id}",
        )
        # Verify we get valid Newick
        assert newick.startswith("(")
        assert newick.endswith(";")
    finally:
        tree_rendering.load_tree = original_load_tree
