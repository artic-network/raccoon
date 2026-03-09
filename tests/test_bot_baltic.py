import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from raccoon.utils import plotly_baltic


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
    monkeypatch.setattr(plotly_baltic, "load_tree", lambda *_args, **_kwargs: FakeTree())
    monkeypatch.setattr(plotly_baltic, "ensure_node_label", lambda node: node.name)
    yield


def test_build_tree_plot_includes_branch_mutations_and_labels(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    html = plotly_baltic.build_tree_plot("tree", branch_snps_path=str(csv_path))

    assert "Branch: root_tipA" in html
    assert ("10: G->A (GA)" in html) or ("10: G-\\u003eA (GA)" in html)
    assert "tipA" in html
    assert "Reset view" in html
    assert "Expand tree" in html


def test_build_tree_plot_without_branch_snps(tmp_path):
    html = plotly_baltic.build_tree_plot("tree", branch_snps_path=str(tmp_path / "missing.csv"))
    assert "Reset view" in html
    assert "Branch:" not in html


def test_build_tree_plot_color_by_includes_location_date_year_and_field_indexes(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA|Loc1|ClusterA|2024"
    fake_tree.tip_b.name = "tipB|Loc2|ClusterB|2024-02"

    original_load_tree = plotly_baltic.load_tree
    original_ensure = plotly_baltic.ensure_node_label
    try:
        plotly_baltic.load_tree = lambda *_args, **_kwargs: fake_tree
        plotly_baltic.ensure_node_label = lambda node: node.name
        html = plotly_baltic.build_tree_plot("tree", branch_snps_path=str(csv_path))
    finally:
        plotly_baltic.load_tree = original_load_tree
        plotly_baltic.ensure_node_label = original_ensure

    assert '"label":"date"' in html
    assert '"label":"year"' in html
    assert '"label":"branch_type"' not in html
    assert '"label":"label"' not in html
    assert '"label":"field_1"' not in html
    assert '"label":"field_-1"' not in html


def test_build_tree_plot_color_by_includes_field_indexes_when_tip_fields_missing(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA|Loc1|ClusterA|2024"
    fake_tree.tip_b.name = "tipB|Loc2|ClusterB|2024-02"

    original_load_tree = plotly_baltic.load_tree
    original_ensure = plotly_baltic.ensure_node_label
    try:
        plotly_baltic.load_tree = lambda *_args, **_kwargs: fake_tree
        plotly_baltic.ensure_node_label = lambda node: node.name
        html = plotly_baltic.build_tree_plot("tree", branch_snps_path=str(csv_path), tip_fields="")
    finally:
        plotly_baltic.load_tree = original_load_tree
        plotly_baltic.ensure_node_label = original_ensure

    assert '"label":"field_1"' in html
    assert '"label":"field_-1"' in html


def test_build_tree_plot_strips_braces_from_tip_fields_in_dropdown(tmp_path):
    csv_path = tmp_path / "branch_snps.csv"
    _write_branch_snps_csv(csv_path)

    fake_tree = FakeTree()
    fake_tree.tip_a.name = "tipA|Loc1|2024"
    fake_tree.tip_b.name = "tipB|Loc2|2024-02"

    original_load_tree = plotly_baltic.load_tree
    original_ensure = plotly_baltic.ensure_node_label
    try:
        plotly_baltic.load_tree = lambda *_args, **_kwargs: fake_tree
        plotly_baltic.ensure_node_label = lambda node: node.name
        html = plotly_baltic.build_tree_plot(
            "tree",
            branch_snps_path=str(csv_path),
            tip_fields="{sample}|{location}|{date}",
        )
    finally:
        plotly_baltic.load_tree = original_load_tree
        plotly_baltic.ensure_node_label = original_ensure

    assert '"label":"sample"' in html
    assert '"label":"location"' in html
    assert '"label":"date"' in html
    assert '"label":"{sample}"' not in html
    assert '"label":"{location}"' not in html
    assert '"label":"{date}"' not in html
