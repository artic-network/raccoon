import csv
import os
import types
import pytest
from pathlib import Path

from raccoon.utils import plotting


class FakeNode:
    def __init__(self, name, height, y, parent=None, branchType='leaf'):
        self.name = name
        self.height = height
        self.y = y
        self.parent = parent
        self.branchType = branchType
        self.traits = {}

    def is_leaf(self):
        return self.branchType == 'leaf'


class FakeTree:
    def __init__(self):
        # simple two-leaf tree
        root = FakeNode('root', height=0.0, y=0.0, parent=None, branchType='internal')
        a = FakeNode('A', height=1.0, y=1.0, parent=root, branchType='leaf')
        b = FakeNode('B', height=1.0, y=2.0, parent=root, branchType='leaf')
        self.Objects = [a, b, root]
        self.treeHeight = 1.0
        self.ySpan = 10.0

    def plotTree(self, ax, x_attr=None):
        # no-op for tests
        pass

    def plotPoints(self, ax, size=None, colour=None, x_attr=None):
        pass

    def addText(self, ax=None, x_attr=None, target=None, text=None):
        pass


class FakeBalticModule:
    def load_newick(self, treefile, treeType='divergence'):
        return FakeTree()

    def load_nexus(self, treefile, treeType='divergence'):
        return FakeTree()
    
    def loadNewick(self, treefile, absoluteTime=False):
        return FakeTree()

    def loadNexus(self, treefile, absoluteTime=False):
        return FakeTree()


@pytest.fixture(autouse=True)
def patch_baltic(monkeypatch):
    """Monkeypatch plotting.bt to a fake baltic module for tests."""
    fake = FakeBalticModule()
    monkeypatch.setattr(plotting, 'bt', fake)
    yield


def write_branch_snps_csv(path):
    rows = [
        {'parent': 'root', 'child': 'A', 'site': '10', 'snp': 'G->A', 'dimer': 'GA'},
        {'parent': 'root', 'child': 'B', 'site': '20', 'snp': 'C->T', 'dimer': 'TC'},
    ]
    with open(path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['parent', 'child', 'site', 'snp', 'dimer'])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def test_make_reversion_tree_figure_writes_files(tmp_path):
    base = tmp_path / 'rev_out'
    csv_path = tmp_path / 'branch_snps.csv'
    write_branch_snps_csv(csv_path)

    branch_reversions = {'root_A': {'10A'}}
    will_be_reverted = {'root_A': {'5G'}}

    plotting.make_reversion_tree_figure(str(base), str(csv_path), branch_reversions, will_be_reverted, treefile='tree', w=None, h=None)

    assert base.with_suffix('.svg').exists()
    assert base.with_suffix('.png').exists()


def test_make_reversion_tree_figure_accepts_dict(tmp_path):
    base = tmp_path / 'rev_out2'
    branch_snps_dict = {
        'root_A': [('10', 'G->A', 'GA')],
        'root_B': [('20', 'C->T', 'TC')],
    }
    plotting.make_reversion_tree_figure(str(base), branch_snps_dict, {}, {}, treefile='tree', w=10, h=4)
    assert base.with_suffix('.svg').exists()
    assert base.with_suffix('.png').exists()


def test_make_convergence_tree_figure(tmp_path):
    base = tmp_path / 'conv_out'
    branch_snps_dict = {
        'root_A': [('10', 'G->A', 'GA')],
    }
    branch_convergence = {'root_A': {'G10A'}}
    plotting.make_convergence_tree_figure(str(base), branch_snps_dict, branch_convergence, treefile='tree', w=None, h=None)
    assert base.with_suffix('.svg').exists()
    assert base.with_suffix('.png').exists()


def test_make_reconstruction_tree_figure_w_labels(tmp_path):
    base = tmp_path / 'rec_out'
    branch_snps_dict = {
        'root_A': [('10', 'G->A', 'GA')],
    }
    plotting.make_reconstruction_tree_figure_w_labels(str(base), branch_snps_dict, treefile='tree', point_style='circle', justification='right', w=10, h=5)
    assert base.with_suffix('.svg').exists()
    assert base.with_suffix('.png').exists()


def test_functions_require_baltic(monkeypatch):
    # ensure functions raise when baltic not available
    monkeypatch.setattr(plotting, 'bt', None)
    with pytest.raises(RuntimeError):
        plotting.make_convergence_tree_figure('x', {}, {}, 'tree', None, None)
    with pytest.raises(RuntimeError):
        plotting.make_reversion_tree_figure('x', {}, {}, {}, 'tree', None, None)
    with pytest.raises(RuntimeError):
        plotting.make_reconstruction_tree_figure_w_labels('x', {}, 'tree', 'circle', 'left', None, None)
