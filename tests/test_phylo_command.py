from raccoon.commands import phylo as phylo_cmd


class _Node:
    def __init__(self, label: str, branch_type: str = "leaf"):
        self.name = label
        self.branchType = branch_type


class _Tree:
    def __init__(self, labels):
        self.Objects = [_Node(label) for label in labels]


def test_validate_tip_label_fields_passes_when_labels_match_template(monkeypatch):
    tree = _Tree(["sample1|locA|2024-01-01", "sample2|locB|2024-02"])

    def fake_load_tree(_treefile, tree_format="auto"):
        return tree

    def fake_ensure_label(node):
        return node.name

    monkeypatch.setattr("raccoon.utils.reconstruction_functions.load_tree", fake_load_tree)
    monkeypatch.setattr("raccoon.utils.reconstruction_functions.ensure_node_label", fake_ensure_label)

    ok, err = phylo_cmd._validate_tip_label_fields("tree.nwk", "auto", "sample|location|date")
    assert ok is True
    assert err is None


def test_validate_tip_label_fields_fails_when_labels_too_short(monkeypatch):
    tree = _Tree(["sample1|locA", "sample2|locB"])

    def fake_load_tree(_treefile, tree_format="auto"):
        return tree

    def fake_ensure_label(node):
        return node.name

    monkeypatch.setattr("raccoon.utils.reconstruction_functions.load_tree", fake_load_tree)
    monkeypatch.setattr("raccoon.utils.reconstruction_functions.ensure_node_label", fake_ensure_label)

    ok, err = phylo_cmd._validate_tip_label_fields("tree.nwk", "auto", "sample|location|date")
    assert ok is False
    assert "expected >= 3 fields" in err


def test_validate_tip_label_fields_requires_non_empty_template():
    ok, err = phylo_cmd._validate_tip_label_fields("tree.nwk", "auto", "")
    assert ok is False
    assert "must define at least one field name" in err
