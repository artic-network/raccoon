"""PearTree-based interactive tree rendering utilities.
This is used for tree visualization in the phylo-qc report.
PearTree is a modern phylogenetic tree viewer maintained by the ARTIC Network.
"""
from __future__ import annotations

from typing import Optional

from .reconstruction_functions import load_tree


def build_tree_plot(
    treefile: str,
    tree_format: str = "auto",
    branch_snps_path: Optional[str] = None,
    tip_fields: str = "sample|location|date",
) -> str:
    """Build a tree plot for embedding in HTML reports using PearTree.
    
    Args:
        treefile: Path to tree file (Newick or NEXUS format)
        tree_format: Tree format ('auto', 'newick', or 'nexus')
        branch_snps_path: Optional path to branch SNPs CSV (not used by PearTree)
        tip_fields: Optional tip field names (not used by PearTree)
    
    Returns:
        Newick format tree string suitable for PearTree embedding, or empty string on error
    """
    try:
        tree = load_tree(treefile, tree_format=tree_format)
    except Exception:
        return ""
    
    # Export tree as Newick format for PearTree
    try:
        return _export_newick(tree)
    except Exception:
        return ""


def _export_newick(tree) -> str:
    """Export a Baltic tree object to Newick format.
    
    Args:
        tree: A baltic.tree.Tree object
    
    Returns:
        Newick format string
    """
    def _node_to_newick(node) -> str:
        """Recursively convert a node to Newick format."""
        if getattr(node, "is_leaflike", lambda: False)():
            # Leaf node
            name = getattr(node, "name", "") or ""
            branch_length = getattr(node, "branchLength", 0.0) or 0.0
            return f"{name}:{branch_length}"
        else:
            # Internal node
            children = getattr(node, "children", []) or []
            if not children:
                # Leaf node without is_leaflike method
                name = getattr(node, "name", "") or ""
                branch_length = getattr(node, "branchLength", 0.0) or 0.0
                return f"{name}:{branch_length}"
            
            child_strings = [_node_to_newick(child) for child in children]
            inner = ",".join(child_strings)
            name = getattr(node, "name", "") or ""
            branch_length = getattr(node, "branchLength", 0.0) or 0.0
            
            if name:
                return f"({inner}){name}:{branch_length}"
            else:
                return f"({inner}):{branch_length}"
    
    # Get root and convert to Newick
    root = getattr(tree, "root", None)
    if root is None:
        return ""
    
    newick = _node_to_newick(root) + ";"
    return newick

