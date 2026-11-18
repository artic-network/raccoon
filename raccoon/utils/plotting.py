"""Plotting helpers for raccoon.

These functions take parsed data (or file paths) and create figures. They are
pure (accept an optional Matplotlib `ax`) and suitable for headless testing.
"""
from typing import Optional
import logging
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

# Try to import baltic lazily when plotting functions are called
try:
    import baltic as bt
except Exception:
    bt = None

plt.switch_backend('Agg')

mpl.rcParams.update({'font.size': 18})
new_rc_params = {'text.usetex': False, 'svg.fonttype': 'none'}
mpl.rcParams.update(new_rc_params)
plt.rcParams['font.family'] = 'Helvetica'


def _read_branch_snps(branch_snps):
    """Read branch snps CSV into a dict if a filepath provided."""
    import csv
    import collections

    if isinstance(branch_snps, dict):
        return branch_snps

    d = collections.defaultdict(list)
    with open(branch_snps, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            branch = f"{row['parent']}_{row['child']}"
            d[branch].append((row['site'], row['snp'], row['dimer']))
    return d


def make_reversion_tree_figure(outfile: str, branch_snps, branch_reversions, will_be_reverted, treefile: str, w: Optional[float], h: Optional[float], ax: Optional[plt.Axes] = None):
    """Plot tree figure showing reversions.

    `branch_snps` can be either a path to a CSV or a pre-parsed dict of branch->snps.
    """
    if bt is None:
        raise RuntimeError('baltic is required for plotting trees')

    branch_snps_dict = _read_branch_snps(branch_snps)

    my_tree = bt.load_newick(treefile, treeType='divergence')

    # heuristic sizing
    r2t = 200000 * my_tree.treeHeight
    if w is None:
        width = int(math.sqrt(r2t) * 3) if r2t < 200 else 25
    else:
        width = float(w)
    if h is None:
        height = int(math.sqrt(my_tree.ySpan) * 2) if my_tree.ySpan < 300 else 40
    else:
        height = float(h)

    if ax is None:
        fig, ax = plt.subplots(figsize=(width, height), facecolor='w')
    else:
        fig = ax.figure

    increment = my_tree.treeHeight / 150

    x_attr = lambda k: k.height
    su_func = lambda k: 50 - 30 * k.height / my_tree.treeHeight
    s_func = lambda k: 50 - 20 * k.height / my_tree.treeHeight
    c_func = lambda k: 'dimgrey'

    my_tree.plotTree(ax, x_attr=x_attr)
    my_tree.plotPoints(ax, size=s_func, colour=c_func, x_attr=x_attr)
    mpl.rcParams['font.family'] = 'sans-serif'

    text_x_attr = lambda k: k.height + (increment * 4)
    my_tree.addText(ax, x_attr=text_x_attr, target=lambda k: k.is_leaf(), text=lambda k: k.name)

    for k in my_tree.Objects:
        current_node = k
        if k.branchType == 'leaf':
            current_node.traits['label'] = k.name

        try:
            parent_name = current_node.parent.traits['label']
        except Exception:
            continue
        node_name = current_node.traits['label']
        branch_name = f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            snps = []
            reversions = []
            tb_reversions = []
            if branch_name in branch_reversions:
                for s in branch_reversions[branch_name]:
                    reversions.append(s)
            if branch_name in will_be_reverted:
                for s in will_be_reverted[branch_name]:
                    tb_reversions.append(s)

            snp_placement = current_node.parent.height + increment
            rev_placement = (current_node.parent.height + current_node.height) / 2
            tb_rev_placement = (current_node.parent.height + current_node.height) / 2

            for s in branch_snps_dict[branch_name]:
                site, snp, dimer = s
                if snp == 'G->A':
                    if dimer in ['GA']:
                        snps.append((1, '#995e62'))
                    else:
                        snps.append((2, '#d9b660'))
                elif snp == 'C->T':
                    if dimer in ['TC']:
                        snps.append((1, '#995e62'))
                    else:
                        snps.append((2, '#d9b660'))
                else:
                    snps.append((2, '#d9b660'))

            for reversion in reversions:
                ax.text(s=reversion, x=rev_placement, y=k.y + 0.5, rotation='vertical')
                ax.scatter([rev_placement], [k.y], color='black', s=150, marker=8)
                rev_placement += 3 * increment
            for tb_reversion in tb_reversions:
                ax.text(s=tb_reversion, x=tb_rev_placement, y=k.y + 0.5, rotation='vertical')
                ax.scatter([tb_rev_placement], [k.y], color='black', s=150, marker=9)
                tb_rev_placement += 3 * increment

            for snp in sorted(snps, key=lambda x: x[0]):
                ax.scatter([snp_placement], [k.y + 0.5], color=snp[1], s=30)
                snp_placement += increment

    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left', 'bottom']]
    ax.tick_params(axis='y', size=0)
    ax.tick_params(axis='x', size=0)

    ax.set_yticklabels([])
    ax.set_xticklabels([])

    plt.savefig(f"{outfile}.svg", bbox_inches='tight')
    plt.savefig(f"{outfile}.png", bbox_inches='tight', transparent=True)


def make_convergence_tree_figure(outfile: str, branch_snps, branch_convergence, treefile: str, w: Optional[float], h: Optional[float], ax: Optional[plt.Axes] = None):
    """Plot tree figure showing convergent SNPs."""
    if bt is None:
        raise RuntimeError('baltic is required for plotting trees')

    branch_snps_dict = _read_branch_snps(branch_snps)

    my_tree = bt.load_newick(treefile, treeType='divergence')

    r2t = 200000 * my_tree.treeHeight
    if w is None:
        width = int(math.sqrt(r2t) * 3) if r2t < 200 else 25
    else:
        width = float(w)
    if h is None:
        height = int(math.sqrt(my_tree.ySpan) * 2) if my_tree.ySpan < 300 else 40
    else:
        height = float(h)

    if ax is None:
        fig, ax = plt.subplots(figsize=(width, height), facecolor='w')
    else:
        fig = ax.figure

    increment = my_tree.treeHeight / 150

    x_attr = lambda k: k.height
    s_func = lambda k: 50 - 20 * k.height / my_tree.treeHeight
    c_func = lambda k: 'dimgrey'

    my_tree.plotTree(ax, x_attr=x_attr)
    my_tree.plotPoints(ax, size=s_func, colour=c_func, x_attr=x_attr)
    mpl.rcParams['font.family'] = 'sans-serif'

    text_x_attr = lambda k: k.height + (increment * 4)
    my_tree.addText(ax, x_attr=text_x_attr, target=lambda k: k.is_leaf(), text=lambda k: k.name)

    for k in my_tree.Objects:
        current_node = k
        if k.branchType == 'leaf':
            current_node.traits['label'] = k.name

        try:
            parent_name = current_node.parent.traits['label']
        except Exception:
            continue
        node_name = current_node.traits['label']
        branch_name = f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            snps = []
            convergent_snps = []
            if branch_name in branch_convergence:
                for s in branch_convergence[branch_name]:
                    convergent_snps.append(s)

            snp_placement = current_node.parent.height + increment
            c_placement = (current_node.parent.height + current_node.height) / 2
            for s in branch_snps_dict[branch_name]:
                site, snp, dimer = s
                if snp == 'G->A':
                    if dimer in ['GA']:
                        snps.append((1, '#995e62'))
                    else:
                        snps.append((2, '#d9b660'))
                elif snp == 'C->T':
                    if dimer in ['TC']:
                        snps.append((1, '#995e62'))
                    else:
                        snps.append((2, '#d9b660'))
                else:
                    snps.append((2, '#d9b660'))

            for c_snp in convergent_snps:
                ax.text(s=c_snp, x=c_placement, y=k.y + 0.5, rotation='vertical')
                ax.scatter([c_placement], [k.y], color='black', s=150, marker='d')
                c_placement += 2 * increment

            for snp in sorted(snps, key=lambda x: x[0]):
                ax.scatter([snp_placement], [k.y + 0.5], color=snp[1], s=30)
                snp_placement += increment

    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left', 'bottom']]
    ax.tick_params(axis='y', size=0)
    ax.tick_params(axis='x', size=0)

    ax.set_yticklabels([])
    ax.set_xticklabels([])

    plt.savefig(f"{outfile}.svg", bbox_inches='tight')
    plt.savefig(f"{outfile}.png", bbox_inches='tight', transparent=True)


def make_reconstruction_tree_figure_w_labels(outfile: str, branch_snps, treefile: str, point_style: str, justification: str, w: Optional[float] = None, h: Optional[float] = None, ax: Optional[plt.Axes] = None):
    """Plot reconstruction tree figure with labelled SNPs."""
    if bt is None:
        raise RuntimeError('baltic is required for plotting trees')

    branch_snps_dict = _read_branch_snps(branch_snps)

    my_tree = bt.load_newick(treefile, treeType='divergence')

    if w is None:
        r2t = 200000 * my_tree.treeHeight
        width = int(math.sqrt(r2t) * 3) if r2t < 200 else 25
    else:
        width = int(w)

    if h is None:
        height = int(math.sqrt(my_tree.ySpan) * 2) if my_tree.ySpan < 300 else 40
    else:
        height = int(h)

    if ax is None:
        fig, ax = plt.subplots(figsize=(width, height), facecolor='w')
    else:
        fig = ax.figure

    x_attr = lambda k: k.height
    s_func = lambda k: 50 - 20 * k.height / my_tree.treeHeight
    c_func = lambda k: 'dimgrey'

    my_tree.plotTree(ax, x_attr=x_attr)
    my_tree.plotPoints(ax, size=s_func, colour=c_func, x_attr=x_attr)
    mpl.rcParams['font.family'] = 'sans-serif'

    text_x_attr = lambda k: k.height + (my_tree.treeHeight / 150 * 4)
    my_tree.addText(ax, x_attr=text_x_attr, target=lambda k: k.is_leaf(), text=lambda k: k.name)

    right_settings = {"apobec": (1, "#995E62"), "non_apobec": (2, "#D9B660")}
    left_settings = {"apobec": (2, "#995E62"), "non_apobec": (1, "#D9B660")}

    for k in my_tree.Objects:
        if k.branchType == 'leaf':
            k.traits['label'] = k.name

        node_name = k.traits.get('label')
        if not node_name:
            continue
        
        try:
            parent_name = k.parent.traits.get('label')
            if not parent_name:
                continue
        except Exception:
            continue
        branch_name = f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            if justification == 'right':
                setting_dict = right_settings
                snp_placement = k.height - my_tree.treeHeight / 150
            else:
                setting_dict = left_settings
                snp_placement = k.parent.height + (my_tree.treeHeight / 150) / 2

            snps = []
            for s in branch_snps_dict[branch_name]:
                site, snp, dimer = s
                if snp == 'G->A':
                    if dimer in ['GA']:
                        snps.append(setting_dict['apobec'])
                    else:
                        snps.append(setting_dict['non_apobec'])
                elif snp == 'C->T':
                    if dimer in ['TC']:
                        snps.append(setting_dict['apobec'])
                    else:
                        snps.append(setting_dict['non_apobec'])
                else:
                    snps.append(setting_dict['non_apobec'])

            for snp in sorted(snps, key=lambda x: x[0]):
                if point_style == 'circle':
                    ax.scatter([snp_placement], [k.y + 0.5], color=snp[1], s=50)
                else:
                    import matplotlib.patches as patches

                    rect = patches.Rectangle((snp_placement, k.y - 0.5), my_tree.treeHeight / 150 / 2, 1, alpha=1, fill=True, edgecolor='none', facecolor=snp[1])
                    ax.add_patch(rect)

                if justification == 'right':
                    snp_placement -= my_tree.treeHeight / 150
                else:
                    snp_placement += my_tree.treeHeight / 150

    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    ax.tick_params(axis='y', size=0)
    ax.tick_params(axis='x', size=0)

    ax.set_yticklabels([])
    ax.set_xticklabels([])

    plt.margins(0.1, 0.1, tight=True)
    plt.savefig(f"{outfile}.svg", bbox_inches='tight')
    plt.savefig(f"{outfile}.png", bbox_inches='tight')