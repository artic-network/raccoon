import os
import logging
try:
    from raccoon.utils.misc import green, cyan
except Exception:
    def green(s):
        return s
    def cyan(s):
        return s
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import date
import datetime as dt

import csv

import baltic as bt
from scipy import stats

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import collections
import pandas as pd
import matplotlib.patches as patches

import math
plt.switch_backend('Agg') 


mpl.rcParams.update({'font.size': 18})

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)
plt.rcParams['font.family'] = 'Helvetica'


def get_node_states_all_sites(state_file,alignment):
    
    #returns a dict keys off 1-based positions in the genome
    #and the value is a list of tuples (id, base) at a given node
    # allows you to look up for a given site what the base is for a
    #given internal node or tip
    
    node_states = collections.defaultdict(list)
    c = 0
    
    ## first the reconstructed nodes
    with open(state_file,"r") as f:
        for l in f:

            if not l.startswith("#"):
                c+=1
                try:
                    node,site,state,probA,probC,probG,probT = l.rstrip("\n").split("\t")
                except Exception:
                    logging.debug(l)
                    break
                if node != "Node":
                    if state not in ["N","-"]:
                        node_states[site].append((node,state))
                    else:
                        node_states[site].append((node,""))
    ## now the tips
    for record in SeqIO.parse(alignment,"fasta"):
        for site in node_states:
            index = int(site)-1
            base = record.seq[index].upper()
            if base in ["T","C","A","G"]:
                node_states[site].append((record.id,base))
            else:
                node_states[site].append((record.id,""))
                
    return node_states

def get_header_str(dict_values):
    header_str = ""
    for i in sorted(dict_values, key = lambda i : i[0]):
        header_str += f"{i[0]},"
    header_str = header_str.rstrip(",")
    return header_str
    
    
def find_what_sites_vary_unambiguously(node_states,state_differences):
    header_str = get_header_str(node_states["1"])
    
    with open(state_differences,"w") as fw:
        fw.write(f"site,{header_str}\n")

        for site in node_states:
            info = node_states[site]
            
            # get the set of unique bases at a given site
            count = set([i[1] for i in info if i[1]])
            
            #if there's more than one
            if len(count)>1:
                
                #needs to be kep consistent with header str
                info = sorted(info, key = lambda i : i[0])
                base_str = ""
                for i in info:
                    base_str += f"{i[1]},"
                    
                base_str = base_str.rstrip(",")
                fw.write(f"{site},{base_str}\n")
    
def load_unambiguous_varying_sites(infile):
    node_states_diff = collections.defaultdict(dict)
    with open(infile,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site = row["site"]
            for col in row:
                if col != "site":
                    node_states_diff[row["site"]][col] = row[col]
    return node_states_diff


def detect_tree_format(treefile, tree_format="auto"):
    if tree_format and tree_format != "auto":
        return tree_format
    lower = treefile.lower()
    if lower.endswith(".nex") or lower.endswith(".nexus"):
        return "nexus"
    try:
        with open(treefile, "r") as handle:
            head = handle.read(200).lower()
        if "#nexus" in head or "begin trees;" in head:
            return "nexus"
    except Exception:
        pass
    return "newick"


def load_tree(treefile, tree_format="auto"):
    fmt = detect_tree_format(treefile, tree_format)
    if fmt == "nexus":
        return bt.load_nexus(treefile, treeType='divergence')
    return bt.load_newick(treefile, treeType='divergence')


def ensure_node_label(node):
    label = node.traits.get("label") if hasattr(node, "traits") else None
    if not label:
        label = getattr(node, "name", None)
    if label and "/" in label:
        label = label.split("/", 1)[0]
    if label:
        node.traits["label"] = label
    return label


def map_site_changes_to_branches(treefile, outfile,node_states,node_states_diff, tree_format="auto"): 
    my_tree=load_tree(treefile, tree_format=tree_format)
    last_node = ""
    current_node = ""

    with open(outfile,"w") as fw:
        fw.write("parent,child,site,snp,dimer\n")

        for k in my_tree.Objects:
            current_node = k
            ensure_node_label(current_node)

            if last_node:
                node_name = current_node.traits.get("label")
                parent_name = current_node.parent.traits.get("label")
                if not node_name or not parent_name:
                    last_node = current_node
                    continue
                snps = []
                for site in node_states_diff:
                    if node_name not in node_states_diff[site] or parent_name not in node_states_diff[site]:
                        continue
                    node_base = node_states_diff[site][node_name]
                    parent_base = node_states_diff[site][parent_name]

                    if node_base != parent_base:
                        if node_base in ["A","C","G","T"] and parent_base in ["A","C","G","T"]:
                            snp = f"{parent_base}->{node_base}"
                            snps.append(snp)
                            if snp == "G->A":
                                dimer_site = f"{int(site)+1}"
                                dimer_base = ""

                                for i in node_states[dimer_site]:
                                    if i[0] == parent_name:
                                        dimer_base = i[1]
                                dimer = f"{parent_base}{dimer_base}"
                            elif snp == "C->T":
                                dimer_site = f"{int(site)-1}"
                                dimer_base = ""

                                for i in node_states[dimer_site]:
                                    if i[0] == parent_name:
                                        dimer_base = i[1]
                                dimer = f"{dimer_base}{parent_base}"
                            else:
                                dimer = ""
                            fw.write(f"{parent_name},{node_name},{site},{snp},{dimer}\n")

            last_node = current_node

def read_in_branch_snps(branch_snps):
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            branch = f"{row['parent']}_{row['child']}"
            branch_snps_dict[branch].append((row['site'],row['snp'],row['dimer'])) 
    return branch_snps_dict

def get_branch_snps_sites(branch_snps):
    all_snps = collections.Counter()
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_snps[int(row["site"])] += 1

            branch_snps_dict[int(row['site'])].append([row['parent'],row['child'],row['snp'],row['dimer']])
    
    homoplasies = {}
    for k in all_snps:
        if all_snps[k] > 1:
            homoplasies[k] = all_snps[k]
            
    # print(len(homoplasies))
    # print(homoplasies)
    return branch_snps_dict,homoplasies
    
def get_acc_to_metadata_map(metadata):
    acc_dict = {}
    with open(metadata,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                acc_dict[row["accession"]] = row
            except:
                try:
                    acc_dict[row["\ufeffaccession"]] = row
                except:
                    pass
                
    return acc_dict


def get_fig_height(alignment):
    seqs = SeqIO.index(alignment, "fasta")
    n_tips = len(seqs)
    max_label_len = 0
    for rec_id in seqs:
        if len(rec_id) > max_label_len:
            max_label_len = len(rec_id)

    # heuristic to keep labels legible without excessive height
    height = (n_tips * 0.25) + (max_label_len * 0.02)
    if height < 8:
        return 8
    if height > 60:
        return 60
    return height


def get_fig_size(treefile, tree_format="auto", alignment=None):
    """Estimate figure (width, height) for tree rendering.

    Width scales with genetic distance and label length; height scales with tip count.
    """
    my_tree = load_tree(treefile, tree_format=tree_format)

    tips = []
    for node in my_tree.Objects:
        if node.branchType == 'leaf':
            label = ensure_node_label(node)
            if label:
                tips.append(label)

    n_tips = len(tips)
    max_label_len = max([len(t) for t in tips], default=0)

    if alignment:
        try:
            seqs = SeqIO.index(alignment, "fasta")
            n_tips = max(n_tips, len(seqs))
            max_label_len = max(max_label_len, max((len(rec_id) for rec_id in seqs), default=0))
        except Exception:
            pass

    # height: driven by number of tips (vertical spacing)
    height = (n_tips * 0.25) + (max_label_len * 0.01)
    height = min(max(height, 8), 60)

    # width: driven by genetic distance (tree height) + label space
    tree_height = getattr(my_tree, "treeHeight", None)
    if tree_height is None:
        width = 25
    else:
        r2t = 200000 * tree_height
        width = (math.sqrt(r2t) * 3) if r2t < 200 else 25

    label_space = (max_label_len * 0.12) + 2
    width = max(width, label_space + (height * 0.1))
    width = min(max(width, 10), 80)

    return width, height


def make_reconstruction_tree_figure_w_labels(outfile,branch_snps,treefile,point_style,justification,w=None,h=None):
    from .plotting import make_reconstruction_tree_figure_w_labels as _plot
    return _plot(outfile, branch_snps, treefile, point_style, justification, w, h)
    
def generate_reconstruction_files(alignment, state_out, state_differences):
    
    node_states = get_node_states_all_sites(state_out,alignment)
        
    find_what_sites_vary_unambiguously(node_states,state_differences)

    return node_states
    
def load_info(directory, alignment, treefile, state_out, state_differences, branch_snps_out, treefigureout,point_style,point_justify, node_states="",width=None,height=None):
    
    if not node_states:
        node_states = get_node_states_all_sites(state_out, alignment)

    node_states_diff = load_unambiguous_varying_sites(state_differences)

    map_site_changes_to_branches(treefile,
                                 branch_snps_out,
                                 node_states,
                                 node_states_diff)

    make_reconstruction_tree_figure_w_labels(treefigureout,
                                    branch_snps_out,
                                    treefile,
                                    point_style,
                                    point_justify,
                                    width,
                                    height)
    
def get_gene_boundaries(gene_boundaries_file):
    genes = {}
    gene_id = 0
    with open(gene_boundaries_file,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id +=1
            name = f"{row['Name'].replace(' ','_')}_{gene_id}"
            start = int(row["Minimum"])
            end = int(row["Maximum"]) + 1
            length = int(row["Length"])
            direction = row["Direction"]
            genes[(start,end)] = (name,length,direction)
    return genes

def get_grantham_scores(grantham_scores_file):
    grantham_scores = {}

    with open(grantham_scores_file,"r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for col in row:
                if col!="FIRST":
                    mutation = f"{row['FIRST']}{col}"

                    if row[col] != "0":
                        grantham_scores[mutation] = int(row[col])                    