#!/usr/bin/env python3
"""Phylogenetic QC subcommand wrapper."""
import os
import logging
from Bio import SeqIO, Phylo
from raccoon.utils import constants as rc

def _parse_tip_field_names(tip_fields, tip_field_delimiter="|"):
    parsed = []
    for field in str(tip_fields or "").split(tip_field_delimiter):
        name = field.strip()
        if not name:
            continue
        if name.startswith("{") and name.endswith("}") and len(name) > 2:
            name = name[1:-1].strip()
        if name:
            parsed.append(name)
    return parsed



def _detect_biophylo_format(treefile: str, tree_format: str = "auto") -> str:
    if tree_format in {"newick", "nexus"}:
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


def _midpoint_root_tree_file(treefile: str, outdir: str, tree_format: str = "auto") -> str:
    fmt = _detect_biophylo_format(treefile, tree_format=tree_format)
    try:
        tree = Phylo.read(treefile, fmt)
        tree.root_at_midpoint()
        base = os.path.splitext(os.path.basename(treefile))[0]
        ext = ".nexus" if fmt == "nexus" else ".nwk"
        rooted_path = os.path.join(outdir, f"{base}.midpoint_rooted{ext}")
        Phylo.write(tree, rooted_path, fmt)
        return rooted_path
    except Exception:
        return treefile


def main(args):
    """Run phylogenetic QC.

    The implementation defers heavy imports until called so tests can import the package without needing all scientific deps.
    """
    logging.info("Starting phylogenetic QC")
    if not hasattr(args, 'phylogeny'):
        raise ValueError("Expected argparse Namespace from raccoon.command.build_parser")

    try:
        # Lazy import
        from raccoon.utils import phylo_functions as pf
        from raccoon.utils import io
        from raccoon.utils.constants import KEY_OUTDIR, KEY_OUTFILENAME, KEY_PHYLOGENY, KEY_RUN_APOBEC3_PHYLO

        outdir = args.outdir or os.getcwd()
        
        # validate output directory
        if not io.ensure_output_directory(outdir):
            return 1
        
        # validate input files
        assembly_refs = getattr(args, 'assembly_refs', None)
        mask_file = getattr(args, 'mask_file', None)
        alignment = getattr(args, 'alignment', None)
        
        for file in [alignment, assembly_refs, mask_file]:
            if file and not io.validate_input_file(file, "Input file"):
                return 1

        if alignment and not io.validate_alignment_file(alignment):
            return 1

        treefile = io.resolve_existing_file(getattr(args, 'phylogeny', None), outdir, "Phylogeny file")
        if not treefile:
            return 1

        state_file = io.resolve_asr_state_file(getattr(args, 'asr_state', None), treefile)
        if getattr(args, 'asr_state', None) and not state_file:
            return 1

        # Determine if midpoint rooting should be applied
        midpoint_root = bool(getattr(args, 'midpoint_root', False))
        midpoint_root_for_report = midpoint_root and not bool(state_file)
        if midpoint_root and state_file:
            logging.info("--midpoint-root ignored because --asr-state was provided")
        
        midpoint_root_for_report = midpoint_root and not args.outgroup_ids
        if midpoint_root and args.outgroup_ids:
            logging.info("--midpoint-root ignored because --outgroup-ids was provided")
        
        # Parse tip label fields
        tip_field_delimiter = getattr(args, 'tip_field_delimiter', rc.DEFAULT_HEADER_SEPARATOR)
        tip_date_field = getattr(args, 'tip_date_field', rc.DEFAULT_DATE_FIELD)
        tip_fields = getattr(args, 'tip_fields', rc.DEFAULT_HEADER_FIELDS)

        outgroup_ids = []
        if args.outgroup_ids:
            outgroup_ids = [x.strip() for x in args.outgroup_ids.split(',') if x.strip()]

        phylogeny_base = os.path.splitext(os.path.basename(treefile))[0]
        mask_file = mask_file or os.path.join(outdir, f"{phylogeny_base}.mask.csv")

        treefile = _midpoint_root_tree_file(treefile, outdir, tree_format=args.tree_format) if midpoint_root_for_report else treefile

        flags_csv = pf.run_phylo_qc(
            treefile=treefile,
            tree_format=args.tree_format,
            outdir=outdir,
            alignment=alignment,
            state_file=state_file,
            assembly_refs=assembly_refs,
            long_branch_sd=args.long_branch_sd,
            include_apobec=args.run_apobec,
            include_adar=args.run_adar,
            adar_window=args.adar_window,
            adar_min_count=args.adar_min_count,
        )
        try:
            from raccoon.utils import reporting
            reporting.generate_phylo_report(
                outdir=outdir,
                treefile=treefile,
                flags_csv=flags_csv,
                tree_format=args.tree_format,
                tip_fields=tip_fields,
                tip_field_delimiter=tip_field_delimiter,
                tip_date_field=tip_date_field,
                midpoint_root=midpoint_root_for_report,
                outgroup_ids=outgroup_ids if outgroup_ids else None,
                input_cmd_line=getattr(args, "input_cmd_line", None),
            )
        except Exception:
            logging.exception("Failed to generate phylo report")
        logging.info("Phylogenetic QC finished")
        return 0
    except Exception:
        logging.exception("Phylo QC failed")
        return 2
