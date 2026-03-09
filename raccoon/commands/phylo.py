#!/usr/bin/env python3
"""Phylogenetic QC subcommand wrapper."""
import os
import logging


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
        if assembly_refs and not io.validate_alignment_file(assembly_refs):
            return 1

        treefile = io.resolve_existing_file(getattr(args, 'phylogeny', None), outdir, "Phylogeny file")
        if not treefile:
            return 1
        
        mask_file = getattr(args, 'mask_file', None)
        if mask_file and not io.validate_input_file(mask_file, "Mask file"):
            return 1

        alignment = getattr(args, 'alignment', None)
        if alignment and not io.validate_alignment_file(alignment):
            return 1

        state_file = io.resolve_asr_state_file(getattr(args, 'asr_state', None), treefile)
        if getattr(args, 'asr_state', None) and not state_file:
            return 1
        
        config = {}
        config[KEY_OUTDIR] = outdir
        phylogeny_base = os.path.splitext(os.path.basename(treefile))[0]
        config[KEY_OUTFILENAME] = phylogeny_base
        config[KEY_PHYLOGENY] = phylogeny_base
        config[KEY_RUN_APOBEC3_PHYLO] = args.run_apobec

        outgroup_ids = []
        if args.outgroup_ids:
            outgroup_ids = [x.strip() for x in args.outgroup_ids.split(',') if x.strip()]

        mask_file = mask_file or os.path.join(outdir, f"{phylogeny_base}.mask.csv")

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
                input_cmd_line=getattr(args, "input_cmd_line", None),
            )
        except Exception:
            logging.exception("Failed to generate phylo report")
        logging.info("Phylogenetic QC finished")
        return 0
    except Exception:
        logging.exception("Phylo QC failed")
        return 2
