#!/usr/bin/env python3
"""Alignment subcommand wrapper."""
import os
import logging

from raccoon.utils import constants as rc

def main(args):
    """Run alignment QC using functions from raccoon.utils.

    Expects an argparse Namespace with attributes added by the top-level parser.
    """
    if not hasattr(args, 'alignment'):
        raise ValueError("Expected argparse Namespace from raccoon.command.build_parser")

    logging.info("Starting alignment QC")
    try:
        # lazy import to keep module lightweight at import-time
        from raccoon.utils import alignment_functions as af
        from raccoon.utils import io

        outdir = args.outdir or os.getcwd()
        
        # validate output directory
        if not io.ensure_output_directory(outdir):
            return 1
        
        # validate input files
        if not io.validate_alignment_file(args.alignment):
            return 1
        
        genbank = getattr(args, rc.KEY_GENBANK, None)
        if genbank and not io.validate_genbank_file(genbank):
            return 1
        
        reference = getattr(args, rc.KEY_REFERENCE_ID, None)
        max_n_content = getattr(args, rc.KEY_MAX_N_CONTENT, rc.DEFAULT_MAX_N_CONTENT)
        cluster_window = getattr(args, rc.KEY_CLUSTER_WINDOW, rc.DEFAULT_CLUSTER_WINDOW)
        cluster_count = getattr(args, rc.KEY_CLUSTER_COUNT, rc.DEFAULT_CLUSTER_COUNT)
        mask_clustered = getattr(args, rc.KEY_MASK_CLUSTERED, True) 
        mask_n_adjacent = getattr(args, rc.KEY_MASK_N_ADJACENT, True) 
        mask_gap_adjacent = getattr(args, rc.KEY_MASK_GAP_ADJACENT, True) 
        mask_frame_break = getattr(args, rc.KEY_MASK_FRAME_BREAK, True)

        summary = af.run_alignment_qc(
            args.alignment,
            outdir=outdir,
            genbank_path=genbank,
            reference_id=reference,
            max_n_content=max_n_content,
            cluster_window=cluster_window,
            cluster_count=cluster_count,
            mask_clustered=mask_clustered,
            mask_n_adjacent=mask_n_adjacent,
            mask_gap_adjacent=mask_gap_adjacent,
            mask_frame_break=mask_frame_break,
        )


        logging.info("Alignment QC completed")
        logging.info(f"High N sequences: {len(summary['high_n_sequences'])}")
        logging.info(f"Sites to mask: {len(summary['sites_to_mask'])}")

        # write a summary file
        summary_file = os.path.join(outdir, 'alignment_qc_summary.txt')
        with open(summary_file, 'w') as fw:
            fw.write(f"issues_found: {summary['issues_found']}\n")
            fw.write(f"high_n_sequences: {len(summary['high_n_sequences'])}\n")
            fw.write(f"sites_to_mask: {len(summary['sites_to_mask'])}\n")
            fw.write(f"mask_file: {summary['mask_file']}\n")

        logging.info(f"Summary written to {summary_file}")

        try:
            from raccoon.utils import reporting
            reporting.generate_alignment_report(
                outdir=outdir,
                alignment_path=args.alignment,
                mask_file=summary.get(rc.KEY_MASK_FILE),
            )
        except Exception:
            logging.exception("Failed to generate alignment report")

        return 0
    except FileNotFoundError as exc:
        logging.error(str(exc))
        return 1
    except Exception:
        logging.exception("Alignment QC failed")
        return 2
