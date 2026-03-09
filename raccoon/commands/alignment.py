#!/usr/bin/env python3
"""Alignment subcommand wrapper."""
import os
import logging

from raccoon.utils import constants as rc
from raccoon.utils.misc import green

def _build_flagged_criteria(
    *,
    no_flag_clustered: bool,
    no_flag_n_adjacent: bool,
    no_flag_gap_adjacent: bool,
    no_flag_frame_break: bool,
    cluster_count: int,
    cluster_window: int,
) -> str:
    parts = []
    if not no_flag_clustered:
        parts.append(f"clustered SNPs (≥{cluster_count} SNPs within {cluster_window} bp)")
    if not no_flag_n_adjacent:
        parts.append("SNPs adjacent to Ns")
    if not no_flag_gap_adjacent:
        parts.append("SNPs adjacent to gaps")
    if not no_flag_frame_break:
        parts.append("frame-breaking indels")
    if not parts:
        return "none (all flagging checks disabled)"
    return "; ".join(parts)


def _build_removal_criteria(flag_removal_threshold: int) -> str:
    """Build a human-readable description of sequence removal criteria."""
    return f"Sequences with more than {flag_removal_threshold} flagged site(s)"

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
        no_flag_clustered = getattr(args, rc.KEY_NO_FLAG_CLUSTERED, False)
        no_flag_n_adjacent = getattr(args, rc.KEY_NO_FLAG_N_ADJACENT, False)
        no_flag_gap_adjacent = getattr(args, rc.KEY_NO_FLAG_GAP_ADJACENT, False)
        no_flag_frame_break = getattr(args, rc.KEY_NO_FLAG_FRAME_BREAK, False)
        flag_removal_threshold = getattr(args, rc.KEY_FLAG_REMOVAL_THRESHOLD, rc.DEFAULT_FLAG_REMOVAL_THRESHOLD)

        if no_flag_clustered:
            print(green("NOTE: Clustered SNP flagging is disabled"))
        if no_flag_n_adjacent:
            print(green("NOTE: SNPs near Ns will not be flagged"))
        if no_flag_gap_adjacent:
            print(green("NOTE: SNPs near gaps will not be flagged"))

        summary = af.run_alignment_qc(
            args.alignment,
            outdir=outdir,
            genbank_path=genbank,
            reference_id=reference,
            max_n_content=max_n_content,
            cluster_window=cluster_window,
            cluster_count=cluster_count,
            no_flag_clustered=no_flag_clustered,
            no_flag_n_adjacent=no_flag_n_adjacent,
            no_flag_gap_adjacent=no_flag_gap_adjacent,
            no_flag_frame_break=no_flag_frame_break,
            flag_removal_threshold=flag_removal_threshold
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
            flagged_criteria = _build_flagged_criteria(
                no_flag_clustered=no_flag_clustered,
                no_flag_n_adjacent=no_flag_n_adjacent,
                no_flag_gap_adjacent=no_flag_gap_adjacent,
                no_flag_frame_break=no_flag_frame_break,
                cluster_count=cluster_count,
                cluster_window=cluster_window,
            )
            flagged_removal_criteria = _build_removal_criteria(flag_removal_threshold)
            reporting.generate_alignment_report(
                outdir=outdir,
                alignment_path=args.alignment,
                mask_file=summary.get(rc.KEY_MASK_FILE),
                flagged_criteria=flagged_criteria,
                flagged_removal_criteria=flagged_removal_criteria,
                input_cmd_line=getattr(args, "input_cmd_line", None),
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
