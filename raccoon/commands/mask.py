#!/usr/bin/env python3
"""Mask subcommand wrapper."""
import os
import logging


def main(args):
    if not hasattr(args, "alignment"):
        raise ValueError("Expected argparse Namespace from raccoon.command.build_parser")

    logging.info("Starting mask")
    try:
        from raccoon.utils import alignment_functions as af
        from raccoon.utils import io

        outdir = args.outdir or os.getcwd()
        if not io.ensure_output_directory(outdir):
            return 1

        alignment_path, temp_path = io.prepare_alignment_input(args.alignment)
        try:
            if not io.validate_alignment_file(alignment_path):
                return 1

            mask_file = getattr(args, "mask_file", None)
            if not io.validate_input_file(mask_file, "Mask file"):
                return 1

            output_path = args.outfile
            if not output_path:
                base = os.path.splitext(os.path.basename(alignment_path))[0]
                output_path = f"{base}.masked.fasta"

            if not os.path.isabs(output_path):
                output_path = os.path.join(outdir, output_path)

            mask_char = "X" if getattr(args, "sequence_type", "nt") == "aa" else "N"
            masked_sites = af.apply_mask_to_alignment(
                alignment_path,
                mask_file,
                output_path,
                mask_char=mask_char,
            )
            logging.info(f"Masked alignment written to {output_path} ({masked_sites} sites masked)")
            try:
                from raccoon.utils import reporting
                reporting.generate_mask_report(
                    outdir=outdir,
                    alignment_path=alignment_path,
                    mask_file=mask_file,
                    output_alignment=output_path,
                )
            except Exception:
                logging.exception("Failed to generate mask report")
            logging.info("Mask finished")
            return 0
        finally:
            io.cleanup_temp_file(temp_path)
    except FileNotFoundError as exc:
        logging.error(str(exc))
        return 1
    except Exception:
        logging.exception("Mask failed")
        return 2
