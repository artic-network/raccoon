#!/usr/bin/env python3
"""Command-line interface for raccoon toolkit."""

import sys
import os
import argparse
import logging
from raccoon import __version__, _program
from raccoon.commands import alignment as alignment_cmd, phylo as phylo_cmd, combine as combine_cmd, mask as mask_cmd
from raccoon.utils import constants as rc

def build_parser():
    parser = argparse.ArgumentParser(prog=_program, description="raccoon toolkit")
    parser.add_argument('-v','--version', action='version', version=f'raccoon {__version__}', help='Show version and exit')
    parser.add_argument('-V','--verbose', action='count', default=0, help='Increase verbosity (can specify multiple times)')
    sub = parser.add_subparsers(dest='command')
    sub.required = True

    c = sub.add_parser('seq-qc', help='sequence QC (combine and harmonise FASTA files)')
    c.add_argument('-f','--fasta', nargs='+', help='Input FASTA files (one or more)')
    c.add_argument('-m','--metadata', nargs='+', action='extend', default=None, help='Input metadata CSV(s) for FASTA header population (one or more; optional)')
    c.add_argument('--metadata-delimiter', default=rc.DEFAULT_METADATA_DELIMITER, help=f'Metadata delimiter (default: {rc.DEFAULT_METADATA_DELIMITER}; auto-detects .csv and .tsv files)')
    c.add_argument('--metadata-id-field', default=rc.DEFAULT_ID_FIELD, help=f'Metadata id column used for matching sequence ID (default: {rc.DEFAULT_ID_FIELD})')
    c.add_argument('--metadata-location-field', default=rc.DEFAULT_LOCATION_FIELD, help=f'Metadata location column (default: {rc.DEFAULT_LOCATION_FIELD})')
    c.add_argument('--metadata-date-field', default=rc.DEFAULT_DATE_FIELD, help=f'Metadata date column (default: {rc.DEFAULT_DATE_FIELD})')
    c.add_argument('--seq-id-delimiter', default=rc.DEFFAULT_ID_DELIMITER, help=f'Delimiter for parsing IDs from input headers (default: {rc.DEFFAULT_ID_DELIMITER})')
    c.add_argument('--seq-id-field-index', type=int, default=rc.DEFAULT_ID_FIELD_INDEX, help=f'0-based field index for parsing IDs from input headers (default: {rc.DEFAULT_ID_FIELD_INDEX})')
    c.add_argument('--min-length', type=int, default=None, help='Minimum sequence length to keep')
    c.add_argument('--max-n-content', type=float, default=None, help='Maximum N content proportion to keep (e.g. 0.1)')
    c.add_argument('--header-fields', default=rc.DEFAULT_HEADER_FIELDS, help=f'Template for constructing harmonised sequence header (e.g. "{rc.DEFAULT_HEADER_FIELDS}")')
    c.add_argument('--header-separator', default=rc.DEFAULT_HEADER_SEPARATOR, help=f'Header separator (default: {rc.DEFAULT_HEADER_SEPARATOR})')
    c.add_argument('-o', '--output', default='combined.fasta', help='Output FASTA file (use - for stdout)')
    c.set_defaults(func=combine_cmd.main)

    a = sub.add_parser('aln-qc', help='alignment QC')
    a.add_argument('alignment', help='Input alignment fasta file')
    a.add_argument('-t','--sequence-type', choices=['nt','aa'], default='nt', dest='sequence_type', help='Sequence type (default: nt)')
    a.add_argument('-d','--output-dir', default='.', dest='output_dir', help='Output directory')
    a.add_argument('--genbank', help='GenBank file for frame-breaking indel checks', dest='genbank', default=None)
    a.add_argument('--reference-id', help='Reference sequence ID in alignment (for GenBank features)', dest='reference_id', default=None)
    a.add_argument('--max-n-content', type=float, default=0.2, help='N content threshold for flagging (default: 0.2)')
    a.add_argument('--cluster-window', type=int, default=10, dest='cluster_window', help='Window size for clustered SNP detection (default: 10bp)')
    a.add_argument('--cluster-count', type=int, default=3, dest='cluster_count', help='Min SNPs in window to flag as clustered (default: 3)')
    a.add_argument('--mask-clustered', default=True, action=argparse.BooleanOptionalAction, dest='mask_clustered', help='Mask clustered SNPs (default: True)')
    a.add_argument('--mask-n-adjacent', default=True, action=argparse.BooleanOptionalAction, dest='mask_n_adjacent', help='Mask SNPs adjacent to Ns (default: True)')
    a.add_argument('--mask-gap-adjacent', default=True, action=argparse.BooleanOptionalAction, dest='mask_gap_adjacent', help='Mask SNPs adjacent to gaps (default: True)')
    a.add_argument('--mask-frame-break', default=True, action=argparse.BooleanOptionalAction, dest='mask_frame_break', help='Mask frame-breaking indels (default: True)')
    a.set_defaults(func=alignment_cmd.main)

    m = sub.add_parser('mask', help='apply a mask file to an alignment')
    m.add_argument('alignment', help='Input alignment fasta file')
    m.add_argument('--mask-file', required=True, dest='mask_file', help='Mask CSV file from aln-qc')
    m.add_argument('-o', '--output', dest='output', default=None, help='Output masked alignment file')
    m.add_argument('-d', '--outdir', dest='outdir', default='.', help='Output directory')
    m.add_argument('-t', '--sequence-type', choices=['nt','aa'], default='nt', dest='sequence_type', help='Sequence type (default: nt)')
    m.set_defaults(func=mask_cmd.main)

    p = sub.add_parser('tree-qc', help='phylogenetic QC')
    p.add_argument('--assembly-refs', help='Assembly references fasta', dest='assembly_refs', default=None)
    p.add_argument('--phylogeny', help='Phylogeny tree file or base name', dest='phylogeny', required=True)
    p.add_argument('--outdir','-d', help='Output directory', dest='outdir', default='.')
    p.add_argument('--mask-file', help='Mask output file', dest='mask_file', default=None)
    p.add_argument('--outgroup-ids', help='Comma-separated outgroup ids', dest='outgroup_ids', default=None)
    p.add_argument('--alignment', help='Alignment fasta used with ASR state file', dest='alignment', default=None)
    p.add_argument('--asr-state', help='Ancestral state reconstruction file in the format output by IQTREE', dest='asr_state', default=None)
    p.add_argument('--tree-format', choices=['auto','newick','nexus'], default='auto', dest='tree_format', help='Tree format (default: auto)')
    p.add_argument('--long-branch-sd', type=float, default=3.0, dest='long_branch_sd', help='Std dev threshold for long branch flagging (default: 3.0)')
    p.add_argument('--adar-window', type=int, default=300, dest='adar_window', help='Max distance (bp) for ADAR cluster window (default: 300)')
    p.add_argument('--adar-min-count', type=int, default=3, dest='adar_min_count', help='Min ADAR sites in window to flag a branch (default: 3)')
    p.add_argument('--run-apobec', action='store_true', dest='run_apobec', help='Run APOBEC3 phylo checks')
    p.add_argument('--run-adar', action='store_true', dest='run_adar', help='Run ADAR phylo checks')
    p.add_argument('--height', help='Figure height', dest='height', default=None)
    p.set_defaults(func=phylo_cmd.main)

    return parser


def main(sysargs=None):
    if sysargs is None:
        sysargs = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(sysargs)

    # configure logging
    level = logging.WARNING
    if getattr(args,'verbose',0) >=2:
        level = logging.DEBUG
    elif getattr(args,'verbose',0) == 1:
        level = logging.INFO
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')

    try:
        return args.func(args)
    except Exception as e:
        logging.exception("Error running raccoon")
        return 2


if __name__ == '__main__':
    sys.exit(main())