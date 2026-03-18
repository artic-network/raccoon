# raccoon CLI

Use the `raccoon` top-level command with subcommands for different QC tasks.

Top-level usage:

```bash
raccoon <subcommand> [options]
```

Run `raccoon <subcommand> --help` for full command-specific help.

## `seq-qc` subcommand

Purpose: combine one or more FASTA files into a single upper-case, unwrapped FASTA, with optional metadata-driven header harmonisation.

Basic usage:

```bash
raccoon seq-qc -f a.fasta b.fasta -o combined.fasta
```

With metadata-driven headers:

```bash
raccoon seq-qc -f a.fasta b.fasta -o combined.fasta \
  -m metadata.csv other_metadata.csv \
  --metadata-id-field sample \
  --metadata-location-field location \
  --metadata-date-field date \
  --header-separator '|'
```

With custom template fields:

```bash
raccoon seq-qc -f a.fasta b.fasta -o combined.fasta \
  -m metadata.csv \
  --header-fields "{id}|{country}|{date}"
```

Key options:

- `-f, --fasta`: input FASTA files (one or more)
- `-o, --outfile`: output FASTA file (use `-` for stdout)
- `-m, --metadata`: one or more metadata CSV files
- `--metadata-delimiter`: metadata delimiter (default `,`; `.tsv` auto-detected)
- `--metadata-id-field`: metadata ID column (default: `sample`)
- `--metadata-location-field`: metadata location column (default: `location`)
- `--metadata-date-field`: metadata date column (default: `date`)
- `--header-fields`: template for custom headers (takes precedence over location/date field args)
- `--header-separator`: separator used for non-template harmonised headers (default: `|`)
- `--seq-id-delimiter`: delimiter for parsing IDs from input headers (default: `|`)
- `--seq-id-field-index`: 0-based field index used for parsed sequence ID (default: `0`)
- `--min-length`: minimum sequence length to keep
- `--max-n-content`: maximum N-content proportion to keep

Notes:

- `--header-fields` supports metadata column names in braces, e.g. `{id}|{region}|{date}`.
- If `--header-fields` is provided, location/date-specific metadata field args are ignored.

## `aln-qc` subcommand

Purpose: run alignment quality-control checks and produce a mask file and summary.

Basic usage:

```bash
raccoon aln-qc <alignment.fasta> -d outdir
```

With GenBank reference for frame-break checks:

```bash
raccoon aln-qc <alignment.fasta> -d outdir \
  --genbank <reference.gb> --reference-id <ref_id>
```

Disable selected flag classes:

```bash
raccoon aln-qc <alignment.fasta> -d outdir \
  --no-flag-n-adjacent --no-flag-gap-adjacent
```

Key options:

- `alignment` (positional): input alignment (FASTA)
- `-d, --outdir`: output directory
- `-t, --sequence-type`: `nt` or `aa` (default: `nt`)
- `--genbank`: GenBank file for frame-breaking indel checks
- `--reference-id`: reference sequence ID used with GenBank features
- `--max-n-content`: N-content threshold for flagging
- `--cluster-window`: clustered SNP detection window size (bp)
- `--cluster-count`: minimum SNP count in window for clustered flagging
- `--no-flag-clustered`: skip clustered SNP flagging
- `--no-flag-n-adjacent`: skip N-adjacent SNP flagging
- `--no-flag-gap-adjacent`: skip gap-adjacent SNP flagging
- `--no-flag-frame-break`: skip frame-breaking indel flagging
- `--flag-removal-threshold`: mark sequence for removal above this flagged-site count

## `mask` subcommand

Purpose: apply an `aln-qc` mask CSV to an alignment and write a masked FASTA.

Basic usage:

```bash
raccoon mask data/alignment.fasta \
  --mask-file results/alignment_qc/mask_sites.csv \
  -d results/alignment_qc
```

Key options:

- `alignment` (positional): input alignment (FASTA)
- `--mask-file`: mask CSV from `aln-qc`
- `--mask-character`: character to use for masking (default: `?`)
- `-o, --outfile`: output masked alignment file name
- `-d, --outdir`: output directory (default: `.`)
- `-t, --sequence-type`: `nt` or `aa` (default: `nt`)

## `tree-qc` subcommand

Purpose: run phylogenetic QC and generate an interactive tree report.

Basic usage:

```bash
raccoon tree-qc --tree <treefile> -d outdir \
  --alignment <alignment.fasta> --asr-state <treefile>.state \
  --run-adar --adar-window 300 --adar-min-count 3
```

Key options:

- `-t, --tree`: input tree file or basename (required)
- `--tree-format`: `auto`, `newick`, or `nexus` (default: `auto`)
- `-d, --outdir`: output directory
- `--outgroup-ids`: comma-separated outgroup IDs
- `--alignment`: alignment FASTA used with ASR state file
- `--asr-state`: ancestral state reconstruction file
- `--mask-file`: optional mask CSV file with sites to ignore
- `--assembly-refs`: assembly/reference FASTA used for mapping
- `--long-branch-sd`: SD threshold for long-branch flagging (default: `3.0`)
- `--run-apobec`: run APOBEC3 phylo checks
- `--run-adar`: run ADAR phylo checks
- `--adar-window`: max distance (bp) for ADAR cluster window (default: `300`)
- `--adar-min-count`: min ADAR sites in window to flag branch (default: `3`)
- `--tip-fields`: template for parsing tip label fields
- `--tip-field-delimiter`: delimiter for tip field parsing
- `--tip-date-field`: field name treated as date in tip parsing
- `--midpoint-root`: midpoint-root for report visualisation (applied only when `--asr-state` is not provided)
- `--adar-window`: max distance (bp) for ADAR cluster window (default: `300`)
- `--adar-min-count`: min ADAR sites in window to flag branch (default: `3`)
- `--height`: optional figure height

## Global options

- `-v, --version`: show version and exit
- `-V, --verbose`: increase logging verbosity (repeat for more detail)
