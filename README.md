# raccoon

<p align="center">
  <img src="docs/raccoon_logo.png" alt="raccoon logo" width="240" />
</p>

<p align="center"><strong>Rigorous Alignment Curation: Cleanup Of Outliers and Noise</strong></p>

Raccoon is a lightweight toolkit for alignment and phylogenetic QC workflows. It identifies problematic sites (e.g., clustered SNPs, SNPs near Ns/gaps, and frame‑breaking indels) and produces mask files and summaries for downstream analyses.

---

## Contents

- [Use cases](#use-cases)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [CLI usage](#cli-usage)
- [Mask notes](#mask-notes)
- [Example data](#example-data)

## Use cases

- Flag clustered SNPs that may indicate contamination, recombination, or misalignment.
- Detect SNPs adjacent to low-coverage regions (Ns) or gaps.
- Identify frame-breaking indels in coding regions using a GenBank reference.
- Generate mask files to exclude suspect sites prior to phylogenetic or evolutionary analyses.

## Installation

From source:

```bash
pip install .
```

For development (editable install):

```bash
pip install -e .
```

## Quickstart

## CLI usage
Show help:

```bash
raccoon --help
```

### Sequence QC (`seq-qc`)

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

With a custom header template:

```bash
raccoon seq-qc -f a.fasta b.fasta -o combined.fasta \
  -m metadata.csv --header-fields "{id}|{country}|{date}"
```

Key options:

- `-m, --metadata`: metadata CSV file(s) for header harmonisation
- `--metadata-delimiter`: metadata delimiter (default `,`; `.tsv` auto-detected)
- `--metadata-id-field`: metadata ID column (default: `sample`)
- `--metadata-location-field`: metadata location column (default: `location`)
- `--metadata-date-field`: metadata date column (default: `date`)
- `--header-fields`: template for custom headers (e.g. `{id}|{country}|{date}`)
- `--header-separator`: separator used for non-template harmonised headers (default: `|`)
- `--seq-id-delimiter`: delimiter for parsing IDs from input headers (default: `|`)
- `--seq-id-field-index`: 0-based field index for parsed sequence ID (default: `0`)
- `--min-length`: minimum sequence length to keep
- `--max-n-content`: maximum N-content proportion to keep

### Alignment QC (`aln-qc`)

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

- `--max-n-content`: N-content threshold for flagging
- `--cluster-window`: window size (bp) for clustered SNP detection
- `--cluster-count`: minimum SNPs in-window to mark as clustered
- `--no-flag-clustered`: skip clustered SNP flagging
- `--no-flag-n-adjacent`: skip N-adjacent SNP flagging
- `--no-flag-gap-adjacent`: skip gap-adjacent SNP flagging
- `--no-flag-frame-break`: skip frame-breaking indel flagging
- `--flag-removal-threshold`: mark sequence for removal above this flagged-site count

### Apply mask (`mask`)

```bash
raccoon mask <alignment.fasta> \
  --mask-file results/alignment_qc/mask_sites.csv \
  -d results/alignment_qc
```

Key options:

- `--mask-file`: mask CSV file from `aln-qc`
- `--mask-character`: character to use for masking (default: `?`)
- `-o, --outfile`: output masked alignment file name
- `-d, --outdir`: output directory
- `-t, --sequence-type`: `nt` or `aa` (default: `nt`)

### Phylogenetic QC (`tree-qc`)

Basic usage:

```bash
raccoon tree-qc --tree <treefile> -d outdir \
  --alignment <alignment.fasta> --asr-state <treefile>.state \
  --run-adar --adar-window 300 --adar-min-count 3
```

Key options:

- `-t, --tree`: input phylogeny file (required)
- `--tree-format`: `auto`, `newick`, or `nexus`
- `--assembly-refs`: assembly/reference FASTA used for mapping
- `--outgroup-ids`: comma-separated outgroup sequence IDs
- `--mask-file`: optional mask CSV with sites to ignore
- `--tip-fields`: template for parsing tip-label fields
- `--tip-field-delimiter`: delimiter used for tip field parsing
- `--tip-date-field`: field name treated as date in tip parsing
- `--midpoint-root`: midpoint-root tree for report visualisation (ignored with `--asr-state`)
- `--long-branch-sd`: SD threshold for long-branch flagging
- `--run-apobec`: run APOBEC3 checks
- `--run-adar`: run ADAR checks
```

See full CLI details in [docs/cli.md](docs/cli.md).

## Mask notes

Mask output uses the following note values:

| Note | Meaning |
| --- | --- |
| clustered_snps | Clustered SNPs within the configured window. |
| N_adjacent | SNPs adjacent to an N run within the configured window. |
| gap_adjacent | SNPs adjacent to a gap within the configured window. |
| frame_break | Gap sites that break the CDS frame length. |

## Example data

The [examples](examples) folder includes a constructed alignment and GenBank reference suitable for quick testing:

- [examples/constructed_alignment.fasta](examples/constructed_alignment.fasta)
- [examples/constructed_reference.gb](examples/constructed_reference.gb)