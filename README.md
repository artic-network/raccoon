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

```bash
raccoon aln-qc examples/constructed_alignment.fasta -d outdir \
	--genbank examples/constructed_reference.gb --reference-id ref
```

Outputs:

- mask_sites.csv
- alignment_qc_summary.txt

## CLI usage

Show help:

```bash
raccoon --help
```

Alignment QC:

```bash
raccoon aln-qc <alignment.fasta> -d outdir
```

With a GenBank reference for frame‑break detection:

```bash
raccoon aln-qc <alignment.fasta> -d outdir \
  --genbank <reference.gb> --reference-id <ref_id>
```

Masking toggles (defaults are enabled):

```bash
raccoon aln-qc <alignment.fasta> -d outdir \
  --no-mask-n-adjacent --no-mask-gap-adjacent
```

Key alignment options:

- `--n-threshold`: fraction of Ns allowed per sequence before flagging.
- `--cluster-window`: window size (bp) for clustered SNP detection.
- `--cluster-count`: minimum SNPs within a window to flag as clustered.
- `--mask-clustered/--no-mask-clustered`: include/exclude clustered SNPs.
- `--mask-n-adjacent/--no-mask-n-adjacent`: include/exclude SNPs adjacent to Ns.
- `--mask-gap-adjacent/--no-mask-gap-adjacent`: include/exclude SNPs adjacent to gaps.
- `--mask-frame-break/--no-mask-frame-break`: include/exclude frame-breaking indels.

Sequence QC:

```bash
raccoon seq-qc a.fasta b.fasta -o combined.fasta
```

With metadata-driven headers:

```bash
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv other_metadata.csv --metadata-id-field id \
  --metadata-location-field location --metadata-date-field date \
  --header-separator '|'
```

With custom header template:

```bash
# Custom field order and custom fields
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv --header-fields "{id}|{country}|{date}"

# Multiple location levels
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv --header-fields "{id}|{region}|{country}|{date}"

# Custom separator
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv --header-fields "{id}_{location}_{date}"
```

**Header field details:**

- `--header-fields` uses a template format: `{field1}|{field2}|{field3}` where field names come from your metadata CSV columns
- Field values are automatically **sanitized** for safe use in phylogenetic tools (Newick parsers):
  - Converted to **lowercase**
  - Spaces, commas, colons, semi-colons, and parentheses are replaced with underscores
  - Hyphens and pipes are preserved (safe for Newick parsers)
  - ISO dates (YYYY-MM-DD) are not further sanitized
  - Multiple consecutive underscores are collapsed to single underscores
  - Example: `"New York, USA"` → `"new_york_usa"`, `"January 15, 2024"` → `"2024-01-15"`
- When `--header-fields` is provided, it **takes precedence** over `--metadata-location-field` and `--metadata-date-field` arguments
- Missing fields in metadata are output as empty values (e.g., `seq1||2024-01-01` if location is missing), allowing the pipeline to continue with issues logged to `seq_qc_metadata_issues.csv`
- The special field `{id}` always maps to the parsed sequence ID (after `--id-field` extraction)
- **CSV parsing note**: Metadata fields containing delimiters (commas) must be properly quoted. If unescaped delimiters are detected, raccoon will exit with an informative error suggesting field quoting (e.g., `"New York, USA"` instead of `New York, USA`)


**Backward compatibility:**
- Default behavior (no `--header-fields`) still works: `{id}|{location}|{date}` format
- Existing `--metadata-location-field` and `--metadata-date-field` args are respected when `--header-fields` is not provided

Key sequence QC options:

- `--metadata`: metadata CSV file(s) for header harmonization
- `--metadata-id-field`: CSV column to match with sequence IDs (default: id)
- `--header-fields`: template for custom header format (e.g., `{id}|{country}|{date}`)
- `--metadata-location-field`: CSV column for location (default: location; overridden by `--header-fields`)
- `--metadata-date-field`: CSV column for date (default: date; overridden by `--header-fields`)
- `--header-separator`: separator between header fields (default: |; ignored if `--header-fields` is used)
- `--id-delimiter`: delimiter for parsing IDs from input headers (default: |)
- `--id-field`: 0-based field index for ID extraction (default: 0)
- `--min-length`: minimum sequence length to keep
- `--max-n-content`: maximum N content proportion to keep (e.g., 0.1 for 10%)


Phylogenetic QC:

```bash
raccoon tree-qc --phylogeny <treefile> -d outdir \
  --alignment <alignment.fasta> --asr-state <treefile>.state \
  --run-adar --adar-window 300 --adar-min-count 3
```

Key phylo options:

- `--phylogeny`: tree file (Newick or Nexus)
- `--alignment`: alignment used for ASR state parsing
- `--asr-state`: ASR state file (defaults to `<treefile>.state` if present)
- `--tree-format`: auto/newick/nexus
- `--run-adar`: enable ADAR-like edit flagging
- `--run-apobec`: enable APOBEC3-like edit flagging
- `--adar-window`: max distance (bp) for ADAR clustering (default: 300)
- `--adar-min-count`: min ADAR sites in window to flag a branch (default: 3)
- `--long-branch-sd`: std dev threshold for long-branch flagging (default: 3.0)
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