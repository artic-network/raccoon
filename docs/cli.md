# raccoon CLI

Use the `raccoon` top-level command with subcommands for different QC tasks.

Examples:

- Alignment QC
  ```bash
  raccoon aln-qc input_alignment.fasta -d outdir --flag-n
  ```

- Phylogenetic QC
  ```bash
  raccoon tree-qc --phylogeny mytree --assembly-refs refs.fasta -d outdir --run-apobec
  ```

- Sequence QC
  ```bash
  raccoon seq-qc a.fasta b.fasta -o combined.fasta
  ```

Each subcommand has its own help available with `raccoon <subcommand> --help`.

# raccoon CLI

Use the `raccoon` top-level command with subcommands for different QC tasks.

Top-level usage

```bash
raccoon <subcommand> [options]
```

Run `raccoon <subcommand> --help` to see subcommand-specific options.

aln-qc subcommand

Purpose: run alignment quality-control checks and produce a mask file and summary.

Basic usage:

```bash
raccoon aln-qc <alignment.fasta> -d outdir
```

Key options

- `alignment` (positional): path to the input alignment (FASTA)
- `-d, --outdir`: directory to write outputs (created if missing)
- `--genbank`: optional GenBank file containing CDS/features for frame-breaking indel checks
- `--reference-id`: optional sequence id in the GenBank file used as reference mapping
- `--n-threshold`: fraction of N allowed per sequence before flagged (default: 0.2)
- `--cluster-window`: window size (bp) to search for clustered SNPs (default: 10)
- `--cluster-count`: minimum number of SNPs within the window to be considered clustered (default: 3)
- `--mask-clustered/--no-mask-clustered`: include/exclude clustered SNPs in mask output (default: include)
- `--mask-n-adjacent/--no-mask-n-adjacent`: include/exclude SNPs adjacent to Ns in mask output (default: include)
- `--mask-gap-adjacent/--no-mask-gap-adjacent`: include/exclude SNPs adjacent to gaps in mask output (default: include)
- `--mask-frame-break/--no-mask-frame-break`: include/exclude frame-breaking indels in mask output (default: include)

Mask output

- The mask CSV now contains `flagged`, `type`, `minimum`, `maximum`, `length`, `present_in`, `note`.
- `type` is either `site` (mask a single column) or `sequence_record` (remove sequence).
- If a sequence has more than 20 flagged sites, it is emitted as a `sequence_record` entry and individual site rows for that sequence are omitted; the `note` column lists contributing sites.

Behavior and exit codes

- The command will try to create `--outdir` if it does not exist and will verify it is writable.
- Optional input files (`--genbank`, `--reference-id`) are validated if provided.
- Exit code `0` on success, `1` on validation/file errors, `2` on unexpected failures.

Example

```bash
raccoon aln-qc data/sequences.fasta -d results/alignment_qc --genbank refs/ref.gb --reference-id NC_000000
```

Example with masking toggles

```bash
raccoon aln-qc data/sequences.fasta -d results/alignment_qc \
  --genbank refs/ref.gb --reference-id NC_000000 \
  --no-mask-n-adjacent --no-mask-gap-adjacent
```

mask subcommand

Purpose: apply an aln-qc mask CSV to an alignment and write a masked FASTA.

Basic usage:

```bash
raccoon mask data/alignment.fasta --mask-file results/alignment_qc/mask_sites.csv -d results/alignment_qc
```

Key options

- `alignment` (positional): path to the input alignment (FASTA)
- `--mask-file`: mask CSV from aln-qc
- `-o, --outfile`: output masked alignment file name (default: <alignment>.masked.fasta)
- `-d, --outdir`: output directory (default: .)
- `-t, --sequence-type`: sequence type (nt/aa); uses N or X (default: nt)

tree-qc subcommand

Purpose: run phylogenetic QC (SNP anomaly checks, apobec analyses, plotting helpers).

Basic usage:

```bash
raccoon tree-qc --phylogeny treefile.newick --assembly-refs refs.fasta -d outdir
```

Key options

- `--phylogeny`: path to Newick tree file
- `--assembly-refs`: path to assembly/reference FASTA used for mapping
- `-d, --outdir`: directory to write outputs (created if missing)
- `--outgroup-ids`: comma-separated list of outgroup sequence ids
- `--mask-file`: optional mask CSV file with sites to ignore
- `--height`: plotting height parameter (optional)
- `--run-apobec`: flag to run APOBEC3-specific analyses
- `--run-adar`: flag to run ADAR-specific analyses
- `--alignment`: alignment FASTA used with ASR state file
- `--asr-state`: ancestral state reconstruction file (defaults to <tree>.state if present)
- `--tree-format`: tree format (auto/newick/nexus)
- `--long-branch-sd`: std dev threshold for long-branch flagging

Behavior and exit codes

- The command validates `--assembly-refs` and optional mask/tree files, and ensures `--outdir` is writable.
- Exit code `0` on success, `1` on validation/file errors, `2` on unexpected failures.

Example

```bash
raccoon tree-qc --phylogeny trees/rep.tree --assembly-refs data/refs.fasta -d results/phylo --run-apobec
```

seq-qc subcommand

Purpose: combine one or more FASTA files into a single upper-case, unwrapped FASTA, with optional metadata-driven header harmonisation.

Basic usage:

```bash
raccoon seq-qc a.fasta b.fasta -o combined.fasta
```

Header harmonisation using metadata (backward compatible format):

```bash
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv \
  --metadata-id-field id \
  --metadata-location-field location \
  --metadata-date-field date \
  --header-separator '|'
```

Header harmonisation with custom template:

```bash
# Custom field order and custom fields from metadata
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv \
  --header-fields "{id}|{country}|{date}"

# Multiple location levels
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv \
  --header-fields "{id}|{region}|{country}|{date}"

# Different separator
raccoon seq-qc a.fasta b.fasta -o combined.fasta \
  --metadata metadata.csv \
  --header-fields "{id}:{location}:{date}"
```

Key options

- `inputs` (positional): one or more input FASTA files
- `-o, --outfile`: output FASTA file (use `-` for stdout)
- `--metadata`: one or more metadata CSVs used to harmonise headers
- `--metadata-delimiter`: metadata delimiter (default: `,`; auto-detects `.tsv` files)
- `--metadata-id-field`: metadata id column (default: `id`)
- `--header-fields`: template for custom header format (e.g. `{id}|{country}|{date}`); **takes precedence** over `--metadata-location-field` and `--metadata-date-field` if provided
- `--metadata-location-field`: metadata location column (default: `location`; ignored if `--header-fields` is used)
- `--metadata-date-field`: metadata date column (default: `date`; ignored if `--header-fields` is used)
- `--header-separator`: deprecated; use `--header-fields` instead for custom separators (default: `|`)
- `--id-delimiter`: delimiter for parsing sequence IDs from input headers (default: `|`)
- `--id-field`: 0-based field index for ID extraction (default: `0`)
- `--min-length`: minimum sequence length to keep
- `--max-n-content`: maximum N content proportion to keep (e.g. `0.1`)

Header field format and sanitization

- `--header-fields` uses template syntax: `{field1}{sep}{field2}{sep}{field3}` where field names match your metadata CSV columns
- Field values are **automatically sanitized**: converted to lowercase and spaces/special characters replaced with underscores
  - Example: `"United Kingdom"` → `"united_kingdom"`, `"2024-01-01"` → `"2024-01-01"`
- Missing metadata columns are output as empty values (e.g. `seq1||2024-01-01` if location is missing), keeping consistent numbers of fields in each header
- The special field `{id}` always refers to the parsed sequence ID (after `--id-field` extraction), not the metadata id column
- When `--header-fields` is provided, `--metadata-location-field` and `--metadata-date-field` are ignored

Output files

When filtering or using metadata:
- `seq_qc_filter_failures.csv`: sequences filtered by `--min-length` or `--max-n-content`
- `seq_qc_metadata_issues.csv`: missing metadata rows or missing field values
- `seq-qc_report.html`: summary report (if metadata provided)
