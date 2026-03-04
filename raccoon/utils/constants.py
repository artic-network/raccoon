# Constants for CSV/field names and notes to avoid typos
KEY_NAME = "Name"
KEY_MINIMUM = "Minimum"
KEY_MAXIMUM = "Maximum"
KEY_LENGTH = "Length"
KEY_PRESENT_IN = "present_in"
KEY_NOTE = "note"
KEY_SITES_TO_MASK = "sites_to_mask"
KEY_MASK_FILE = "mask_file"
KEY_ISSUES_FOUND = "issues_found"
KEY_CLUSTERED_SNP_COUNT = "clustered_snp_count"

# Note values
NOTE_CLUSTERED_SNPS = "clustered_snps"
NOTE_N_ADJACENT = "N_adjacent"
NOTE_GAP_ADJACENT = "gap_adjacent"
NOTE_REVERSION = "reversion"
NOTE_CONVERGENT = "convergent_snp"
NOTE_FRAME_BREAK = "frame_break"
NOTE_UNIQUE_SNP = "unique_snp"

# General config keys (migrated from utils/config.py)
KEY_INPUT_FASTA = "fasta"
KEY_REFERENCE_FASTA = "reference_fasta"
KEY_TO_MASK = "to_mask"
KEY_GENE_BOUNDARIES = "gene_boundaries"
KEY_OUTDIR = "outdir"
KEY_EXCLUDE_FILE = "exclude_file"
KEY_OUTFILE_STEM = "outfile_stem"
KEY_OUTFILENAME = "outfilename"
KEY_VERBOSE = "verbose"
KEY_THREADS = "threads"
KEY_PHYLO_THREADS = "phylo_threads"
KEY_NO_MASK = "no_mask"
KEY_SEQUENCE_MASK = "sequence_mask"
KEY_TRIM_END = "trim_end"
KEY_EXTRACT_CDS = "extract_cds"
KEY_SEQ_QC = "seq_qc"
KEY_ASSEMBLY_REFERENCES = "assembly_references"

KEY_OUTGROUPS = "outgroups"
KEY_PHYLOGENY = "phylogeny"
KEY_PHYLOGENY_SVG = "phylogeny_svg"
KEY_OUTGROUP_STRING = "outgroup_string"
KEY_GRANTHAM_SCORES = "grantham_scores"

KEY_TREE = "tree"
KEY_BRANCH_RECONSTRUCTION = "branch_reconstruction"
KEY_ASR_TREE = "asr_tree"
KEY_ASR_STATE = "asr_state"
KEY_ASR_ALIGNMENT = "asr_alignment"

KEY_FIG_HEIGHT = "fig_height"
KEY_FIG_WIDTH = "fig_width"
KEY_POINT_STYLE = "point_style"
KEY_POINT_JUSTIFY = "point_justify"

# Additional feature flags used in code
KEY_RUN_APOBEC3_PHYLO = "run_apobec3_phylo"

KEY_GENBANK = "genbank"
KEY_REFERENCE_ID = "reference_id"
KEY_N_THRESHOLD = "n_threshold"
KEY_CLUSTER_WINDOW = "cluster_window"
KEY_CLUSTER_COUNT = "cluster_count"
KEY_MASK_CLUSTERED = "mask_clustered"
KEY_MASK_N_ADJACENT = "mask_n_adjacent"
KEY_MASK_GAP_ADJACENT = "mask_gap_adjacent"
KEY_MASK_FRAME_BREAK = "mask_frame_break"

DEFAULT_ID_FIELD = "sample"
DEFAULT_LOCATION_FIELD = "location"
DEFAULT_DATE_FIELD = "date"

DEFAULT_ID_FIELD_INDEX = 0
DEFFAULT_ID_DELIMITER = "|"
DEFAULT_METADATA_DELIMITER = ","

DEFAULT_HEADER_SEPARATOR = "|"
DEFAULT_HEADER_FIELDS = "sample|location|date"

DEFAULT_N_THRESHOLD = 0.2
DEFAULT_CLUSTER_WINDOW = 10
DEFAULT_CLUSTER_COUNT = 3