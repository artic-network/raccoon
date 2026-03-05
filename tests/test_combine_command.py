import csv
import textwrap
from pathlib import Path

from raccoon.commands import combine


class MockArgs:
    def __init__(self, **kwargs):
        self.fasta = kwargs.get("fasta", [])
        self.output = kwargs.get("output", None)
        self.metadata = kwargs.get("metadata", None)
        self.metadata_delimiter = kwargs.get("metadata_delimiter", ",")
        self.metadata_id_field = kwargs.get("metadata_id_field", "sample")
        self.metadata_location_field = kwargs.get("metadata_location_field", "location")
        self.metadata_date_field = kwargs.get("metadata_date_field", "date")
        self.header_separator = kwargs.get("header_separator", "|")
        self.header_fields = kwargs.get("header_fields", None)
        self.seq_id_delimiter = kwargs.get("seq_id_delimiter", "|")
        self.seq_id_field_index = kwargs.get("seq_id_field_index", 0)
        self.min_length = kwargs.get("min_length", None)
        self.max_n_content = kwargs.get("max_n_content", None)


def _write_fasta(path, entries):
    lines = []
    for header, seq in entries:
        lines.append(f">{header}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")


def test_combine_uppercase_and_unwrapped(tmp_path):
    f1 = tmp_path / "a.fasta"
    f2 = tmp_path / "b.fasta"

    _write_fasta(f1, [("seq1", "acgtacgt"), ("seq2", "aaaa")])
    _write_fasta(f2, [("seq3", "tttt")])

    out = tmp_path / "combined.fasta"
    args = MockArgs(fasta=[str(f1), str(f2)], output=str(out))

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1"
    assert lines[1] == "ACGTACGT"
    assert lines[2] == ">seq2"
    assert lines[3] == "AAAA"
    assert lines[4] == ">seq3"
    assert lines[5] == "TTTT"


def test_combine_harmonises_headers_with_metadata(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|uk|2024-01-01"
    assert lines[1] == "ACGT"
    assert lines[2] == ">seq2|us|2024-02-02"
    assert lines[3] == "GGGG"


def test_combine_examples_same_headers(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    examples = repo_root / "tests" / "test_data" / "combine"
    out = tmp_path / "combined_same_headers.fasta"

    args = MockArgs(
        fasta=[
            str(examples / "inputs" / "set_a.fasta"),
            str(examples / "inputs" / "set_b.fasta"),
        ],
        output=str(out),
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert len(lines) == 16
    headers = lines[0::2]
    expected_headers = [
        ">A001",
        ">A002",
        ">A003",
        ">A004",
        ">B001",
        ">B002",
        ">B003",
        ">B004",
    ]
    assert headers == expected_headers
    for sequence in lines[1::2]:
        assert sequence == sequence.upper()
        assert " " not in sequence


def test_combine_examples_harmonised_headers(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    examples = repo_root / "tests" / "test_data" / "combine"
    out = tmp_path / "combined_harmonised_headers.fasta"

    args = MockArgs(
        fasta=[
            str(examples / "inputs" / "set_a.fasta"),
            str(examples / "inputs" / "set_b.fasta"),
        ],
        output=str(out),
        metadata=[str(examples / "metadata.csv")],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert len(lines) == 16
    headers = lines[0::2]
    expected_headers = [
        ">A001|sitea|2024-01-15",
        ">A002|sitea|2024-01-20",
        ">A003|sitec|2024-02-01",
        ">A004|sitec|2024-02-15",
        ">B001|locx|2023-11-05",
        ">B002|locx|2023-12-12",
        ">B003|locz|2024-03-03",
        ">B004|locz|2024-03-10",
    ]
    assert headers == expected_headers


def test_combine_multiple_metadata_files(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    meta1 = tmp_path / "meta1.csv"
    meta1.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            """
        )
    )
    meta2 = tmp_path / "meta2.csv"
    meta2.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(meta1), str(meta2)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|uk|2024-01-01"
    assert lines[1] == "ACGT"
    assert lines[2] == ">seq2|us|2024-02-02"
    assert lines[3] == "GGGG"


def test_combine_parses_id_from_header(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("sample1|Loc1|2024-01-01", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            sample1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        id_delimiter="|",
        id_field=0,
    )

    result = combine.main(args)
    assert result == 0
    text = out.read_text()
    assert text.startswith(">sample1|uk|2024-01-01")


def test_combine_id_field_out_of_range_keeps_full_id(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("sample1|Loc1|2024-01-01", "acgt")])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        id_delimiter="|",
        id_field=10,
    )

    result = combine.main(args)
    assert result == 0
    text = out.read_text()
    assert text.startswith(">sample1|Loc1|2024-01-01")


def _read_csv_rows(path: Path):
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(row)
    return rows


def test_combine_writes_filter_and_metadata_issue_csvs(tmp_path):
    fasta_path = tmp_path / "inputs.fasta"
    _write_fasta(
        fasta_path,
        [
            ("seq1", "AAAA"),
            ("seq2", "NNNN"),
            ("seq3", "ACGTACGT"),
        ],
    )

    metadata_path = tmp_path / "meta.tsv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample\tlocation\tdate
            seq1\tUSA\t2024-01-01
            seq2\t\t2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_delimiter=",",
        min_length=5,
        max_n_content=0.2,
    )

    result = combine.main(args)
    assert result == 0

    filter_csv = tmp_path / "seq_qc_filter_failures.csv"
    metadata_csv = tmp_path / "seq_qc_metadata_issues.csv"
    assert filter_csv.exists()
    assert metadata_csv.exists()

    filter_rows = _read_csv_rows(filter_csv)
    filter_by_id = {row["id"]: row for row in filter_rows}
    assert set(filter_by_id) == {"seq1", "seq2"}
    assert filter_by_id["seq1"]["reason"] == "length < 5"
    assert filter_by_id["seq2"]["reason"] == "length < 5; N content > 0.2"

    metadata_rows = _read_csv_rows(metadata_csv)
    issues = {(row["id"], row["issue"], row["status"]) for row in metadata_rows}
    assert ("seq2", "missing location", "filtered") in issues
    assert ("seq3", "missing metadata row", "kept") in issues


def test_header_template_basic(tmp_path):
    """Test the new header template feature with default-like template."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        header_fields="{sample}|{location}|{date}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|uk|2024-01-01"
    assert lines[1] == "ACGT"
    assert lines[2] == ">seq2|us|2024-02-02"
    assert lines[3] == "GGGG"


def test_header_template_different_order(tmp_path):
    """Test header template with fields in different order."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        header_fields="{date}_{location}_{sample}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">2024-01-01_uk_seq1"
    assert lines[1] == "ACGT"
    assert lines[2] == ">2024-02-02_us_seq2"
    assert lines[3] == "GGGG"


def test_header_template_custom_separator(tmp_path):
    """Test header template with custom separator."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        header_fields="{sample}:{location}:{date}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1:uk:2024-01-01"


def test_header_template_subset_fields(tmp_path):
    """Test header template with only a subset of available fields."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        header_fields="{sample}|{location}",  # No date
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|uk"


def test_header_template_mismatch_id_fields(tmp_path):
    """Test header template with mismatching custom metadata id names."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample_id,country,collection_date
            seq1,United Kingdom,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="sample",
        header_fields="{sample_id}|{country}|{collection_date}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1"


def test_header_template_custom_metadata_fields(tmp_path):
    """Test header template with custom metadata field names."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample_id,country,collection_date
            seq1,United Kingdom,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="sample_id",
        header_fields="{sample_id}|{country}|{collection_date}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|united_kingdom|2024-01-01"


def test_backward_compatibility_no_header_fields(tmp_path):
    """Test that existing behavior works without --header-fields."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    # Note: NOT providing header_fields, should use default
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should produce default format: id|location|date
    assert lines[0] == ">seq1|uk|2024-01-01"
    assert lines[2] == ">seq2|us|2024-02-02"


# ========== EDGE CASE / CONFLICTING SIGNAL TESTS ==========


def test_metadata_args_without_metadata_file(tmp_path):
    """Test providing metadata args (delimiter, id-field) but no metadata file.
    Should work fine - just treated as if no metadata provided."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    out = tmp_path / "combined.fasta"
    # Provide metadata args but NOT the actual metadata file
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=None,  # No metadata file
        metadata_delimiter="\t",  # These args are ignored
        metadata_id_field="sample_id",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should output just the original header (no metadata applied)
    assert lines[0] == ">seq1"


def test_metadata_missing_location_field(tmp_path):
    """Test metadata file missing the specified location field.
    Should still work, but location values will be empty."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    # No 'location' column, only 'sample' and 'date'
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,date
            seq1,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_location_field="location",  # This field doesn't exist
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # location is empty, but should still include it in header
    assert lines[0] == ">seq1||2024-01-01"


def test_metadata_missing_date_field(tmp_path):
    """Test metadata file missing the specified date field.
    Should still work, but date values will be empty."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    # No 'date' column, only 'sample' and 'location'
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location
            seq1,UK
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_date_field="date",  # This field doesn't exist
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # date is empty, but location is present
    assert lines[0] == ">seq1|uk|"


def test_metadata_both_location_and_date_missing(tmp_path):
    """Test metadata file missing both location and date fields."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    # Only id column present
    metadata_path.write_text(
        textwrap.dedent(
            """\
            id
            seq1
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="id",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Both location and date are empty
    assert lines[0] == ">seq1||"


def test_header_fields_takes_precedence_over_location_date_args(tmp_path):
    """Test that header-fields template takes precedence over metadata-location-field
    and metadata-date-field arguments."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date,country,collection_date
            seq1,UK,2024-01-01,Great Britain,Jan 2024
            """
        )
    )

    out = tmp_path / "combined.fasta"
    # Provide both old args AND new header-fields arg
    # header-fields should take precedence
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_location_field="location",
        metadata_date_field="date",
        header_fields="{sample}|{country}|{collection_date}",  # Uses different fields
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should use country and collection_date, not location and date
    assert lines[0] == ">seq1|great_britain|jan_2024"


def test_header_fields_with_nonexistent_metadata_field(tmp_path):
    """Test header-fields template referencing a field that doesn't exist in metadata.
    Should exit with an error and informative message."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location
            seq1,UK
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        # Template references 'date' field that doesn't exist
        header_fields="{sample}|{location}|{date}",
    )

    result = combine.main(args)
    # Should exit with error because 'date' field doesn't exist in metadata
    assert result == 1
def test_header_fields_invalid_template_syntax(tmp_path):
    """Test that invalid template syntax is caught."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location
            seq1,UK
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        header_fields="invalid template with no placeholders",  # No {} placeholders
    )

    result = combine.main(args)
    # Should error
    assert result == 1


def test_header_fields_with_mixed_custom_and_standard_fields(tmp_path):
    """Test header-fields using both standard fields (sample, date) and custom fields."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date,region,host
            seq1,UK,2024-01-01,Europe,human
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        header_fields="{host}_{region}_{sample}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">human_europe_seq1"


def test_conflicting_metadata_delimiter_and_tsv_extension(tmp_path):
    """Test metadata-delimiter arg vs auto-detection from .tsv extension.
    Auto-detection should take precedence for .tsv files."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.tsv"  # .tsv extension triggers auto-detection
    metadata_path.write_text(
        "sample\tlocation\tdate\nseq1\tUK\t2024-01-01\n"
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_delimiter=",",  # User specifies comma, but file is TampleSV
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should correctly parse TSV despite comma being specified
    assert lines[0] == ">seq1|uk|2024-01-01"


def test_metadata_id_field_with_pipe_parsing(tmp_path):
    """Test that metadata-id-field works correctly in conjunction with id-field
    and id-delimiter for header parsing."""
    fasta_path = tmp_path / "a.fasta"
    # Header has pipes, we'll extract field 0
    _write_fasta(fasta_path, [("sample1|extra|data", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample_id,location,date
            sample1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="sample_id",  # Metadata uses different column name
        id_delimiter="|",
        id_field=0,  # Extract first field from header
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should have matched via parsed ID
    assert lines[0] == ">sample1|uk|2024-01-01"


def test_header_fields_id_special_mapping(tmp_path):
    """Test that {sample} in header template always maps to parsed ID,
    not to metadata id column."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("MYSEQ", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            seq_name,location
            MYSEQ,UK
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="seq_name",
        header_fields="{seq_name}|{location}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # {seq_name} should be the parsed ID (MYSEQ), and location from metadata
    assert lines[0] == ">myseq|uk"


def test_empty_metadata_field_values(tmp_path):
    """Test handling of explicitly empty values in metadata."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    # Empty location for seq1, empty date for seq2
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,,2024-01-01
            seq2,UK,
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1||2024-01-01"  # Empty location
    assert lines[2] == ">seq2|uk|"  # Empty date


def test_header_fields_only_id(tmp_path):
    """Test header template with only id field (no metadata refs)."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=None,  # No metadata
        header_fields="{sample}",  # Just the sample ID, no metadata fields
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1"
    assert lines[2] == ">seq2"


def test_metadata_args_present_but_empty_list(tmp_path):
    """Test that empty metadata list doesn't cause issues."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[],  # Empty list, not None
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1"


def test_date_harmonization_with_header_fields(tmp_path):
    """Test that dates are harmonized to ISO format when using header-fields template."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,sample_date
            seq1,UK,15/01/2024
            seq2,US,"January 20, 2024"
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        # Template with custom date field name - final field should be treated as date
        header_fields="{sample}|{location}|{sample_date}",
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Dates should be harmonized to ISO YYYY-MM-DD format, then sanitized (lowercase)
    # 15/01/2024 -> 2024-01-15, January 20, 2024 -> 2024-01-20
    assert lines[0] == ">seq1|uk|2024-01-15"
    assert lines[2] == ">seq2|us|2024-01-20"


def test_date_harmonization_default_date_field(tmp_path):
    """Test that dates are harmonized using the default 'date' field when --header-fields is not used."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024/01/15
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        # No header_fields template - use defaults
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # 2024/01/15 should be harmonized to 2024-01-15
    assert lines[0] == ">seq1|uk|2024-01-15"


def test_csv_unescaped_delimiters_error(tmp_path):
    """Test that unescaped delimiters in metadata produce informative error."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    # This CSV has an unescaped comma in the location field
    # It will be interpreted as 4 columns instead of 3
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,New York, USA,2024-01-15
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    # Should exit with error due to inconsistent field count
    assert result == 1


def test_csv_properly_quoted_delimiters(tmp_path):
    """Test that properly quoted fields with delimiters work correctly."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    # Proper CSV with quoted fields containing delimiters
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,"New York, USA",2024-01-15
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # "New York, USA" -> "new_york_usa" (spaces and special chars to underscores)
    assert lines[0] == ">seq1|new_york_usa|2024-01-15"


def test_sanitization_preserves_hyphens(tmp_path):
    """Test that hyphens in metadata are preserved during sanitization."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,New-York,2024-01-15
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Hyphens should be preserved, not converted to underscores
    assert lines[0] == ">seq1|new-york|2024-01-15"


def test_sanitization_special_chars_to_underscores(tmp_path):
    """Test that special characters (commas, colons, etc.) are converted to underscores."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,"New York; USA: Region (Northeast)",2024-01-15
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Special chars should be converted to underscores, but hyphens in dates preserved
    assert lines[0] == ">seq1|new_york_usa_region_northeast|2024-01-15"


# ============================================================================
# NEW EDGE CASE TESTS - Additional coverage for missing scenarios
# ============================================================================


def test_min_length_filtering(tmp_path):
    """Test --min-length filtering excludes short sequences."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgtacgtacgt"),  # 12 bp
        ("seq2", "acgt"),           # 4 bp
        ("seq3", "acgtacgt"),       # 8 bp
    ])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        min_length=10,  # Only keep sequences >= 10 bp
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should only have seq1
    assert ">seq1" in lines
    assert ">seq2" not in lines  # 4 bp < 10 bp
    assert ">seq3" not in lines  # 8 bp < 10 bp


def test_max_n_content_filtering(tmp_path):
    """Test --max-n-content filtering excludes sequences with too many Ns."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgtacgtacgt"),      # 0% N
        ("seq2", "nnnnnnnnacgt"),       # 80% N (8 of 10)
        ("seq3", "acgtnacgtnacgt"),     # 20% N (2 of 10)
    ])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        max_n_content=0.5,  # Only keep sequences with <= 50% N
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should have seq1 and seq3, but not seq2
    assert ">seq1" in lines
    assert ">seq2" not in lines  # 80% N > 50%
    assert ">seq3" in lines      # 20% N <= 50%


def test_min_length_and_max_n_content_combined(tmp_path):
    """Test that both length and N-content filters work together."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgtacgtacgtacgt"),      # 16 bp, 0% N - KEEP
        ("seq2", "acgt"),                  # 4 bp, 0% N - filtered by length
        ("seq3", "nnnnnnnnnnnnnnactivgt"),  # ~80% N - filtered by N-content
        ("seq4", "acgtnacgt"),             # 8 bp, 12.5% N - filtered by length
    ])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        min_length=10,
        max_n_content=0.5,
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert ">seq1" in lines
    assert ">seq2" not in lines
    assert ">seq3" not in lines
    assert ">seq4" not in lines


def test_id_not_in_metadata_is_included_with_empty_values(tmp_path):
    """Test that sequences with IDs not in metadata are included with empty field values."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgt"),
        ("seq2", "gggg"),  # Not in metadata
        ("seq3", "tttt"),
    ])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq3,US,2024-01-15
            """
        )
    )

    out = tmp_path / "combined.fasta"
    metadata_issues = tmp_path / "seq_qc_metadata_issues.csv"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # seq2 should have empty location and date fields
    assert ">seq2||" in lines


def test_duplicate_ids_in_metadata_uses_first_match(tmp_path):
    """Test that when multiple metadata rows have the same ID, the first match is used."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq1,US,2024-02-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should use first match (UK, not US)
    assert lines[0] == ">seq1|uk|2024-01-01"


def test_case_sensitivity_in_id_matching(tmp_path):
    """Test that ID matching is case-sensitive by default."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgt"),    # lowercase
        ("SEQ2", "gggg"),    # uppercase
    ])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-01-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # seq1 matches (case-sensitive), SEQ2 doesn't match seq2
    assert lines[0] == ">seq1|uk|2024-01-01"
    # SEQ2 should have empty fields since it doesn't match seq2
    assert lines[2] == ">seq2||"


def test_whitespace_in_headers_preserved(tmp_path):
    """Test handling of whitespace in FASTA headers."""
    fasta_path = tmp_path / "a.fasta"
    # FASTA headers with leading/trailing whitespace
    fasta_path.write_text(">  seq1  \nacgt\n")

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Whitespace should be preserved in the header value
    assert ">  seq1  " in lines


def test_whitespace_in_metadata_values(tmp_path):
    """Test that leading/trailing whitespace in metadata is handled."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    # Metadata with whitespace (CSV parser should handle this)
    metadata_path.write_text(
        'sample,location,date\nseq1," UK ",2024-01-01\n'
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Whitespace around value should be stripped by CSV parser, then sanitized
    assert lines[0] == ">seq1|uk|2024-01-01"


def test_missing_metadata_id_column(tmp_path):
    """Test error handling when metadata file doesn't have the specified ID column."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            id,location,date
            seq1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="sample",  # Asking for 'sample' column which doesn't exist
    )

    result = combine.main(args)
    # Should error because 'sample' column not found
    assert result == 1


def test_empty_fasta_file(tmp_path):
    """Test handling of empty FASTA file."""
    fasta_path = tmp_path / "empty.fasta"
    fasta_path.write_text("")

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
    )

    result = combine.main(args)
    assert result == 0

    # Output should be empty
    assert out.read_text() == ""


def test_single_sequence(tmp_path):
    """Test combining with just a single sequence."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgtacgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            sample,location,date
            seq1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert len(lines) == 2
    assert lines[0] == ">seq1|uk|2024-01-01"
    assert lines[1] == "ACGTACGT"


def test_duplicate_sequence_ids_in_same_file(tmp_path):
    """Test handling of duplicate sequence IDs within a single input file."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgt"),
        ("seq1", "gggg"),  # Duplicate ID
        ("seq2", "tttt"),
    ])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Both seq1 entries should be included
    header_count = sum(1 for line in lines if line.startswith(">"))
    assert header_count == 3


def test_unicode_characters_in_metadata(tmp_path):
    """Test handling of Unicode characters in metadata fields."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt")])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        'sample,location,date\nseq1,São Paulo,2024-01-01\n',
        encoding='utf-8'
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Unicode characters should be replaced during sanitization
    assert lines[0] == ">seq1|s_o_paulo|2024-01-01"


def test_sequence_id_with_special_delimiter_characters(tmp_path):
    """Test parsing sequence IDs when header contains the delimiter character."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("sample1|extra|info", "acgt"),  # Header with pipes
        ("sample2|data", "gggg"),
    ])

    metadata_path = tmp_path / "meta.csv"
    metadata_path.write_text(
        textwrap.dedent(
            """\
            id,location,date
            sample1,UK,2024-01-01
            sample2,US,2024-01-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        metadata_id_field="id",
        seq_id_delimiter="|",
        seq_id_field_index=0,  # Take first field
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">sample1|uk|2024-01-01"
    assert lines[2] == ">sample2|us|2024-01-02"


def test_zero_length_sequence(tmp_path):
    """Test handling of zero-length sequences."""
    fasta_path = tmp_path / "a.fasta"
    fasta_path.write_text(">seq1\n>seq2\nacgt\n")  # seq1 has no sequence data

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    # Should include both sequences, even empty one
    assert ">seq1" in lines
    assert ">seq2" in lines


def test_all_sequences_filtered_out(tmp_path):
    """Test behavior when all sequences are filtered by length/N-content."""
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [
        ("seq1", "acgt"),      # 4 bp
        ("seq2", "nnnnnnnn"),  # 100% N
    ])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        fasta=[str(fasta_path)],
        output=str(out),
        min_length=10,
        max_n_content=0.5,
    )

    result = combine.main(args)
    assert result == 0

    # Output file should be empty or just have headers
    content = out.read_text().strip()
    assert content == "" or content.count(">") == 0
