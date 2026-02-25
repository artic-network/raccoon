import csv
import textwrap
from pathlib import Path

from raccoon.commands import combine


class MockArgs:
    def __init__(self, **kwargs):
        self.inputs = kwargs.get("inputs", [])
        self.output = kwargs.get("output", None)
        self.metadata = kwargs.get("metadata", None)
        self.metadata_delimiter = kwargs.get("metadata_delimiter", ",")
        self.metadata_id_field = kwargs.get("metadata_id_field", "id")
        self.metadata_location_field = kwargs.get("metadata_location_field", "location")
        self.metadata_date_field = kwargs.get("metadata_date_field", "date")
        self.header_separator = kwargs.get("header_separator", "|")
        self.id_delimiter = kwargs.get("id_delimiter", "|")
        self.id_field = kwargs.get("id_field", 0)
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
    args = MockArgs(inputs=[str(f1), str(f2)], output=str(out))

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
            id,location,date
            seq1,UK,2024-01-01
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
    )

    result = combine.main(args)
    assert result == 0

    lines = out.read_text().strip().splitlines()
    assert lines[0] == ">seq1|UK|2024-01-01"
    assert lines[1] == "ACGT"
    assert lines[2] == ">seq2|US|2024-02-02"
    assert lines[3] == "GGGG"


def test_combine_examples_same_headers(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    examples = repo_root / "tests" / "test_data" / "combine"
    out = tmp_path / "combined_same_headers.fasta"

    args = MockArgs(
        inputs=[
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
        inputs=[
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
        ">A001|SiteA|2024-01-15",
        ">A002|SiteA|2024-01-20",
        ">A003|SiteC|2024-02-01",
        ">A004|SiteC|2024-02-15",
        ">B001|LocX|2023-11-05",
        ">B002|LocX|2023-12-12",
        ">B003|LocZ|2024-03-03",
        ">B004|LocZ|2024-03-10",
    ]
    assert headers == expected_headers


def test_combine_multiple_metadata_files(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("seq1", "acgt"), ("seq2", "gggg")])

    meta1 = tmp_path / "meta1.csv"
    meta1.write_text(
        textwrap.dedent(
            """\
            id,location,date
            seq1,UK,2024-01-01
            """
        )
    )
    meta2 = tmp_path / "meta2.csv"
    meta2.write_text(
        textwrap.dedent(
            """\
            id,location,date
            seq2,US,2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
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
            id,location,date
            sample1,UK,2024-01-01
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
        output=str(out),
        metadata=[str(metadata_path)],
        id_delimiter="|",
        id_field=0,
    )

    result = combine.main(args)
    assert result == 0
    text = out.read_text()
    assert text.startswith(">sample1|UK|2024-01-01")


def test_combine_id_field_out_of_range_keeps_full_id(tmp_path):
    fasta_path = tmp_path / "a.fasta"
    _write_fasta(fasta_path, [("sample1|Loc1|2024-01-01", "acgt")])

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
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
            id\tlocation\tdate
            seq1\tUSA\t2024-01-01
            seq2\t\t2024-02-02
            """
        )
    )

    out = tmp_path / "combined.fasta"
    args = MockArgs(
        inputs=[str(fasta_path)],
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
