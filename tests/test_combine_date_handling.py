#!/usr/bin/env python3
"""Tests for date handling and --require-date functionality in combine command."""

import csv
import os
import tempfile
import pytest
from raccoon.commands import combine


class TestHarmonizeDate:
    """Test date harmonization with various formats."""
    
    def test_iso_date_already_harmonized(self):
        """Test that valid ISO dates are recognized."""
        result, issue = combine.harmonize_date("2024-01-15")
        assert result == "2024-01-15"
        assert issue is None
    
    def test_year_month_day_with_slashes(self):
        """Test YYYY/MM/DD format."""
        result, issue = combine.harmonize_date("2024/01/15")
        assert result == "2024-01-15"
        assert issue is None
    
    def test_day_month_year_with_slashes(self):
        """Test DD/MM/YYYY format."""
        result, issue = combine.harmonize_date("15/01/2024")
        assert result == "2024-01-15"
        assert issue is None
    
    def test_year_month_only_with_hyphen(self):
        """Test YYYY-MM format (year-month only)."""
        result, issue = combine.harmonize_date("2024-01")
        assert result == "2024-01-01"  # Should use first day of month
        assert issue is None
    
    def test_year_month_only_with_slash(self):
        """Test YYYY/MM format."""
        result, issue = combine.harmonize_date("2024/01")
        assert result == "2024-01-01"
        assert issue is None
    
    def test_year_only_four_digits(self):
        """Test year-only format YYYY."""
        result, issue = combine.harmonize_date("2024")
        assert result == "2024-01-01"  # Should use January 1st
        assert issue is None
    
    def test_compact_date_format_yyyymmdd(self):
        """Test compact YYYYMMDD format."""
        result, issue = combine.harmonize_date("20240115")
        assert result == "2024-01-15"
        assert issue is None
    
    def test_compact_date_format_yyyymm(self):
        """Test compact YYYYMM format."""
        result, issue = combine.harmonize_date("202401")
        assert result == "2024-01-01"
        assert issue is None
    
    def test_invalid_date(self):
        """Test that invalid dates return None."""
        result, issue = combine.harmonize_date("invalid-date")
        assert result is None
        assert issue is not None
    
    def test_empty_string(self):
        """Test that empty strings return None."""
        result, issue = combine.harmonize_date("")
        assert result is None
        assert issue is None
    
    def test_none_input(self):
        """Test that None input returns None."""
        result, issue = combine.harmonize_date(None)
        assert result is None
        assert issue is None
    
    def test_date_with_dots(self):
        """Test YYYY.MM.DD format."""
        result, issue = combine.harmonize_date("2024.01.15")
        assert result == "2024-01-15"
        assert issue is None


class TestExtractDateFromText:
    """Test date extraction from arbitrary text (e.g., headers)."""
    
    def test_extract_iso_date_from_header(self):
        """Test extracting ISO date from FASTA header."""
        header = "sample_001|2024-01-15|location"
        result = combine.extract_date_from_text(header)
        assert result == "2024-01-15"
    
    def test_extract_compact_date_from_header(self):
        """Test extracting compact YYYYMMDD date."""
        header = "sample_20240115_extra"
        result = combine.extract_date_from_text(header)
        assert result == "2024-01-15"
    
    def test_extract_year_only_from_header(self):
        """Test extracting year-only date from header."""
        header = "sequence_2024_isolate"
        result = combine.extract_date_from_text(header)
        assert result == "2024-01-01"
    
    def test_no_date_in_text(self):
        """Test that text without dates returns None."""
        header = "sample_no_date_here"
        result = combine.extract_date_from_text(header)
        assert result is None
    
    def test_empty_text(self):
        """Test that empty text returns None."""
        result = combine.extract_date_from_text("")
        assert result is None
    
    def test_none_text(self):
        """Test that None input returns None."""
        result = combine.extract_date_from_text(None)
        assert result is None
    
    def test_extract_first_date_multiple_dates(self):
        """Test that first valid date is extracted."""
        header = "2024-01-01_and_2025-06-15"
        result = combine.extract_date_from_text(header)
        # Should return the first parseable date found
        assert result is not None
        assert result == "2024-01-01"


class TestParseHeaderTemplate:
    """Test header template parsing with field names."""
    
    def test_simple_template(self):
        """Test parsing simple template."""
        separator, fields = combine.parse_header_template("{id}|{location}|{date}")
        assert separator == "|"
        assert fields == ["id", "location", "date"]
    
    def test_template_with_spaces_in_field_names(self):
        """Test that field names with spaces are parsed correctly."""
        separator, fields = combine.parse_header_template("{sample}|{Country}|{Travel History}|{DOC}")
        assert separator == "|"
        assert fields == ["sample", "Country", "Travel History", "DOC"]
    
    def test_template_with_hyphen_separator(self):
        """Test template with different separator."""
        separator, fields = combine.parse_header_template("{id}-{location}-{date}")
        assert separator == "-"
        assert fields == ["id", "location", "date"]
    
    def test_template_with_spaces_in_fields(self):
        """Test complex field names."""
        separator, fields = combine.parse_header_template("{Sample ID}|{Collection Date}|{Geographic Location}")
        assert separator == "|"
        assert fields == ["Sample ID", "Collection Date", "Geographic Location"]
    
    def test_template_missing_placeholder(self):
        """Test that template without placeholders raises error."""
        with pytest.raises(ValueError, match="must contain at least one field placeholder"):
            combine.parse_header_template("no_placeholders_here")
    
    def test_template_inconsistent_separators(self):
        """Test that template with inconsistent separators raises error."""
        with pytest.raises(ValueError, match="single consistent separator"):
            combine.parse_header_template("{id}|{location}-{date}")


class TestRequireDateFiltering:
    """Test --require-date flag behavior in main function."""
    
    def test_require_date_with_valid_metadata(self):
        """Test sequences pass with valid dates in metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create FASTA file
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, "w") as f:
                f.write(">seq1\nACGT\n>seq2\nTGCA\n")
            
            # Create metadata with dates
            metadata_path = os.path.join(tmpdir, "metadata.csv")
            with open(metadata_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["id", "location", "date"])
                writer.writeheader()
                writer.writerow({"id": "seq1", "location": "USA", "date": "2024-01-15"})
                writer.writerow({"id": "seq2", "location": "Canada", "date": "2024-06"})
            
            # Mock args
            class Args:
                fasta = [fasta_path]
                metadata = [metadata_path]
                metadata_delimiter = ","
                metadata_id_field = "id"
                metadata_location_field = "location"
                metadata_date_field = "date"
                seq_id_delimiter = None
                seq_id_field_index = 0
                min_length = None
                max_n_content = None
                require_date = True
                header_fields = None
                header_separator = "|"
                outfile = os.path.join(tmpdir, "output.fasta")
                input_cmd_line = "test"
            
            result = combine.main(Args())
            assert result == 0
            
            # Check output contains both sequences
            with open(Args.outfile) as f:
                content = f.read()
                assert count_fasta_records(content) == 2
    
    def test_require_date_filters_missing_dates(self):
        """Test sequences are filtered when date is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, "w") as f:
                f.write(">seq1\nACGT\n>seq2\nTGCA\n")
            
            metadata_path = os.path.join(tmpdir, "metadata.csv")
            with open(metadata_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["id", "location", "date"])
                writer.writeheader()
                writer.writerow({"id": "seq1", "location": "USA", "date": "2024-01-15"})
                writer.writerow({"id": "seq2", "location": "Canada", "date": ""})  # Missing date
            
            class Args:
                fasta = [fasta_path]
                metadata = [metadata_path]
                metadata_delimiter = ","
                metadata_id_field = "id"
                metadata_location_field = "location"
                metadata_date_field = "date"
                seq_id_delimiter = None
                seq_id_field_index = 0
                min_length = None
                max_n_content = None
                require_date = True
                header_fields = None
                header_separator = "|"
                outfile = os.path.join(tmpdir, "output.fasta")
                input_cmd_line = "test"
            
            result = combine.main(Args())
            assert result == 0
            
            # Check output contains only seq1
            with open(Args.outfile) as f:
                content = f.read()
                assert count_fasta_records(content) == 1
                assert "seq1" in content
    
    def test_require_date_fallback_to_header(self):
        """Test that sequences without metadata row can pass if date is in header."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, "w") as f:
                # seq1 has no metadata, but has date in header
                f.write(">seq1_2024-01-15\nACGT\n>seq2\nTGCA\n")
            
            metadata_path = os.path.join(tmpdir, "metadata.csv")
            with open(metadata_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["id", "location", "date"])
                writer.writeheader()
                writer.writerow({"id": "seq2", "location": "Canada", "date": "2024-06"})
            
            class Args:
                fasta = [fasta_path]
                metadata = [metadata_path]
                metadata_delimiter = ","
                metadata_id_field = "id"
                metadata_location_field = "location"
                metadata_date_field = "date"
                seq_id_delimiter = None
                seq_id_field_index = 0
                min_length = None
                max_n_content = None
                require_date = True
                header_fields = None
                header_separator = "|"
                outfile = os.path.join(tmpdir, "output.fasta")
                input_cmd_line = "test"
            
            result = combine.main(Args())
            assert result == 0
            
            # Both sequences should pass
            with open(Args.outfile) as f:
                content = f.read()
                assert count_fasta_records(content) == 2
    
    def test_require_date_filters_no_date_anywhere(self):
        """Test sequences are filtered when no date in metadata or header."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, "w") as f:
                # seq1 has no date in header
                f.write(">seq1_no_date\nACGT\n")
            
            metadata_path = os.path.join(tmpdir, "metadata.csv")
            with open(metadata_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["id", "location", "date"])
                writer.writeheader()
                # seq1 not in metadata
            
            class Args:
                fasta = [fasta_path]
                metadata = [metadata_path]
                metadata_delimiter = ","
                metadata_id_field = "id"
                metadata_location_field = "location"
                metadata_date_field = "date"
                seq_id_delimiter = None
                seq_id_field_index = 0
                min_length = None
                max_n_content = None
                require_date = True
                header_fields = None
                header_separator = "|"
                outfile = os.path.join(tmpdir, "output.fasta")
                input_cmd_line = "test"
            
            result = combine.main(Args())
            assert result == 0
            
            # No sequences should pass
            with open(Args.outfile) as f:
                content = f.read()
                assert count_fasta_records(content) == 0
    
    def test_require_date_unparsable_date_filtered(self):
        """Test sequences are filtered with unparsable dates."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "test.fasta")
            with open(fasta_path, "w") as f:
                f.write(">seq1\nACGT\n")
            
            metadata_path = os.path.join(tmpdir, "metadata.csv")
            with open(metadata_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["id", "location", "date"])
                writer.writeheader()
                writer.writerow({"id": "seq1", "location": "USA", "date": "invalid-date-format"})
            
            class Args:
                fasta = [fasta_path]
                metadata = [metadata_path]
                metadata_delimiter = ","
                metadata_id_field = "id"
                metadata_location_field = "location"
                metadata_date_field = "date"
                seq_id_delimiter = None
                seq_id_field_index = 0
                min_length = None
                max_n_content = None
                require_date = True
                header_fields = None
                header_separator = "|"
                outfile = os.path.join(tmpdir, "output.fasta")
                input_cmd_line = "test"
            
            result = combine.main(Args())
            assert result == 0
            
            # Sequence should be filtered
            with open(Args.outfile) as f:
                content = f.read()
                assert count_fasta_records(content) == 0


def count_fasta_records(fasta_content: str) -> int:
    """Count the number of FASTA records in a string."""
    return fasta_content.count(">")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
