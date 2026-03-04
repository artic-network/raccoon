#!/usr/bin/env python3
"""Combine FASTA files and optionally harmonise headers from metadata."""
import csv
import logging
import os
import re
import sys
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO
from raccoon.utils import constants as rc


def parse_header_template(template: str) -> Tuple[str, List[str]]:
    """Parse a header template like '{id}|{location}|{date}' into separator and field names.
    
    Args:
        template: Header template string with placeholders like {fieldname}
        
    Returns:
        Tuple of (separator, field_names_list)
        e.g. ("|", ["id", "location", "date"])
    """
    # Find all {field} placeholders
    field_pattern = r'\{(\w+)\}'
    fields = re.findall(field_pattern, template)
    
    if not fields:
        raise ValueError(f"Header template must contain at least one field placeholder like {{id}}: {template}")
    
    # Determine separator by removing field placeholders
    separator = re.sub(field_pattern, '', template)
    
    if len(set(separator)) > 1:
        raise ValueError(f"Header template must use a single consistent separator between fields (e.g. '|'): {template}")

    return separator, fields


def format_header_from_template(
    template: str,
    field_values: Dict[str, str],
) -> str:
    """Format a header using a template string and field values.
    
    Args:
        template: Header template like '{id}|{location}|{date}'
        field_values: Dictionary mapping field names to their values
        
    Returns:
        Formatted header string with sanitized values
        (special characters replaced with underscores, except for ISO dates)
    """
    def _sanitize(value: str) -> str:
        if value is None or value == "":
            return ""
        cleaned = str(value).strip()
        lowered = cleaned.lower()
        
        # Skip sanitization for ISO dates (YYYY-MM-DD format)
        if re.match(r'^\d{4}-\d{2}-\d{2}$', lowered):
            return lowered
        
        # Replace specific characters that break Newick parsers: spaces, commas, colons, semi-colons, parentheses
        # Preserve hyphens as they are safe
        sanitized = re.sub(r'[ ,;:\(\)]', '_', lowered)
        # Replace multiple consecutive underscores with single underscore
        sanitized = re.sub(r'_+', '_', sanitized)
        # Remove leading/trailing underscores
        sanitized = sanitized.strip('_')
        return sanitized
    
    # Create sanitized field values
    sanitized = {k: _sanitize(v) for k, v in field_values.items()}
    
    # Format using template
    try:
        return template.format(**sanitized)
    except KeyError as e:
        raise ValueError(f"Missing field {e} in template data")



def harmonize_date(date_str: str) -> Tuple[Optional[str], Optional[str]]:
    """Harmonize a date string to ISO YYYY-MM-DD format.
    
    Args:
        date_str: date string to harmonize
        
    Returns:
        Tuple of (harmonized_date, issue_description or None)
        - If successful: (ISO_date_string, None)
        - If parsing failed: (None, error_message)
        - If ambiguous: (best_guess_date, "ambiguous: ...")
    """
    if not date_str or not isinstance(date_str, str):
        return None, None
    
    date_str = date_str.strip()
    if not date_str:
        return None, None
    
    # Already in ISO format?
    iso_pattern = r'^\d{4}-\d{2}-\d{2}$'
    if re.match(iso_pattern, date_str):
        try:
            datetime.strptime(date_str, '%Y-%m-%d')
            return date_str, None
        except ValueError:
            return None, f"invalid ISO date: {date_str}"
    
    # Try common formats
    common_formats = [
        '%Y/%m/%d',      # 2024/01/15
        '%d/%m/%Y',      # 15/01/2024
        '%m/%d/%Y',      # 01/15/2024
        '%Y-%m-%d',      # Already handled above, but keep for completeness
        '%d-%m-%Y',      # 15-01-2024
        '%m-%d-%Y',      # 01-15-2024
        '%Y.%m.%d',      # 2024.01.15
        '%d.%m.%Y',      # 15.01.2024
        '%B %d, %Y',     # January 15, 2024
        '%b %d, %Y',     # Jan 15, 2024
        '%d %B %Y',      # 15 January 2024
        '%d %b %Y',      # 15 Jan 2024
        '%Y%m%d',        # 20240115
    ]
    
    for fmt in common_formats:
        try:
            parsed = datetime.strptime(date_str, fmt)
            return parsed.strftime('%Y-%m-%d'), None
        except ValueError:
            continue
    
    # Couldn't parse with unambiguous formats, try ambiguous ones
    # Check for patterns like DD/MM or MM/DD where both are valid
    parts = re.split(r'[-/\.]', date_str)
    
    if len(parts) == 3:
        # Try to infer: if any part > 12, we know its position
        try:
            nums = [int(p) for p in parts]
        except ValueError:
            return None, f"failed to parse: {date_str}"
        
        # Pattern: num/num/YYYY
        if len(parts[2]) == 4:  # Last part is year
            year = nums[2]
            first, second = nums[0], nums[1]
            
            # If first or second > 12, we know it's the day
            if 12 < first <= 31:
                # first is day, second is month
                try:
                    parsed = datetime(year, second, first)
                    return parsed.strftime('%Y-%m-%d'), None
                except ValueError:
                    return None, f"invalid date values: {date_str}"
            elif 12 < second <= 31:
                # second is day, first is month
                try:
                    parsed = datetime(year, first, second)
                    return parsed.strftime('%Y-%m-%d'), None
                except ValueError:
                    return None, f"invalid date values: {date_str}"
            else:
                # Both <= 12, ambiguous - assume first is month, second is day
                try:
                    parsed = datetime(year, first, second)
                    return parsed.strftime('%Y-%m-%d'), f"ambiguous: interpreted {date_str} as {year}-{first:02d}-{second:02d}"
                except ValueError:
                    return None, f"invalid date values: {date_str}"
        
        # Pattern: YYYY/num/num (assume month/day)
        elif len(parts[0]) == 4:  # First part is year
            year = nums[0]
            month, day = nums[1], nums[2]
            try:
                parsed = datetime(year, month, day)
                return parsed.strftime('%Y-%m-%d'), None
            except ValueError:
                return None, f"invalid date values: {date_str}"
    
    return None, f"failed to parse: {date_str}"


def detect_date_field(
    header_fields_template: Optional[str],
    metadata_date_field_arg: Optional[str],
    template_field_names: List[str],
    metadata_columns: Optional[List[str]] = None,
) -> Optional[str]:
    """Detect which field is the date field.
    
    Priority:
    1. If --metadata-date-field was explicitly provided and not the default, use that
    2. If --header-fields provided, assume the final field is the date field
    3. Check if that field exists in metadata columns
    
    Args:
        header_fields_template: The --header-fields argument (or None)
        metadata_date_field_arg: The --metadata-date-field argument
        template_field_names: Parsed field names from header_fields_template
        metadata_columns: Available columns in metadata (for validation)
        
    Returns:
        The name of the date field to use, or None if no date field
    """
    # If header_fields_template is provided, assume the last field is the date
    if header_fields_template is not None and template_field_names:
        last_field = template_field_names[-1]
        # Validate the field exists in metadata if columns provided
        if metadata_columns is not None and last_field not in metadata_columns:
            return None  # Field doesn't exist, will be caught in validation
        return last_field
    
    # Otherwise, use the provided metadata_date_field
    return metadata_date_field_arg




def get_field(row: Dict[str, str], field: str) -> str:
    """Get a field from a CSV row, handling BOM markers."""
    if field in row:
        return row.get(field, "") or ""
    bom_field = f"\ufeff{field}"
    return row.get(bom_field, "") or ""


def load_metadata_map(metadata_path: str, id_field: str, delimiter: str) -> Optional[Tuple[Dict[str, Dict[str, str]], Optional[str]]]:
    """Load metadata from CSV file.
    
    Returns:
        Tuple of (metadata_dict, error_message)
        - If successful: (metadata_dict, None)
        - If parse error with unescaped delimiters: (None, error_message)
    """
    metadata = {}
    
    with open(metadata_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames
        expected_field_count = len(fieldnames) if fieldnames else 0
        
        for row_num, row in enumerate(reader, start=2):  # Start at 2 (after header)
            # Check if row has more fields than expected
            # This happens when there are unescaped delimiters in the data
            actual_field_count = len(row)
            if actual_field_count > expected_field_count:
                return None, (
                    f"Row {row_num} has {actual_field_count} fields but header has {expected_field_count}. "
                    f"This may indicate unescaped delimiters in the metadata. "
                    f"Please quote fields containing '{delimiter}' characters."
                )
            
            key = get_field(row, id_field)
            if not key:
                continue
            metadata[key] = row
    
    return metadata, None



def _infer_delimiter(path: str, default: str) -> str:
    lowered = path.lower()
    if lowered.endswith(".tsv") or lowered.endswith(".tab"):
        return "\t"
    return default


def load_metadata_maps(metadata_paths: Iterable[str], id_field: str, delimiter: str) -> Tuple[Optional[Dict[str, Dict[str, str]]], Optional[str]]:
    """Load metadata from multiple CSV files.
    
    Returns:
        Tuple of (metadata_dict, error_message)
        - If successful: (metadata_dict, None)
        - If error: (None, error_message)
    """
    merged = {}
    for path in metadata_paths:
        effective_delimiter = _infer_delimiter(path, delimiter)
        meta_map, error = load_metadata_map(path, id_field, effective_delimiter)
        if error:
            return None, f"Error in {path}: {error}"
        merged.update(meta_map)
    return merged, None


def format_header(
    record_id: str,
    row: Dict[str, str],
    location_field: str,
    date_field: str,
    sep: str,
) -> str:
    def _sanitize(value: str) -> str:
        if value is None:
            return ""
        cleaned = str(value).strip()
        lowered = cleaned.lower()
        return "_".join(lowered.split())

    location = _sanitize(get_field(row, location_field))
    date = _sanitize(get_field(row, date_field))
    return f"{record_id}{sep}{location}{sep}{date}"


def write_fasta_record(handle, header: str, sequence: str) -> None:
    sequence = "".join(sequence.split()).upper()
    handle.write(f">{header}\n{sequence}\n")


def n_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return seq.count("N") / len(seq)


def parse_record_id(record_id: str, delimiter: str, field_index: int) -> str:
    if field_index < 0:
        return record_id
    parts = record_id.split(delimiter) if delimiter else [record_id]
    if not parts:
        return record_id
    if field_index >= len(parts):
        return record_id
    return parts[field_index] or record_id


def main(args):
    """Combine fasta files into a single upper-case, unwrapped FASTA. Headers can be harmonised using metadata and field templates."""

    try:
        from raccoon.utils import io

        input_fastas = args.fasta or []
        if not input_fastas:
            logging.error("No input FASTA files provided")
            return 1

        for path in input_fastas:
            if not io.validate_input_file(path, "FASTA file"):
                return 1

        metadata_map: Optional[Dict[str, Dict[str, str]]] = None
        metadata_paths = list(getattr(args, "metadata", []) or [])
        metadata_delimiter = getattr(args, "metadata_delimiter", rc.DEFAULT_METADATA_DELIMITER)
        metadata_id_field = getattr(args, "metadata_id_field", rc.DEFAULT_ID_FIELD)
        metadata_location_field = getattr(args, "metadata_location_field", rc.DEFAULT_LOCATION_FIELD)
        metadata_date_field = getattr(args, "metadata_date_field", rc.DEFAULT_DATE_FIELD)
        seq_id_delimiter = getattr(args, "seq_id_delimiter", rc.DEFAULT_ID_DELIMITER)
        seq_id_field_index = getattr(args, "seq_id_field_index", rc.DEFAULT_ID_FIELD_INDEX)
        min_length = getattr(args, "min_length", None)
        max_n_content = getattr(args, "max_n_content", None)
        header_fields = getattr(args, "header_fields", None)
        header_separator = getattr(args, "header_separator", rc.DEFAULT_HEADER_SEPARATOR)

        # Determine header template and field names
        if header_fields_template is None:
            # onstruct template from existing parameters
            header_fields_template = f"{metadata_id_field}{header_separator}{metadata_location_field}{header_separator}{metadata_date_field}"
            template_field_names = [metadata_id_field, metadata_location_field, metadata_date_field]
        else:
            # Parse user-provided template
            try:
                _, template_field_names = parse_header_template(header_fields_template)
            except ValueError as e:
                logging.error(f"Invalid header template: {str(e)}")
                return 1

        if metadata_paths:
            for path in metadata_paths:
                if not io.validate_input_file(path, "Metadata file"):
                    return 1
            metadata_map, metadata_error = load_metadata_maps(metadata_paths, metadata_id_field, metadata_delimiter)
            if metadata_error:
                logging.error(f"Failed to load metadata: {metadata_error}")
                return 1
            
            # Get available metadata columns from the first row (for validation)
            metadata_columns = None
            if metadata_map:
                first_row = next(iter(metadata_map.values()))
                metadata_columns = list(first_row.keys())
                
                # If header_fields template is provided, validate all fields exist in metadata
                if header_fields_template is not None and getattr(args, "header_fields", None) is not None:
                    for field_name in template_field_names:
                        if field_name != "id" and field_name not in metadata_columns:
                            # Check for BOM-prefixed field name
                            bom_field = f"\ufeff{field_name}"
                            if bom_field not in metadata_columns:
                                logging.error(
                                    f"Field '{field_name}' in header template not found in metadata columns: {', '.join(metadata_columns)}"
                                )
                                return 1
            
            # Detect the date field
            detected_date_field = detect_date_field(
                getattr(args, "header_fields", None),
                metadata_date_field,
                template_field_names,
                metadata_columns,
            )
        else:
            detected_date_field = None

        output_path = args.output or "combined.fasta"
        if output_path == "-":
            out_handle = sys.stdout
            close_handle = False
        else:
            if not io.ensure_parent_directory(output_path):
                return 1
            out_handle = open(output_path, "w")
            close_handle = True

        filtered_count = 0
        kept_count = 0
        filter_failures = []
        metadata_issues = []
        try:
            for path in input_fastas:
                for rec in SeqIO.parse(path, "fasta"):
                    seq = str(rec.seq)
                    seq_len = len(seq)
                    n_prop = n_content(seq)
                    parsed_id = parse_record_id(rec.id, id_delimiter, id_field)
                    metadata_row = None
                    location_value = ""
                    date_value = ""
                    if metadata_map is not None:
                        metadata_row = metadata_map.get(parsed_id)
                        if metadata_row:
                            location_value = get_field(metadata_row, metadata_location_field).strip()
                            date_value = get_field(metadata_row, metadata_date_field).strip()
                        else:
                            logging.warning(f"No metadata row found for {parsed_id}")
                    reasons = []
                    if min_length is not None and seq_len < min_length:
                        reasons.append(f"length < {min_length}")
                    if max_n_content is not None and n_prop > max_n_content:
                        reasons.append(f"N content > {max_n_content}")
                    status = "filtered" if reasons else "kept"
                    if metadata_map is not None:
                        if metadata_row:
                            if not location_value:
                                metadata_issues.append({
                                    "file": os.path.basename(path),
                                    "id": rec.id,
                                    "parsed_id": parsed_id,
                                    "status": status,
                                    "issue": "missing location",
                                    "location": location_value,
                                    "date": date_value,
                                })
                            if not date_value:
                                metadata_issues.append({
                                    "file": os.path.basename(path),
                                    "id": rec.id,
                                    "parsed_id": parsed_id,
                                    "status": status,
                                    "issue": "missing date",
                                    "location": location_value,
                                    "date": date_value,
                                })
                        else:
                            metadata_issues.append({
                                "file": os.path.basename(path),
                                "id": rec.id,
                                "parsed_id": parsed_id,
                                "status": status,
                                "issue": "missing metadata row",
                                "location": "",
                                "date": "",
                            })
                    if reasons:
                        filter_failures.append({
                            "file": os.path.basename(path),
                            "id": rec.id,
                            "parsed_id": parsed_id,
                            "length": seq_len,
                            "n_content": round(n_prop, 6),
                            "reason": "; ".join(reasons),
                        })
                        filtered_count += 1
                        continue
                    header = rec.id
                    if metadata_map is not None:
                        row = metadata_map.get(parsed_id)
                        if row:
                            # Build field values dictionary for template formatting
                            field_values = {"id": parsed_id}
                            for field_name in template_field_names:
                                if field_name == "id":
                                    continue  # Already added
                                
                                # Get the raw value from metadata
                                raw_value = get_field(row, field_name)
                                
                                # If this is the date field, try to harmonize it first (before sanitization)
                                if detected_date_field and field_name == detected_date_field:
                                    harmonized, _ = harmonize_date(raw_value)
                                    field_values[field_name] = harmonized or raw_value
                                else:
                                    field_values[field_name] = raw_value
                            
                            try:
                                header = format_header_from_template(header_fields_template, field_values)
                            except ValueError as e:
                                logging.error(f"Failed to format header for {parsed_id}: {str(e)}")
                                continue
                    write_fasta_record(out_handle, header, seq)
                    kept_count += 1
        finally:
            if close_handle:
                out_handle.close()

        report_outdir = os.path.dirname(output_path) or os.getcwd()
        if not io.ensure_output_directory(report_outdir):
            return 1
        if min_length is not None or max_n_content is not None:
            filter_csv = os.path.join(report_outdir, "seq_qc_filter_failures.csv")
            with open(filter_csv, "w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["file", "id", "parsed_id", "length", "n_content", "reason"],
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(filter_failures)
        if metadata_map is not None:
            metadata_csv = os.path.join(report_outdir, "seq_qc_metadata_issues.csv")
            with open(metadata_csv, "w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["file", "id", "parsed_id", "status", "issue", "location", "date"],
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(metadata_issues)

        try:
            from raccoon.utils import reporting
            reporting.generate_combine_report(
                outdir=report_outdir,
                output_fasta=output_path if output_path != "-" else "",
                input_fastas=input_fastas,
                metadata_paths=metadata_paths or None,
                metadata_id_field=metadata_id_field,
                metadata_location_field=metadata_location_field,
                metadata_date_field=metadata_date_field,
                header_separator=header_separator,
                min_length=min_length,
                max_n_content=max_n_content,
                filter_failures=filter_failures,
                metadata_issues=metadata_issues,
            )
        except Exception:
            logging.exception("Failed to generate combine report")

        logging.info("Combined %d input FASTA files", len(input_fastas))
        if min_length is not None or max_n_content is not None:
            logging.info("Filtered %d sequences; kept %d", filtered_count, kept_count)
        return 0
    except Exception:
        logging.exception("Combine failed")
        return 2
