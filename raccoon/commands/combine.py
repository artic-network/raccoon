#!/usr/bin/env python3
"""Combine FASTA files and optionally harmonise headers from metadata."""
import csv
import logging
import os
import sys
from typing import Dict, Iterable, Optional

from Bio import SeqIO


def get_field(row: Dict[str, str], field: str) -> str:
    if field in row:
        return row.get(field, "") or ""
    bom_field = f"\ufeff{field}"
    return row.get(bom_field, "") or ""


def load_metadata_map(metadata_path: str, id_field: str, delimiter: str) -> Dict[str, Dict[str, str]]:
    metadata = {}
    with open(metadata_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            key = get_field(row, id_field)
            if not key:
                continue
            metadata[key] = row
    return metadata


def _infer_delimiter(path: str, default: str) -> str:
    lowered = path.lower()
    if lowered.endswith(".tsv") or lowered.endswith(".tab"):
        return "\t"
    return default


def load_metadata_maps(metadata_paths: Iterable[str], id_field: str, delimiter: str) -> Dict[str, Dict[str, str]]:
    merged = {}
    for path in metadata_paths:
        effective_delimiter = _infer_delimiter(path, delimiter)
        merged.update(load_metadata_map(path, id_field, effective_delimiter))
    return merged


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
    """Combine fasta files into a single upper-case, unwrapped FASTA."""
    if not hasattr(args, "inputs"):
        raise ValueError("Expected argparse Namespace from raccoon.command.build_parser")

    try:
        from raccoon.utils import io

        inputs = args.inputs or []
        if not inputs:
            logging.error("No input FASTA files provided")
            return 1

        for path in inputs:
            if not io.validate_input_file(path, "FASTA file"):
                return 1

        metadata_map: Optional[Dict[str, Dict[str, str]]] = None
        metadata_paths = list(getattr(args, "metadata", []) or [])
        metadata_id_field = getattr(args, "metadata_id_field", "id")
        metadata_location_field = getattr(args, "metadata_location_field", "location")
        metadata_date_field = getattr(args, "metadata_date_field", "date")
        metadata_delimiter = getattr(args, "metadata_delimiter", ",")
        header_separator = getattr(args, "header_separator", "|")
        id_delimiter = getattr(args, "id_delimiter", "|")
        id_field = getattr(args, "id_field", 0)
        min_length = getattr(args, "min_length", None)
        max_n_content = getattr(args, "max_n_content", None)

        if metadata_paths:
            for path in metadata_paths:
                if not io.validate_input_file(path, "Metadata file"):
                    return 1
            metadata_map = load_metadata_maps(metadata_paths, metadata_id_field, metadata_delimiter)

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
            for path in inputs:
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
                            logging.warning("No metadata row found for %s", parsed_id)
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
                            header = format_header(
                                parsed_id,
                                row,
                                metadata_location_field,
                                metadata_date_field,
                                header_separator,
                            )
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
                input_fastas=inputs,
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

        logging.info("Combined %d input FASTA files", len(inputs))
        if min_length is not None or max_n_content is not None:
            logging.info("Filtered %d sequences; kept %d", filtered_count, kept_count)
        return 0
    except Exception:
        logging.exception("Combine failed")
        return 2
