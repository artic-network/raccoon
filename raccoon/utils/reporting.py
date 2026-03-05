"""HTML report generation for raccoon commands."""
from __future__ import annotations

import base64
import math
import numpy as np
import csv
import os
import platform
import sys
from datetime import datetime
from typing import Iterable, Optional, Dict, Any, List

from Bio import SeqIO
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import plot
from jinja2 import Environment, FileSystemLoader, select_autoescape

from .reconstruction_functions import load_tree, ensure_node_label
from .plotly_baltic import build_tree_plot


def _svg_data_uri(path: str) -> str:
    try:
        with open(path, "rb") as handle:
            data = handle.read()
        encoded = base64.b64encode(data).decode("ascii")
        return f"data:image/svg+xml;base64,{encoded}"
    except Exception:
        return ""


def _logo_html(data_uri: str, css_class: str, alt_text: str) -> str:
    if data_uri:
        return f'<img class="{css_class}" src="{data_uri}" alt="{alt_text}" />'
    return f'<div class="{css_class}" aria-label="{alt_text}"></div>'


def _render_html(template_name: str, context: Dict[str, Any]) -> str:
    templates_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "templates"))
    env = Environment(
        loader=FileSystemLoader(templates_dir),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template(template_name)
    return template.render(**context)


def _write_html(outfile: str, title: str, template_name: str, context: Dict[str, Any]) -> None:
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    raccoon_logo = _svg_data_uri(os.path.join(base_dir, "raccoon", "assets", "raccoon_logo.svg"))
    artic_logo = _svg_data_uri(os.path.join(base_dir, "raccoon", "assets", "artic-logo-small.svg"))
    raccoon_logo_html = _logo_html(raccoon_logo, "logo", "Raccoon logo")
    artic_logo_html = _logo_html(artic_logo, "logo-small", "ARTIC Network logo")
    base_context = {
        "title": title,
        "generated_stamp": context.get("generated_stamp", datetime.now().strftime("%Y-%m-%d %H:%M")),
        "raccoon_logo": raccoon_logo,
        "raccoon_logo_html": raccoon_logo_html,
        "artic_logo_html": artic_logo_html,
    }
    merged_context = {**base_context, **context}
    html = _render_html(template_name, merged_context)
    with open(outfile, "w") as handle:
        handle.write(html)


def _table_context(df: Optional[pd.DataFrame]) -> Optional[Dict[str, Any]]:
    if df is None or df.empty:
        return None
    headers = list(df.columns)
    rows = df.astype(str).values.tolist()
    return {"headers": headers, "rows": rows}


def _n_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return seq.count("N") / len(seq)


def _apply_plot_style(fig: go.Figure) -> None:
    fig.update_layout(
        font=dict(family="Helvetica Neue", size=12, color="#111"),
        colorway=["#4BA3A8", "#7A6BB1", "#D08BA8"],
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=True,
        hoverlabel=dict(font=dict(color="#ffffff"), bordercolor="#ffffff"),
    )
    fig.update_xaxes(showline=True, linecolor="black", linewidth=1, gridcolor="rgba(0,0,0,0.05)")
    fig.update_yaxes(showline=True, linecolor="black", linewidth=1, gridcolor="rgba(0,0,0,0.05)")


def _plot_div(fig: go.Figure, div_id: Optional[str] = None) -> str:
    return pio.to_html(
        fig,
        include_plotlyjs="cdn",
        full_html=False,
        config={
            "displayModeBar": False,
            "responsive": True,
            "showTips": False,
            "doubleClick": False,
        },
        div_id=div_id,
    )




def _safe_mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def _safe_min(values: Iterable[int]) -> int:
    values = list(values)
    return min(values) if values else 0


def _safe_max(values: Iterable[int]) -> int:
    values = list(values)
    return max(values) if values else 0


def generate_combine_report(
    outdir: str,
    output_fasta: str,
    input_fastas: Iterable[str],
    metadata_paths: Optional[Iterable[str]] = None,
    metadata_id_field: str = "id",
    metadata_location_field: str = "location",
    metadata_date_field: str = "date",
    header_separator: str = "|",
    min_length: Optional[int] = None,
    max_n_content: Optional[float] = None,
    filter_failures: Optional[List[Dict[str, Any]]] = None,
    metadata_issues: Optional[List[Dict[str, Any]]] = None,
) -> str:
    records_summary = []
    ids_by_file = {}
    seq_details_by_file: Dict[str, List[Dict[str, Any]]] = {}
    filtered_rows: List[Dict[str, Any]] = []
    for path in input_fastas:
        lengths = []
        n_contents = []
        ids = []
        seq_details = []
        for rec in SeqIO.parse(path, "fasta"):
            seq = str(rec.seq)
            seq_len = len(seq)
            n_content = _n_content(seq)
            lengths.append(seq_len)
            n_contents.append(n_content)
            ids.append(rec.id)
            seq_details.append({
                "id": rec.id,
                "length": seq_len,
                "n_content": n_content,
            })
            reasons = []
            if min_length is not None and seq_len < min_length:
                reasons.append(f"length < {min_length}")
            if max_n_content is not None and n_content > max_n_content:
                reasons.append(f"N content > {max_n_content}")
            if reasons:
                filtered_rows.append({
                    "file": os.path.basename(path),
                    "id": rec.id,
                    "length": seq_len,
                    "n_content": round(n_content, 4),
                    "reason": "; ".join(reasons),
                })
        records_summary.append({
            "file": os.path.basename(path),
            "sequences": len(lengths),
            "len_min": _safe_min(lengths),
            "len_max": _safe_max(lengths),
            "len_mean": round(_safe_mean(lengths), 2),
            "n_min": round(_safe_min(n_contents), 4),
            "n_max": round(_safe_max(n_contents), 4),
            "n_mean": round(_safe_mean(n_contents), 4),
        })
        ids_by_file[os.path.basename(path)] = ids
        seq_details_by_file[os.path.basename(path)] = seq_details

    metadata_summary = "No metadata used."
    metadata_tables: List[Dict[str, Any]] = []
    metadata_locations = pd.Series(dtype=str)
    metadata_dates = pd.Series(dtype="datetime64[ns]")
    def _infer_delimiter(path: str, fallback: str = ",") -> str:
        lowered = path.lower()
        if lowered.endswith(".tsv") or lowered.endswith(".tab"):
            return "\t"
        return fallback

    if metadata_paths:
        try:
            frames = []
            for path in metadata_paths:
                frame = pd.read_csv(path, sep=_infer_delimiter(path))
                frames.append(frame)
                cols = list(frame.columns)
                metadata_tables.append({
                    "title": os.path.basename(path),
                    "row_count": len(frame),
                    "headers": cols,
                    "rows": frame.astype(str).values.tolist(),
                })
            meta = pd.concat(frames, ignore_index=True)
            metadata_locations = meta.get(metadata_location_field, pd.Series(dtype=str)).dropna().astype(str)
            metadata_dates = pd.to_datetime(meta.get(metadata_date_field, pd.Series(dtype=str)), errors="coerce")
            date_min = metadata_dates.min()
            date_max = metadata_dates.max()
            date_range = "n/a"
            if pd.notna(date_min) and pd.notna(date_max):
                date_range = f"{date_min.date().isoformat()} → {date_max.date().isoformat()}"
            metadata_summary = (
                f"Metadata files: {', '.join([os.path.basename(p) for p in metadata_paths])}. "
                f"ID field: {metadata_id_field}. "
                f"Location field: {metadata_location_field} (unique: {metadata_locations.nunique()}). "
                f"Date field: {metadata_date_field} (range: {date_range})."
            )
        except Exception:
            metadata_summary = "Metadata provided, but summary could not be parsed."

    header_stats = {"locations": 0, "dates": ""}
    try:
        if output_fasta:
            locs = set()
            dates = []
            for rec in SeqIO.parse(output_fasta, "fasta"):
                parts = rec.id.split(header_separator)
                if len(parts) >= 2 and parts[1]:
                    locs.add(parts[1])
                if len(parts) >= 3:
                    try:
                        parsed = pd.to_datetime(parts[2])
                        if pd.notna(parsed):
                            dates.append(parsed)
                    except Exception:
                        pass
            header_stats["locations"] = len(locs)
            if dates:
                dmin = min(dates)
                dmax = max(dates)
                header_stats["dates"] = f"{dmin.date().isoformat()} → {dmax.date().isoformat()}"
        elif metadata_paths and not metadata_locations.empty:
            header_stats["locations"] = int(metadata_locations.nunique())
            if not metadata_dates.empty:
                date_min = metadata_dates.min()
                date_max = metadata_dates.max()
                if pd.notna(date_min) and pd.notna(date_max):
                    header_stats["dates"] = f"{date_min.date().isoformat()} → {date_max.date().isoformat()}"
    except Exception:
        pass

    dataset_plot_html = ""
    date_location_records = []
    if output_fasta:
        try:
            for rec in SeqIO.parse(output_fasta, "fasta"):
                parts = rec.id.split(header_separator)
                if len(parts) < 3:
                    continue
                location = parts[1].strip()
                date_value = parts[2].strip()
                if not location or not date_value:
                    continue
                parsed_date = pd.to_datetime(date_value, errors="coerce")
                if pd.isna(parsed_date):
                    continue
                date_location_records.append({
                    "date": parsed_date,
                    "location": location,
                    "id": parts[0],
                })
        except Exception:
            date_location_records = []
    elif output_fasta:
        try:
            for rec in SeqIO.parse(output_fasta, "fasta"):
                parts = rec.id.split(header_separator)
                if len(parts) >= 3:
                    date_value = pd.to_datetime(parts[2], errors="coerce")
                    location = parts[1]
                    if pd.notna(date_value) and location:
                        date_location_records.append({
                            "date": date_value,
                            "location": location,
                            "id": parts[0],
                        })
        except Exception:
            date_location_records = []

    if date_location_records:
        df = pd.DataFrame(date_location_records)
        fig = go.Figure(data=[
            go.Scatter(
                x=df["date"],
                y=df["location"],
                mode="markers",
                marker=dict(size=8, opacity=0.7),
                customdata=df.get("id"),
                hovertemplate="Date: %{x}<br>Location: %{y}<br>ID: %{customdata}<extra></extra>",
            )
        ])
        # Ensure all unique locations are shown as y-axis ticks
        unique_locations = sorted(df["location"].unique())
        # Dynamically set plot height based on number of locations (min 400px)
        plot_height = max(400, 30 * len(unique_locations))
        fig.update_layout(
            title="Sampling dates by location",
            xaxis_title="Date",
            yaxis_title="Location",
            height=plot_height,
        )
        fig.update_yaxes(categoryorder="array", categoryarray=unique_locations)
        _apply_plot_style(fig)
        dataset_plot_html = _plot_div(fig)

    filter_summary = []
    if min_length is not None:
        filter_summary.append(f"Minimum length: {min_length}")
    if max_n_content is not None:
        filter_summary.append(f"Maximum N content: {max_n_content}")
    filter_summary_text = ", ".join(filter_summary) if filter_summary else "No filters applied."

    filtered_lookup = {(row["file"], row["id"]) for row in filtered_rows}

    total_sequences = sum(row["sequences"] for row in records_summary)
    filtered_count = len(filtered_rows)

    length_plot_html = ""
    lengths = []
    length_ids = []
    for path in input_fastas:
        for rec in SeqIO.parse(path, "fasta"):
            lengths.append(len(rec.seq))
            length_ids.append(rec.id)
    if lengths:
        min_len = min(lengths)
        max_len = max(lengths)
        median_len = float(np.median(lengths))
        min_ids = [seq_id for seq_id, seq_len in zip(length_ids, lengths) if seq_len == min_len]
        max_ids = [seq_id for seq_id, seq_len in zip(length_ids, lengths) if seq_len == max_len]
        span = max_len - min_len
        pad = max(50, int(span * 0.03)) if span else 50
        y_min = max(0, min_len - pad)
        y_max = max_len + pad

        def _format_length(value: float) -> str:
            if float(value).is_integer():
                return str(int(value))
            return f"{value:.1f}"

        min_ids_text = ", ".join(min_ids) if min_ids else "n/a"
        max_ids_text = ", ".join(max_ids) if max_ids else "n/a"

        fig = go.Figure(data=[
            go.Box(
                x=lengths,
                y=[0] * len(lengths),
                boxpoints="all",
                jitter=0.3,
                pointpos=0,
                marker=dict(color="#4BA3A8", opacity=0.6),
                line=dict(color="#4BA3A8"),
                hoverinfo="skip",
                name="Lengths",
                orientation="h",
            )
        ])
        fig.add_trace(
            go.Scatter(
                x=[min_len],
                y=[0],
                mode="markers",
                marker=dict(color="#c77c8a", size=8, symbol="circle"),
                hovertemplate=f"Min length: {_format_length(min_len)} bp<br>IDs: {min_ids_text}<extra></extra>",
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[max_len],
                y=[0],
                mode="markers",
                marker=dict(color="#7A6BB1", size=8, symbol="circle"),
                hovertemplate=f"Max length: {_format_length(max_len)} bp<br>IDs: {max_ids_text}<extra></extra>",
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[median_len],
                y=[0.0],
                mode="markers",
                marker=dict(color="#4BA3A8", size=10, symbol="circle"),
                hovertemplate=f"Median length: {_format_length(median_len)} bp<extra></extra>",
                showlegend=False,
            )
        )
        fig.update_layout(
            xaxis_title="Sequence length (bp)",
            yaxis_title="",
            showlegend=False,
            dragmode=False,
            hovermode="closest",
        )
        fig.update_xaxes(range=[y_min, y_max])
        fig.update_yaxes(range=[-1, 1], showticklabels=False, showgrid=False, zeroline=False)
        if min_length is not None:
            fig.add_trace(
                go.Scatter(
                    x=[min_length, min_length],
                    y=[-0.9, 0.9],
                    mode="lines",
                    line=dict(color="#c77c8a", width=2, dash="dash"),
                    hovertemplate=f"Minimum length filter: {_format_length(min_length)} bp<extra></extra>",
                    showlegend=False,
                )
            )
        _apply_plot_style(fig)
        length_plot_html = _plot_div(fig)

    provenance_inputs = ", ".join([os.path.basename(p) for p in input_fastas])
    provenance_metadata = ", ".join([os.path.basename(p) for p in metadata_paths]) if metadata_paths else "None"
    provenance_output = os.path.basename(output_fasta) if output_fasta else "stdout"
    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from raccoon import __version__ as raccoon_version
    except Exception:
        raccoon_version = "unknown"

    final_rows = []
    if output_fasta:
        try:
            for rec in SeqIO.parse(output_fasta, "fasta"):
                seq = str(rec.seq)
                final_rows.append({
                    "id": rec.id,
                    "length": len(seq),
                    "n_content": round(_n_content(seq), 4),
                })
        except Exception:
            final_rows = []

    cmd_parts = ["raccoon", "seq-qc"]
    cmd_parts.extend([os.path.basename(p) for p in input_fastas])
    if output_fasta:
        cmd_parts.extend(["--output", os.path.basename(output_fasta)])
    if metadata_paths:
        cmd_parts.append("--metadata")
        cmd_parts.extend([os.path.basename(p) for p in metadata_paths])
    if min_length is not None:
        cmd_parts.extend(["--min-length", str(min_length)])
    if max_n_content is not None:
        cmd_parts.extend(["--max-n-content", str(max_n_content)])
    cmd_line = " ".join(cmd_parts)

    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from raccoon import __version__ as raccoon_version
    except Exception:
        raccoon_version = "unknown"
    outpath = os.path.join(outdir, "seq-qc_report.html")
    filter_failures_table = _table_context(pd.DataFrame(filter_failures)) if filter_failures else None
    metadata_issues_table = _table_context(pd.DataFrame(metadata_issues)) if metadata_issues else None

    context = {
        "summary": {
            "total_sequences": total_sequences,
            "filtered_count": filtered_count,
            "locations": header_stats.get("locations", 0),
            "date_range": header_stats.get("dates", ""),
            "cmd_line": cmd_line,
            "filters": filter_summary_text,
        },
        "records_summary": records_summary,
        "seq_details_by_file": [
            {"file": fname, "details": details}
            for fname, details in seq_details_by_file.items()
        ],
        "filtered_keys": {f"{f}::{i}" for f, i in filtered_lookup},
        "filter_failures_table": filter_failures_table,
        "metadata_issues_table": metadata_issues_table,
        "length_plot_html": length_plot_html,
        "metadata_summary": metadata_summary,
        "metadata_tables": metadata_tables,
        "dataset_plot_html": dataset_plot_html,
        "final_rows": final_rows,
        "datafiles": {
            "inputs": provenance_inputs,
            "metadata": provenance_metadata,
            "output": provenance_output,
        },
        "report_metadata": {
            "generated_stamp": generated_stamp,
            "raccoon_version": raccoon_version,
            "python_version": sys.version.split()[0],
            "platform": f"{platform.system()} {platform.release()}",
        },
        "generated_stamp": generated_stamp,
    }
    _write_html(outpath, "Raccoon seq-qc report", "seq_qc.html", context)
    return outpath


def generate_alignment_report(outdir: str, alignment_path: str, mask_file: Optional[str] = None) -> str:
    lengths = []
    n_contents = []
    completeness = []
    seq_ids = []
    seq_strings = []
    for rec in SeqIO.parse(alignment_path, "fasta"):
        seq = str(rec.seq)
        lengths.append(len(seq))
        n_contents.append(_n_content(seq))
        seq_ids.append(rec.id)
        seq_strings.append(seq)
        if seq:
            valid = sum(1 for c in seq.upper() if c not in ["N", "-"])
            completeness.append(valid / len(seq))

    aln_len = _safe_max(lengths)
    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from raccoon import __version__ as raccoon_version
    except Exception:
        raccoon_version = "unknown"
    site_to_ids: Dict[int, List[str]] = {}
    if mask_file and os.path.exists(mask_file):
        try:
            with open(mask_file, "r") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    row_type = (row.get("type") or "site").strip().lower()
                    if row_type != "site":
                        continue
                    site = row.get("flagged") or row.get("Name") or row.get("site")
                    present_in = row.get("present_in", "")
                    if site is None:
                        continue
                    try:
                        site_int = int(site)
                    except Exception:
                        continue
                    if site_int < 1 or (aln_len and site_int > aln_len):
                        continue
                    ids = [v.strip() for v in str(present_in).split(",") if v.strip()]
                    if ids:
                        site_to_ids[site_int] = ids
        except Exception:
            site_to_ids = {}
    sites_table: Optional[Dict[str, Any]] = None
    sequence_removals_table: Optional[Dict[str, Any]] = None
    if mask_file and os.path.exists(mask_file):
        site_rows = []
        sequence_rows = []
        with open(mask_file, "r") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                row_type = (row.get("type") or "site").strip().lower()
                if row_type == "site":
                    site_rows.append(row)
                elif row_type == "sequence_record":
                    flagged = row.get("flagged", "")
                    note = row.get("note", "")
                    sites = [s for s in note.split(",") if s.strip()]
                    sequence_rows.append({
                        "sequence": flagged,
                        "problematic_sites": len(sites),
                        "sites": note,
                    })
        if site_rows:
            headers = list(site_rows[0].keys())
            rows = [[row.get(h, "") for h in headers] for row in site_rows]
            sites_table = {"headers": headers, "rows": rows}
        if sequence_rows:
            sequence_removals_table = _table_context(pd.DataFrame(sequence_rows))
    n_blocks_plot_html = ""
    if seq_strings and aln_len:
        z = []
        text = []
        y_positions = list(range(len(seq_ids)))
        for seq in seq_strings:
            row = [1 if c.upper() == "N" else 0 for c in seq]
            if len(row) < aln_len:
                row.extend([0] * (aln_len - len(row)))
            z.append(row[:aln_len])
            text.append(list(seq[:aln_len]))
        fig = go.Figure(data=[go.Heatmap(
            z=z,
            x=list(range(1, aln_len + 1)),
            y=y_positions,
            colorscale=[[0, "#ffffff"], [1, "#bfc3c8"]]
            # showscale=False,
            # text=text,
            # texttemplate="",
            # textfont=dict(size=8),
            # hovertemplate="ID: %{customdata}<br>Position: %{x}<br>Base: %{text}<br>N: %{z}<extra></extra>",
            # customdata=[[seq_ids[i]] * aln_len for i in range(len(seq_ids))],
        )])
        # shapes = []
        # for i in range(len(seq_ids)):
        #     if i % 2 == 1:
        #         shapes.append(dict(
        #             type="rect",
        #             xref="x",
        #             yref="y",
        #             x0=0.5,
        #             x1=aln_len + 0.5,
        #             y0=i - 0.5,
        #             y1=i + 0.5,
        #             fillcolor="#ede8f3",
        #             opacity=0.4,
        #             line_width=0,
        #             layer="below",
        #         ))
        height = max(400, len(seq_ids) * 14)
        tick_size = 12 if len(seq_ids) <= 40 else 8
        # fig.update_layout(
        #     xaxis_title="Position (bp)",
        #     yaxis_title="Sequence",
        #     yaxis=dict(tickmode="array", tickvals=y_positions, ticktext=seq_ids, tickfont=dict(size=tick_size)),
        #     showlegend=False,
        #     # shapes=shapes,
        #     height=height
        # )
        _apply_plot_style(fig)
        n_blocks_plot_html = _plot_div(fig, div_id="n-blocks-plot")

    flagged_plot_html = ""
    if mask_file and os.path.exists(mask_file):
        site_rows = []
        with open(mask_file, "r") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                row_type = (row.get("type") or "site").strip().lower()
                if row_type != "site":
                    continue
                notes = row.get("note", "").split(";")
                site_val = row.get("flagged") or row.get("Name")
                try:
                    site_int = int(site_val)
                except Exception:
                    continue
                for note in notes:
                    if note:
                        site_rows.append({"site": site_int, "note": note})
        if site_rows:
            df = pd.DataFrame(site_rows)
            fig = go.Figure()
            def _pretty_label(value: str) -> str:
                return value.replace("_", " ").title()

            for note, subset in df.groupby("note"):
                pretty = _pretty_label(note)
                fig.add_trace(go.Scatter(
                    x=subset["site"],
                    y=[pretty] * len(subset),
                    mode="markers",
                    name=pretty,
                    marker=dict(size=8),
                    customdata=[
                        "many" if len(site_to_ids.get(int(site), [])) > 5
                        else "<br>".join(site_to_ids.get(int(site), []))
                        for site in subset["site"]
                    ],
                    hovertemplate="Site: %{x}<br>Category: %{y}<br>IDs: %{customdata}<extra></extra>",
                ))
            fig.update_layout(xaxis_title="Position (bp)", yaxis_title="Category")
            fig.update_xaxes(range=[0, aln_len])
            _apply_plot_style(fig)
            flagged_plot_html = _plot_div(fig)

    diversity_plot_html = ""
    if seq_strings and aln_len:
        diversities = []
        for pos in range(aln_len):
            counts = {}
            total = 0
            for seq in seq_strings:
                if pos >= len(seq):
                    continue
                base = seq[pos].upper()
                if base in ["-", "N"]:
                    continue
                counts[base] = counts.get(base, 0) + 1
                total += 1
            if total == 0:
                diversities.append(0.0)
            else:
                h = 0.0
                for c in counts.values():
                    p = c / total
                    h -= p * math.log2(p)
                diversities.append(h)
        fig = go.Figure(data=[go.Scatter(x=list(range(1, aln_len + 1)), y=diversities, mode="lines")])
        fig.update_layout(xaxis_title="Position (bp)", yaxis_title="Shannon diversity", showlegend=False)
        fig.update_yaxes(range=[0, max(diversities) if diversities else 0])
        _apply_plot_style(fig)
        diversity_plot_html = _plot_div(fig)

    outpath = os.path.join(outdir, "aln-qc_report.html")
    context = {
        "summary": {
            "sequences": len(lengths),
            "alignment_length": aln_len,
            "mean_n_content": round(_safe_mean(n_contents), 4),
            "mean_completeness": round(_safe_mean(completeness), 4),
        },
        # "n_blocks_plot_html": n_blocks_plot_html,
        # "has_n_blocks_plot": bool(n_blocks_plot_html),
        "sites_table": sites_table,
        "sequence_removals_table": sequence_removals_table,
        "flagged_plot_html": flagged_plot_html,
        "diversity_plot_html": diversity_plot_html,
        "datafiles": {
            "alignment": os.path.basename(alignment_path),
            "mask_file": os.path.basename(mask_file) if mask_file else "None",
            "outdir": os.path.basename(outdir) if outdir else "",
        },
        "report_metadata": {
            "generated_stamp": generated_stamp,
            "raccoon_version": raccoon_version,
            "python_version": sys.version.split()[0],
            "platform": f"{platform.system()} {platform.release()}",
        },
        "generated_stamp": generated_stamp,
    }
    _write_html(outpath, "Raccoon aln-qc report", "aln_qc.html", context)
    return outpath


def generate_mask_report(
    outdir: str,
    alignment_path: str,
    mask_file: Optional[str] = None,
    output_alignment: Optional[str] = None,
) -> str:
    lengths = []
    seq_ids = []
    for rec in SeqIO.parse(alignment_path, "fasta"):
        seq = str(rec.seq)
        lengths.append(len(seq))
        seq_ids.append(rec.id)

    aln_len = _safe_max(lengths)
    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from raccoon import __version__ as raccoon_version
    except Exception:
        raccoon_version = "unknown"

    positions: List[int] = []
    sequences_to_remove: List[str] = []
    if mask_file and os.path.exists(mask_file):
        try:
            from raccoon.utils import alignment_functions as af
            positions, sequences_to_remove = af.parse_mask_rows(mask_file)
        except Exception:
            positions = []
            sequences_to_remove = []

    valid_positions = [pos for pos in positions if 1 <= pos <= aln_len]
    masked_count = len(valid_positions)
    masked_pct = (masked_count / aln_len) if aln_len else 0.0

    masked_table = None
    if seq_ids:
        rows = []
        for seq_id in seq_ids:
            rows.append({
                "sequence": seq_id,
                "masked_sites": masked_count,
                "masked_pct": round(masked_pct, 4),
            })
        masked_table = _table_context(pd.DataFrame(rows))

    mask_sites_table: Optional[Dict[str, Any]] = None
    if mask_file and os.path.exists(mask_file):
        site_rows = []
        with open(mask_file, "r") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                site_rows.append(row)
        if site_rows:
            headers = list(site_rows[0].keys())
            rows = [[row.get(h, "") for h in headers] for row in site_rows]
            mask_sites_table = {"headers": headers, "rows": rows}

    outpath = os.path.join(outdir, "mask_report.html")
    context = {
        "summary": {
            "sequences": len(seq_ids),
            "alignment_length": aln_len,
            "masked_sites": masked_count,
            "masked_pct": round(masked_pct, 4),
            "sequences_removed": len(sequences_to_remove),
        },
        "masked_table": masked_table,
        "mask_sites_table": mask_sites_table,
        "datafiles": {
            "alignment": os.path.basename(alignment_path),
            "mask_file": os.path.basename(mask_file) if mask_file else "None",
            "output_alignment": os.path.basename(output_alignment) if output_alignment else "None",
            "outdir": os.path.basename(outdir) if outdir else "",
        },
        "report_metadata": {
            "generated_stamp": generated_stamp,
            "raccoon_version": raccoon_version,
            "python_version": sys.version.split()[0],
            "platform": f"{platform.system()} {platform.release()}",
        },
        "generated_stamp": generated_stamp,
    }
    _write_html(outpath, "Raccoon mask report", "mask.html", context)
    return outpath


def generate_phylo_report(outdir: str, treefile: str, flags_csv: Optional[str] = None, tree_format: str = "auto") -> str:
    my_tree = load_tree(treefile, tree_format=tree_format)
    tip_names = []
    tip_heights = []
    tip_dates = []
    for node in my_tree.Objects:
        if node.branchType == 'leaf':
            label = ensure_node_label(node)
            if label:
                tip_names.append(label)
                tip_heights.append(node.height)
                parts = label.split("|")
                date_val = None
                if len(parts) >= 4:
                    date_val = parts[-1]
                elif len(parts) >= 3:
                    date_val = parts[-1]
                tip_dates.append(date_val)

    flags_df = None
    if flags_csv and os.path.exists(flags_csv):
        try:
            flags_df = pd.read_csv(flags_csv)
        except Exception:
            flags_df = None

    convergent_table = None
    reversion_table = None
    immune_editing_table = None
    if flags_df is not None and not flags_df.empty and "present_in" in flags_df.columns:
        def _merge_present_in(frame: pd.DataFrame) -> pd.DataFrame:
            if frame.empty:
                return frame
            cols = [c for c in ["site", "mutation_type"] if c in frame.columns]
            def _join_unique(values: pd.Series) -> str:
                items = [str(v) for v in values if pd.notna(v) and str(v).strip()]
                return ";".join(sorted(set(items)))
            aggregations = {
                "present_in": _join_unique,
            }
            if "mask_boolean" in frame.columns:
                aggregations["mask_boolean"] = "any"
            merged = frame.groupby(cols, dropna=False, as_index=False).agg(aggregations)
            return merged

        convergent = flags_df[flags_df["mutation_type"].str.contains("convergent", case=False, na=False)]
        convergent = _merge_present_in(convergent)
        convergent_table = _table_context(convergent)
        reversion = flags_df[flags_df["mutation_type"].str.contains("reversion", case=False, na=False)]
        reversion = _merge_present_in(reversion)
        reversion_table = _table_context(reversion)
        immune = flags_df[flags_df["mutation_type"].str.contains("adar|apobec", case=False, na=False)]
        immune_editing_table = _table_context(immune)

    generated_stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from raccoon import __version__ as raccoon_version
    except Exception:
        raccoon_version = "unknown"

    root_to_tip_plot = "<p>No root-to-tip distances available.</p>"
    if tip_heights and tip_dates:
        date_series = pd.to_datetime(pd.Series(tip_dates), errors="coerce")
        mask = date_series.notna()
        if mask.any():
            x_dates = date_series[mask]
            y_heights = np.array(tip_heights)[mask.values]
            x_num = x_dates.map(pd.Timestamp.toordinal).astype(float).values
            if len(x_num) >= 2:
                slope, intercept = np.polyfit(x_num, y_heights, 1)
                y_hat = slope * x_num + intercept
                resid = y_heights - y_hat
                n = len(x_num)
                s_err = np.sqrt(np.sum(resid ** 2) / max(n - 2, 1))
                x_mean = np.mean(x_num)
                s_xx = np.sum((x_num - x_mean) ** 2) or 1.0
                ci = 3.0 * s_err * np.sqrt(1 / n + (x_num - x_mean) ** 2 / s_xx)
                upper = y_hat + ci
                lower = y_hat - ci
                outside = (y_heights > upper) | (y_heights < lower)
                order = np.argsort(x_num)
                x_dates = x_dates.iloc[order]
                y_heights = y_heights[order]
                y_hat = y_hat[order]
                upper = upper[order]
                lower = lower[order]

                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=x_dates,
                    y=y_heights,
                    mode="markers",
                    marker=dict(color="#4BA3A8", size=8),
                    text=[tip_names[i] for i, m in enumerate(mask.values) if m],
                    hovertemplate="%{text}<br>Date: %{x|%Y-%m-%d}<br>Distance: %{y:.4f}<extra></extra>",
                    name="Tips",
                ))
                x_min = float(x_num.min())
                x_max = float(x_num.max())
                x_line = [x_dates.min(), x_dates.max()]
                y_line = [slope * x_min + intercept, slope * x_max + intercept]
                fig.add_trace(go.Scatter(
                    x=x_line,
                    y=y_line,
                    mode="lines",
                    line=dict(color="#7A6BB1"),
                    hoverinfo="skip",
                    name="Regression",
                ))
                fig.add_trace(go.Scatter(
                    x=x_dates,
                    y=upper,
                    mode="lines",
                    line=dict(color="rgba(0,0,0,0)"),
                    hoverinfo="skip",
                    showlegend=False,
                ))
                fig.add_trace(go.Scatter(
                    x=x_dates,
                    y=lower,
                    mode="lines",
                    line=dict(color="rgba(0,0,0,0)"),
                    fill="tonexty",
                    fillcolor="rgba(182,170,201,0.25)",
                    hoverinfo="skip",
                    showlegend=False,
                ))
                fig.update_layout(
                    xaxis_title="Date",
                    yaxis_title="Root-to-tip distance",
                    yaxis_tickformat=".1e",
                    showlegend=False,
                )
                _apply_plot_style(fig)
                root_to_tip_plot = _plot_div(fig)

    mutation_types_plot = "<p>No mutation types available.</p>"
    if flags_df is not None and not flags_df.empty and "mutation_type" in flags_df:
        counts = flags_df["mutation_type"].value_counts().reset_index()
        counts.columns = ["mutation_type", "count"]
        fig = go.Figure(data=[go.Bar(x=counts["mutation_type"], y=counts["count"])])
        fig.update_layout(xaxis_title="Type", yaxis_title="Count", showlegend=False)
        _apply_plot_style(fig)
        mutation_types_plot = _plot_div(fig)

    phylogeny_base = os.path.splitext(os.path.basename(treefile))[0]
    branch_snps_path = os.path.join(outdir, f"{phylogeny_base}.branch_snps.reconstruction.csv")
    if not os.path.exists(branch_snps_path):
        branch_snps_path = None
    tree_plot_html = build_tree_plot(
        treefile,
        tree_format=tree_format,
        branch_snps_path=branch_snps_path,
    )

    outpath = os.path.join(outdir, "tree-qc_report.html")
    context = {
        "summary": {
            "tips": len(tip_names),
            "tree_height": getattr(my_tree, "treeHeight", "n/a"),
            "y_span": getattr(my_tree, "ySpan", "n/a"),
        },
        "tree_plot_html": tree_plot_html,
        "root_to_tip_plot_html": root_to_tip_plot,
        "convergent_table": convergent_table,
        "reversion_table": reversion_table,
        "immune_editing_table": immune_editing_table,
        "mutation_types_plot_html": mutation_types_plot,
        "datafiles": {
            "treefile": os.path.basename(treefile),
            "flags_csv": os.path.basename(flags_csv) if flags_csv else "None",
            "outdir": os.path.basename(outdir) if outdir else "",
        },
        "report_metadata": {
            "generated_stamp": generated_stamp,
            "raccoon_version": raccoon_version,
            "python_version": sys.version.split()[0],
            "platform": f"{platform.system()} {platform.release()}",
        },
        "generated_stamp": generated_stamp,
    }
    _write_html(outpath, "Raccoon tree-qc report", "tree_qc.html", context)
    return outpath
