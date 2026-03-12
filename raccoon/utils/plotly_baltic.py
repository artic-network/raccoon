"""Plotly-based interactive tree rendering utilities.
This is used for tree visualization in the phylo-qc report and may be useful for other applications.
The code was inspired by the baltic library by Gytis Dudas (https://github.com/evogytis/baltic)
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional
from datetime import datetime

import csv
import os
import re

import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd

from .reconstruction_functions import load_tree, ensure_node_label


def _scale_values(values: List[Optional[float]], factor: float) -> List[Optional[float]]:
    return [None if v is None else v * factor for v in values]


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
        config={"displayModeBar": False, "showTips": False, "doubleClick": False, "responsive": True},
        auto_play=False,
        div_id=div_id,
    )


def _format_traits(traits: Dict[str, Any]) -> str:
    lines = []
    for key in sorted(traits):
        value = traits.get(key)
        if value is None or value == "":
            continue
        lines.append(f"{key}: {value}")
    return "<br>".join(lines)


def _read_branch_snps(branch_snps_path: Optional[str]) -> Dict[str, List[Dict[str, str]]]:
    if not branch_snps_path or not os.path.exists(branch_snps_path):
        return {}
    branch_snps: Dict[str, List[Dict[str, str]]] = {}
    with open(branch_snps_path, "r") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parent = row.get("parent", "")
            child = row.get("child", "")
            if not parent or not child:
                continue
            branch_key = f"{parent}_{child}"
            branch_snps.setdefault(branch_key, []).append({
                "site": row.get("site", ""),
                "snp": row.get("snp", ""),
                "dimer": row.get("dimer", ""),
            })
    return branch_snps


def _format_branch_snps(branch_key: str, snps: List[Dict[str, str]], max_rows: int = 20) -> str:
    if not snps:
        return f"Branch: {branch_key}<br>No mutations"
    lines = [f"Branch: {branch_key}", f"Mutations: {len(snps)}"]
    for item in snps[:max_rows]:
        site = item.get("site", "")
        snp = item.get("snp", "")
        dimer = item.get("dimer", "")
        suffix = f" ({dimer})" if dimer else ""
        lines.append(f"{site}: {snp}{suffix}")
    if len(snps) > max_rows:
        lines.append(f"… and {len(snps) - max_rows} more")
    return "<br>".join(lines)


def _parse_tip_field_names(tip_fields: Optional[str]) -> List[str]:
    if not tip_fields:
        return []
    parsed: List[str] = []
    for field in str(tip_fields).split("|"):
        name = field.strip()
        if not name:
            continue
        if name.startswith("{") and name.endswith("}") and len(name) > 2:
            name = name[1:-1].strip()
        if name:
            parsed.append(name)
    return parsed


def _parse_label_traits(label: str, tip_field_names: Optional[List[str]] = None) -> Dict[str, Any]:
    traits: Dict[str, Any] = {}
    if not label:
        return traits
    parts = [p.strip() for p in label.split("|")]
    date_field_index: Optional[int] = None

    include_index_traits = not (tip_field_names or [])
    if include_index_traits:
        for idx, value in enumerate(parts):
            if value:
                traits.setdefault(f"field_{idx}", value)

        if parts and parts[-1]:
            traits.setdefault("field_-1", parts[-1])

    for idx, field_name in enumerate(tip_field_names or []):
        if idx < len(parts) and parts[idx]:
            traits.setdefault(field_name, parts[idx])
        if field_name.lower() in {"date", "dates"}:
            date_field_index = idx

    date_candidates: List[str] = []
    if date_field_index is not None and 0 <= date_field_index < len(parts):
        candidate = parts[date_field_index].strip()
        if candidate:
            date_candidates.append(candidate)

    # Preserve legacy behaviour of trying the final field as date-like when:
    # - no template fields were supplied, OR
    # - labels contain multiple fields and named-date parsing did not succeed.
    # This restores year/date derivation for labels such as "id|loc|cluster|2024"
    # while still avoiding synthetic date/year traits for single-field IDs.
    if parts and (include_index_traits or len(parts) > 1):
        trailing = parts[-1].strip()
        if trailing and trailing not in date_candidates:
            date_candidates.append(trailing)

    parsed_date = None
    parsed_raw = ""
    for raw in date_candidates:
        try:
            parsed = pd.to_datetime(raw, errors="coerce")
            if pd.notna(parsed):
                parsed_date = parsed.to_pydatetime()
                parsed_raw = raw
                break
        except Exception:
            continue

    if parsed_date is not None:
        if re.match(r"^\d{4}$", parsed_raw):
            traits.setdefault("date", f"{parsed_date.year:04d}")
        elif re.match(r"^\d{4}[-/]\d{2}$", parsed_raw):
            traits.setdefault("date", f"{parsed_date.year:04d}-{parsed_date.month:02d}")
        else:
            traits.setdefault("date", parsed_date.date().isoformat())
        traits.setdefault("year", str(parsed_date.year))
    return traits


def build_tree_plot(
    treefile: str,
    tree_format: str = "auto",
    branch_snps_path: Optional[str] = None,
    tip_fields: str = "sample|location|date",
) -> str:
    try:
        tree = load_tree(treefile, tree_format=tree_format)
    except Exception:
        return ""

    try:
        tree.drawTree()
    except Exception:
        return ""

    line_x: List[Optional[float]] = []
    line_y: List[Optional[float]] = []

    try:
        internal_nodes = tree.getInternal()
    except Exception:
        internal_nodes = []

    for node in internal_nodes:
        if not getattr(node, "children", None):
            continue
        child_ys = []
        for child in node.children:
            if child.x is None or child.y is None or node.x is None:
                continue
            line_x.extend([node.x, child.x, None])
            line_y.extend([child.y, child.y, None])
            child_ys.append(child.y)
        if child_ys and node.x is not None:
            line_x.extend([node.x, node.x, None])
            line_y.extend([min(child_ys), max(child_ys), None])

    tip_x: List[float] = []
    tip_y: List[float] = []
    tip_hover: List[str] = []
    tip_traits: List[Dict[str, Any]] = []
    tip_labels: List[str] = []

    node_x: List[float] = []
    node_y: List[float] = []
    node_hover: List[str] = []

    branch_hover_x: List[float] = []
    branch_hover_y: List[float] = []
    branch_hover_text: List[str] = []

    branch_snps = _read_branch_snps(branch_snps_path)
    tip_field_names = _parse_tip_field_names(tip_fields)

    for obj in getattr(tree, "Objects", []):
        label = ensure_node_label(obj) or getattr(obj, "name", "") or ""
        traits = dict(getattr(obj, "traits", {}) or {})
        traits.setdefault("label", label)
        traits.setdefault("branch_type", getattr(obj, "branchType", ""))
        traits.update(_parse_label_traits(label, tip_field_names=tip_field_names))

        if getattr(obj, "is_leaflike", lambda: False)():
            if obj.x is None or obj.y is None:
                continue
            tip_x.append(obj.x)
            tip_y.append(obj.y)
            tip_traits.append(traits)
            tip_hover.append(_format_traits(traits))
            tip_labels.append(label)
        elif getattr(obj, "is_node", lambda: False)():
            if obj.x is None or obj.y is None:
                continue
            node_x.append(obj.x)
            node_y.append(obj.y)
            node_hover.append(_format_traits(traits))

    for node in internal_nodes:
        parent_label = ensure_node_label(node) or getattr(node, "name", "") or ""
        for child in getattr(node, "children", []) or []:
            if child.x is None or child.y is None or node.x is None:
                continue
            child_label = ensure_node_label(child) or getattr(child, "name", "") or ""
            branch_key = f"{parent_label}_{child_label}"
            snps = branch_snps.get(branch_key)
            if not snps:
                continue
            branch_hover_x.append((node.x + child.x) / 2)
            branch_hover_y.append(child.y)
            branch_hover_text.append(_format_branch_snps(branch_key, snps))

    color_keys = set()
    for traits in tip_traits:
        for key, value in traits.items():
            if value is not None and value != "":
                color_keys.add(key)
    excluded_color_keys = {"label", "branch_type"}
    preferred = [key for key in ["date", "year", *tip_field_names] if key in color_keys]
    if not tip_field_names:
        preferred.extend([key for key in ["field_1", "field_-1"] if key in color_keys])
    color_keys = preferred + sorted([k for k in color_keys if k not in preferred and k not in excluded_color_keys])
    default_key = color_keys[0] if color_keys else "__none__"

    palette = [
        "#4BA3A8",
        "#7A6BB1",
        "#D08BA8",
        "#5A9BD4",
        "#B9A44C",
        "#8BC59A",
        "#E26D5C",
        "#6C8EA5",
    ]

    color_maps: Dict[str, Dict[str, str]] = {}
    for key in color_keys:
        values = sorted({str(t.get(key, "")) for t in tip_traits if t.get(key, "") != ""})
        mapping = {}
        for idx, value in enumerate(values):
            mapping[value] = palette[idx % len(palette)]
        color_maps[key] = mapping

    def _colors_for(key: str) -> List[str]:
        mapping = color_maps.get(key, {})
        return [mapping.get(str(t.get(key, "")), "#4d4d4d") for t in tip_traits]

    x_values = [v for v in line_x if v is not None]
    x_values.extend(tip_x)
    x_values.extend(node_x)
    x_min = min(x_values) if x_values else 0.0
    x_max = max(x_values) if x_values else 1.0
    x_pad = (x_max - x_min) * 0.02 if x_values else 0.1
    x_min_pad = x_min - x_pad
    x_max_pad = x_max + x_pad

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=line_x,
        y=line_y,
        mode="lines",
        line=dict(color="#222", width=1),
        hoverinfo="skip",
        showlegend=False,
        name="Branches",
    ))
    fig.add_trace(go.Scatter(
        x=tip_x,
        y=tip_y,
        mode="markers",
        marker=dict(size=8, color=_colors_for(default_key)),
        hovertemplate="%{text}<extra></extra>",
        text=tip_hover,
        name="Tips",
    ))
    fig.add_trace(go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        marker=dict(size=4, color="#999999"),
        hovertemplate="%{text}<extra></extra>",
        text=node_hover,
        name="Nodes",
    ))
    if tip_x:
        tip_span = max(tip_x) - min(tip_x)
        x_offset = tip_span * 0.01 if tip_span else 0.01
        tip_label_x = [x + x_offset for x in tip_x]
    else:
        tip_label_x = []

    label_max_x = max(tip_label_x) if tip_label_x else None
    if label_max_x is not None:
        label_span = (x_max - x_min) if x_values else 1.0
        max_label_len = max((len(label) for label in tip_labels), default=0)
        label_extra = label_span * max_label_len * 0.006
        if label_max_x + label_extra > x_max_pad:
            x_max_pad = label_max_x + label_extra + x_pad
    x_range_default = (x_min_pad, x_max_pad)

    fig.add_trace(go.Scatter(
        x=tip_label_x,
        y=tip_y,
        mode="text",
        text=tip_labels,
        textposition="middle right",
        textfont=dict(size=10, color="#222"),
        hoverinfo="skip",
        showlegend=False,
        name="Tip labels",
        visible=True,
    ))
    fig.add_trace(go.Scatter(
        x=branch_hover_x,
        y=branch_hover_y,
        mode="markers",
        marker=dict(size=10, color="rgba(0,0,0,0)"),
        hovertemplate="%{text}<extra></extra>",
        text=branch_hover_text,
        showlegend=False,
        name="Branch mutations",
    ))

    scale_values = [round(1.0 + (i * 0.01), 2) for i in range(101)]
    y_min = min(tip_y) if tip_y else 0
    y_max = max(tip_y) if tip_y else 0
    y_span = y_max - y_min
    y_pad = max(1.0, y_span * 0.03) if tip_y else 1.0
    y_min_pad = y_min - y_pad
    y_max_pad = y_max + y_pad
    base_height = max(500, int((y_max_pad if tip_y else 0) * 15)) if tip_y else 500
    y_range_default = (y_min_pad, y_max_pad)
    height_default = base_height

    slider_steps = []
    for scale in scale_values:
        slider_steps.append({
            "label": "",
            "method": "relayout",
            "args": [{
                "height": int(base_height * scale),
                "yaxis.range": list(y_range_default),
            }],
        })

    color_buttons = []
    for key in color_keys:
        color_buttons.append({
            "label": key,
            "method": "restyle",
            "args": [{"marker.color": [_colors_for(key)]}, [1]],
        })

    fig.update_layout(
        xaxis_title="Substitutions per site",
        showlegend=False,
        height=height_default,
        xaxis=dict(range=list(x_range_default)),
        yaxis=dict(range=list(y_range_default) if tip_y else None),
        sliders=[{
            "active": 0,
            "pad": {"t": 0, "b": 0},
            "ticklen": 0,
            "tickwidth": 0,
            "len": 0.28,
            "x": 0.62,
            "y": 1.12,
            "xanchor": "left",
            "yanchor": "top",
            "steps": slider_steps,
        }],
        updatemenus=[{
            "buttons": color_buttons,
            "direction": "down",
            "showactive": True,
            "x": 0.12,
            "y": 1.12,
            "xanchor": "left",
            "yanchor": "top",
        }, {
            "buttons": [
                {"label": "On", "method": "restyle", "args": [{"visible": True}, [3]]},
                {"label": "Off", "method": "restyle", "args": [{"visible": False}, [3]]},
            ],
            "direction": "down",
            "showactive": False,
            "x": 0.36,
            "y": 1.12,
            "xanchor": "left",
            "yanchor": "top",
        }, {
            "buttons": [
                {
                    "label": "Reset view",
                    "method": "relayout",
                    "args": [{
                        "xaxis.range[0]": x_range_default[0],
                        "xaxis.range[1]": x_range_default[1],
                        "yaxis.autorange": False,
                        "yaxis.range[0]": y_range_default[0],
                        "yaxis.range[1]": y_range_default[1],
                        "height": height_default,
                        "sliders[0].active": 0,
                    }],
                }
            ],
            "type": "buttons",
            "direction": "right",
            "showactive": False,
            "x": 0.98,
            "y": 1.12,
            "xanchor": "right",
            "yanchor": "top",
        }],
        annotations=[
            {
                "text": "Colour by",
                "xref": "paper",
                "yref": "paper",
                "x": 0.02,
                "y": 1.12,
                "showarrow": False,
                "xanchor": "left",
                "yanchor": "top",
                "font": {"size": 12, "color": "#111"},
            },
            {
                "text": "Tip labels",
                "xref": "paper",
                "yref": "paper",
                "x": 0.28,
                "y": 1.12,
                "showarrow": False,
                "xanchor": "left",
                "yanchor": "top",
                "font": {"size": 12, "color": "#111"},
            },
            {
                "text": "Expand tree",
                "xref": "paper",
                "yref": "paper",
                "x": 0.56,
                "y": 1.12,
                "showarrow": False,
                "xanchor": "left",
                "yanchor": "top",
                "font": {"size": 12, "color": "#111"},
            },
        ],
        margin=dict(l=60, r=20, t=60, b=40),
    )
    _apply_plot_style(fig)
    fig.update_yaxes(showticklabels=False, ticks="", title=None, showline=False)
    return _plot_div(fig)
