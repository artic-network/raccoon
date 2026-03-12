"""Tests for raccoon.commands.phylo helper functions."""
from raccoon.commands import phylo as phylo_cmd
from raccoon.utils.reporting import _extract_date_from_tip_label


# ---------------------------------------------------------------------------
# _parse_tip_field_names
# ---------------------------------------------------------------------------

def test_parse_tip_field_names_plain_delimited():
    assert phylo_cmd._parse_tip_field_names("sample|location|date") == ["sample", "location", "date"]


def test_parse_tip_field_names_brace_wrapped():
    assert phylo_cmd._parse_tip_field_names("{id}|{location}|{date}") == ["id", "location", "date"]


def test_parse_tip_field_names_single_field():
    assert phylo_cmd._parse_tip_field_names("{id}") == ["id"]


def test_parse_tip_field_names_empty_returns_empty():
    assert phylo_cmd._parse_tip_field_names("") == []
    assert phylo_cmd._parse_tip_field_names(None) == []


# ---------------------------------------------------------------------------
# _extract_date_from_tip_label  — cascade logic
# ---------------------------------------------------------------------------

def test_extract_date_parses_named_date_field():
    """Named template date field used when present and parseable."""
    raw, dt, prec = _extract_date_from_tip_label(
        "sampleA|London|2024-03-15",
        tip_fields="{id}|{location}|{date}",
    )
    assert raw == "2024-03-15"
    assert dt is not None
    assert dt.year == 2024


def test_extract_date_falls_back_to_last_field_when_named_field_missing():
    """No 'date' named field in template — last field tried."""
    raw, dt, prec = _extract_date_from_tip_label(
        "sampleA|2024-06-01",
        tip_fields="{id}|{location}",
    )
    # named 'date' field absent; last field "2024-06-01" should parse
    assert raw == "2024-06-01"
    assert dt is not None
    assert dt.year == 2024


def test_extract_date_falls_back_to_last_field_when_named_field_out_of_range():
    """Label is short — named date field index has no corresponding part; fall back to last."""
    raw, dt, prec = _extract_date_from_tip_label(
        "sampleA|2024-01-01",
        tip_fields="{id}|{location}|{date}",
    )
    # date field is index 2 but label only has 2 parts; last field (index 1) is "2024-01-01"
    assert raw == "2024-01-01"
    assert dt is not None


def test_extract_date_single_field_id_no_template_returns_none():
    """Plain ID label with no parseable date anywhere → (None, None, 'day')."""
    raw, dt, prec = _extract_date_from_tip_label("sampleID123")
    assert raw is None
    assert dt is None


def test_extract_date_single_field_id_with_default_template_returns_none():
    """Single-field label with the default 3-field template → graceful None."""
    raw, dt, prec = _extract_date_from_tip_label(
        "sampleID123",
        tip_fields="{id}|{location}|{date}",
    )
    assert raw is None
    assert dt is None


def test_extract_date_partial_tree_some_tips_have_dates():
    """Confirms mixed-date trees return None only for undateable tips."""
    cases = [
        ("sampleA|Loc|2024-01-01", True),
        ("sampleB", False),
        ("sampleC|Loc|2023-12-31", True),
    ]
    for label, expect_date in cases:
        _, dt, _ = _extract_date_from_tip_label(label, tip_fields="{id}|{location}|{date}")
        if expect_date:
            assert dt is not None, f"Expected date for {label}"
        else:
            assert dt is None, f"Expected no date for {label}"


def test_extract_date_year_only_precision():
    _, dt, prec = _extract_date_from_tip_label("seq|Loc|2024", tip_fields="{id}|{location}|{date}")
    assert dt is not None
    assert prec == "year"


def test_extract_date_month_precision():
    _, dt, prec = _extract_date_from_tip_label("seq|Loc|2024-06", tip_fields="{id}|{location}|{date}")
    assert dt is not None
    assert prec == "month"


def test_extract_date_missing_middle_field_uses_last_field_fallback():
    """Missing middle field should still allow date extraction from final field."""
    raw, dt, _ = _extract_date_from_tip_label(
        "seq||2024-06-15",
        tip_fields="{id}|{location}|{date}",
    )
    assert raw == "2024-06-15"
    assert dt is not None


def test_extract_date_incorrect_format_returns_none():
    raw, dt, _ = _extract_date_from_tip_label(
        "seq|Loc|not-a-date",
        tip_fields="{id}|{location}|{date}",
    )
    assert raw is None
    assert dt is None


def test_extract_date_named_date_invalid_but_last_field_valid_falls_back():
    """If configured date field is invalid, parser should attempt the final label field."""
    raw, dt, _ = _extract_date_from_tip_label(
        "seq|not-a-date|2025-01",
        tip_fields="{sample}|{date}|{other}",
    )
    assert raw == "2025-01"
    assert dt is not None
