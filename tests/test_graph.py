from unittest.mock import MagicMock, patch

import pytest

from graph.queries import _normalise, fetch_available_diseases, fetch_disease_subgraph

# ---------------------------------------------------------------------------
# _normalise
# ---------------------------------------------------------------------------

_BASE_RECORD = {
    "source_gene_id": "100",
    "target_gene_id": "200",
    "source_gene_name": "ALPHA",
    "target_gene_name": "BETA",
    "relationship_type": "Association",
    "statistical_weight": 0.5,
    "publication_date": "2024-01-01",
    "gwas_significant": False,
    "disease_tag": "Schizophrenia",
}


def test_normalise_lower_id_is_source():
    r = _normalise({**_BASE_RECORD, "source_gene_id": "999", "target_gene_id": "1"})
    assert r["src_id"] == "1"
    assert r["tgt_id"] == "999"


def test_normalise_names_follow_flipped_ids():
    r = _normalise({**_BASE_RECORD, "source_gene_id": "999", "target_gene_id": "1",
                    "source_gene_name": "HIGH", "target_gene_name": "LOW"})
    assert r["src_name"] == "LOW"
    assert r["tgt_name"] == "HIGH"


def test_normalise_symmetric():
    """Forward and reverse records must produce identical output."""
    fwd = _normalise(_BASE_RECORD)
    rev = _normalise({
        **_BASE_RECORD,
        "source_gene_id": _BASE_RECORD["target_gene_id"],
        "target_gene_id": _BASE_RECORD["source_gene_id"],
        "source_gene_name": _BASE_RECORD["target_gene_name"],
        "target_gene_name": _BASE_RECORD["source_gene_name"],
    })
    assert fwd == rev


def test_normalise_already_in_order():
    r = _normalise(_BASE_RECORD)  # "100" < "200"
    assert r["src_id"] == "100"
    assert r["tgt_id"] == "200"


def test_normalise_includes_disease_tag():
    r = _normalise(_BASE_RECORD)
    assert r["disease_tag"] == "Schizophrenia"


# ---------------------------------------------------------------------------
# fetch_available_diseases
# ---------------------------------------------------------------------------

def test_fetch_available_diseases_returns_list():
    mock_driver = MagicMock()
    with patch("graph.queries.run_query", return_value=[
        {"disease_tag": "Alzheimers"},
        {"disease_tag": "Schizophrenia"},
    ]):
        result = fetch_available_diseases(mock_driver)
    assert result == ["Alzheimers", "Schizophrenia"]


def test_fetch_available_diseases_empty_graph():
    mock_driver = MagicMock()
    with patch("graph.queries.run_query", return_value=[]):
        assert fetch_available_diseases(mock_driver) == []


# ---------------------------------------------------------------------------
# fetch_disease_subgraph
# ---------------------------------------------------------------------------

def test_fetch_disease_subgraph_passes_disease_tag():
    mock_driver = MagicMock()
    with patch("graph.queries.run_query", return_value=[]) as mock_rq:
        fetch_disease_subgraph(mock_driver, "Schizophrenia")
    _, call_kwargs = mock_rq.call_args
    assert call_kwargs.get("disease_tag") == "Schizophrenia"
