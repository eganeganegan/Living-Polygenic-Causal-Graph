import pytest
from ingestion.parsers import extract_gene_relations

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_PUBLICATION_ROLE_KEYED = {
    "id": "12345678",
    "date": "2024-01-15",
    "passages": [
        {
            "annotations": [
                {"id": "T1", "infons": {"type": "Gene", "identifier": "7157", "name": "TP53"}, "text": "TP53"},
                {"id": "T2", "infons": {"type": "Gene", "identifier": "675", "name": "BRCA2"}, "text": "BRCA2"},
            ]
        }
    ],
    "relations": [
        {
            "id": "R1",
            "infons": {
                "type": "Association",
                "role1": {"type": "Gene", "identifier": "7157", "name": "TP53", "refid": "T1"},
                "role2": {"type": "Gene", "identifier": "675", "name": "BRCA2", "refid": "T2"},
            },
        }
    ],
}

_PUBLICATION_NODES_FALLBACK = {
    "id": "99999999",
    "date": "2024-02-01",
    "passages": [
        {
            "annotations": [
                {"id": "T1", "infons": {"type": "Gene", "identifier": "7157", "name": "TP53"}, "text": "TP53"},
                {"id": "T2", "infons": {"type": "Gene", "identifier": "675", "name": "BRCA2"}, "text": "BRCA2"},
            ]
        }
    ],
    "relations": [
        {
            "id": "R1",
            "infons": {"type": "Binding"},
            "nodes": [{"refid": "T1", "role": "subject"}, {"refid": "T2", "role": "object"}],
        }
    ],
}

# ---------------------------------------------------------------------------
# Parser tests
# ---------------------------------------------------------------------------

def test_extract_role_keyed_returns_correct_fields():
    records = extract_gene_relations(_PUBLICATION_ROLE_KEYED)
    assert len(records) == 1
    r = records[0]
    assert r["source_gene_id"] == "7157"
    assert r["target_gene_id"] == "675"
    assert r["source_gene_name"] == "TP53"
    assert r["target_gene_name"] == "BRCA2"
    assert r["relationship_type"] == "Association"
    assert r["pmid"] == "12345678"


def test_extract_nodes_fallback_path():
    """Nodes-keyed format (no role1/role2) should parse correctly."""
    records = extract_gene_relations(_PUBLICATION_NODES_FALLBACK)
    assert len(records) == 1
    assert records[0]["source_gene_id"] == "7157"
    assert records[0]["relationship_type"] == "Binding"


def test_extract_skips_self_loop():
    pub = {
        "id": "1", "passages": [],
        "relations": [{"id": "R1", "infons": {
            "type": "Association",
            "role1": {"type": "Gene", "identifier": "7157"},
            "role2": {"type": "Gene", "identifier": "7157"},
        }}],
    }
    assert extract_gene_relations(pub) == []


def test_extract_skips_non_gene_role():
    pub = {
        "id": "1", "passages": [],
        "relations": [{"id": "R1", "infons": {
            "type": "Association",
            "role1": {"type": "Gene", "identifier": "7157"},
            "role2": {"type": "Disease", "identifier": "MESH:D012559"},
        }}],
    }
    assert extract_gene_relations(pub) == []


def test_extract_skips_missing_identifier():
    pub = {
        "id": "1", "passages": [],
        "relations": [{"id": "R1", "infons": {
            "type": "Association",
            "role1": {"type": "Gene", "identifier": ""},
            "role2": {"type": "Gene", "identifier": "675"},
        }}],
    }
    assert extract_gene_relations(pub) == []


def test_extract_empty_publication():
    assert extract_gene_relations({"id": "1", "passages": [], "relations": []}) == []
