import io
from pathlib import Path

import pandas as pd
import pytest

from transformation.bayesian_updater import (
    GWAS_SIGNIFICANCE_THRESHOLD,
    _bayesian_update,
    _gwas_significant_genes,
    update_interaction_weights,
)
from transformation.gwas_integrator import load_gwas

# ---------------------------------------------------------------------------
# _bayesian_update
# ---------------------------------------------------------------------------

def test_bayesian_update_increases_prior():
    assert _bayesian_update(0.1, 0.8, 0.1) > 0.1


def test_bayesian_update_bounded():
    # Result must stay strictly within (0, 1)
    assert 0 < _bayesian_update(0.999, 0.8, 0.1) < 1.0
    assert 0 < _bayesian_update(0.001, 0.8, 0.1) < 1.0


def test_gwas_boosted_likelihood_exceeds_base():
    prior = 0.2
    base = _bayesian_update(prior, 0.80, 0.10)
    boosted = _bayesian_update(prior, 0.92, 0.10)
    assert boosted > base


# ---------------------------------------------------------------------------
# _gwas_significant_genes
# ---------------------------------------------------------------------------

def test_gwas_significant_uses_min_per_gene():
    """A gene is significant only if its MINIMUM p-value is below threshold."""
    df = pd.DataFrame({
        "gene_id": ["7157", "7157", "675"],
        "p_value": [0.001, 1e-12, 0.5],  # 7157 min = 1e-12 (sig), 675 = 0.5 (not sig)
    })
    sig = _gwas_significant_genes(df)
    assert "7157" in sig
    assert "675" not in sig


def test_gwas_significant_empty_df():
    assert _gwas_significant_genes(pd.DataFrame(columns=["gene_id", "p_value"])) == frozenset()


# ---------------------------------------------------------------------------
# update_interaction_weights
# ---------------------------------------------------------------------------

def test_update_weights_no_gwas_uses_base_likelihood():
    interactions = [{"source_gene_id": "A", "target_gene_id": "B", "statistical_weight": 0.1}]
    update_interaction_weights(interactions, pd.DataFrame(columns=["gene_id", "p_value"]))
    assert 0.1 < interactions[0]["statistical_weight"] < 1.0
    assert interactions[0]["gwas_significant"] is False


def test_update_weights_gwas_hit_flagged_and_higher():
    base_interactions = [{"source_gene_id": "A", "target_gene_id": "B", "statistical_weight": 0.1}]
    gwas_interactions = [{"source_gene_id": "A", "target_gene_id": "B", "statistical_weight": 0.1}]
    gwas_df = pd.DataFrame({"gene_id": ["A"], "p_value": [1e-10]})

    update_interaction_weights(base_interactions, pd.DataFrame(columns=["gene_id", "p_value"]))
    update_interaction_weights(gwas_interactions, gwas_df)

    assert gwas_interactions[0]["statistical_weight"] > base_interactions[0]["statistical_weight"]
    assert gwas_interactions[0]["gwas_significant"] is True


def test_prior_clamping_prevents_degenerate_output():
    """Priors at 0 or 1 must still produce a valid posterior."""
    for extreme in [0.0, 1.0]:
        interactions = [{"source_gene_id": "X", "target_gene_id": "Y", "statistical_weight": extreme}]
        update_interaction_weights(interactions, pd.DataFrame(columns=["gene_id", "p_value"]))
        w = interactions[0]["statistical_weight"]
        assert 0 < w < 1


def test_update_returns_same_list():
    interactions = [{"source_gene_id": "A", "target_gene_id": "B"}]
    result = update_interaction_weights(interactions, pd.DataFrame(columns=["gene_id", "p_value"]))
    assert result is interactions


# ---------------------------------------------------------------------------
# load_gwas
# ---------------------------------------------------------------------------

def test_load_gwas_simple_format(tmp_path: Path):
    p = tmp_path / "test.tsv"
    p.write_text("gene_id\tp_value\n7157\t1e-10\n675\t0.05\n")
    df = load_gwas(p)
    assert list(df.columns) == ["gene_id", "p_value"]
    assert len(df) == 2
    assert df.loc[df["gene_id"] == "7157", "p_value"].iloc[0] == pytest.approx(1e-10)


def test_load_gwas_missing_file(tmp_path: Path):
    df = load_gwas(tmp_path / "nonexistent.tsv")
    assert df.empty
    assert list(df.columns) == ["gene_id", "p_value"]


def test_load_gwas_drops_invalid_pvalues(tmp_path: Path):
    p = tmp_path / "test.tsv"
    p.write_text("gene_id\tp_value\n7157\t1e-10\n675\tNA\n")
    df = load_gwas(p)
    assert len(df) == 1
    assert df.iloc[0]["gene_id"] == "7157"


def test_load_gwas_catalog_format(tmp_path: Path):
    p = tmp_path / "catalog.tsv"
    p.write_text("MAPPED_GENE\tP-VALUE\nTP53\t1e-9\nBRCA1, BRCA2\t2e-8\n")
    df = load_gwas(p)
    # BRCA1, BRCA2 should be exploded into two rows
    assert len(df) == 3
    assert set(df["gene_id"]) == {"TP53", "BRCA1", "BRCA2"}
