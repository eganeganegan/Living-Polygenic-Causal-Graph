import logging

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Genome-wide significance threshold (Bonferroni-corrected for ~1M SNPs)
# ---------------------------------------------------------------------------
GWAS_SIGNIFICANCE_THRESHOLD = 5e-8

# ---------------------------------------------------------------------------
# Likelihood parameters
# P(co-mention | true interaction),       base and GWAS-boosted variants
# P(co-mention | no true interaction),    background co-occurrence noise
# ---------------------------------------------------------------------------
_L_INTERACT = 0.80       # base sensitivity
_L_INTERACT_GWAS = 0.92  # raised when ≥1 gene is genome-wide significant
_L_NO_INTERACT = 0.10    # specificity complement (false positive rate)

_DEFAULT_PRIOR = 0.10    # cold-start prior for edges absent from the graph


def _bayesian_update(prior: float, likelihood: float, fp_rate: float) -> float:
    """
    One step of Bayes' theorem for a binary interaction hypothesis.

    posterior = P(H|E) = P(E|H)*P(H) / [P(E|H)*P(H) + P(E|¬H)*P(¬H)]
    """
    p = likelihood * prior
    return p / (p + fp_rate * (1.0 - prior))


def _gwas_significant_genes(gwas_df: pd.DataFrame) -> frozenset[str]:
    """Return gene IDs whose minimum GWAS p-value is below the GW threshold."""
    if gwas_df.empty:
        return frozenset()
    # Take per-gene minimum p-value to be robust to multi-SNP DataFrames
    min_p = gwas_df.groupby("gene_id")["p_value"].min()
    return frozenset(min_p[min_p < GWAS_SIGNIFICANCE_THRESHOLD].index.astype(str))


def update_interaction_weights(
    interactions: list[dict],
    gwas_df: pd.DataFrame,
    default_prior: float = _DEFAULT_PRIOR,
) -> list[dict]:
    """
    Apply a Bayesian weight update to each gene-gene interaction record.

    Each record's existing ``statistical_weight`` is used as the prior.
    Records lacking this key (new edges) fall back to ``default_prior``.
    The posterior is boosted when either gene in the pair carries a
    genome-wide significant GWAS association (p < 5 × 10⁻⁸).

    ``gwas_df`` must have columns:
        gene_id  (str),   NCBI Entrez Gene ID
        p_value  (float), association p-value (smallest per gene is used)

    Mutates each dict in-place by setting / overwriting:
        statistical_weight  (float), updated posterior probability
        gwas_significant    (bool), whether a GWAS hit was present

    Returns the same list for convenience.
    """
    if not interactions:
        return interactions

    sig_genes = _gwas_significant_genes(gwas_df)
    logger.info(
        "%d genome-wide significant gene(s) found in GWAS DataFrame.",
        len(sig_genes),
    )

    updated = 0
    for record in interactions:
        prior = float(record.get("statistical_weight", default_prior))

        # Clamp prior to open interval (0, 1), log-odds undefined at boundaries
        prior = max(1e-6, min(prior, 1.0 - 1e-6))

        gwas_hit = bool(
            str(record.get("source_gene_id", "")) in sig_genes
            or str(record.get("target_gene_id", "")) in sig_genes
        )

        likelihood = _L_INTERACT_GWAS if gwas_hit else _L_INTERACT
        posterior = _bayesian_update(prior, likelihood, _L_NO_INTERACT)

        # Clamp after rounding: round() can promote 0.9999999… to exactly 1.0
        record["statistical_weight"] = min(round(posterior, 6), 1.0 - 1e-6)
        record["gwas_significant"] = gwas_hit
        updated += 1

    logger.info(
        "Updated weights for %d interaction(s). GWAS-boosted: %d.",
        updated,
        sum(r["gwas_significant"] for r in interactions),
    )
    return interactions
