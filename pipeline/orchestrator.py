"""
End-to-end ETL pipeline orchestrator.

Call run_pipeline() directly or invoke via CLI:

    python -m pipeline.orchestrator
    python -m pipeline.orchestrator --diseases "Schizophrenia,Bipolar Disorder" --lookback 48
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

import config
from graph.neo4j_client import get_driver
from graph.queries import fetch_existing_weights, upsert_interactions
from graph.schema import init_schema
from ingestion.fetch_articles import fetch_incremental
from transformation.bayesian_updater import update_interaction_weights
from transformation.gwas_integrator import load_gwas

logger = logging.getLogger(__name__)


def _gwas_path_for(disease: str, gwas_dir: Path) -> Path:
    """Derive a conventional filename from a disease name."""
    slug = disease.lower().replace(" ", "_").replace("/", "_")
    for suffix in (".tsv", ".csv", ".txt"):
        candidate = gwas_dir / f"{slug}{suffix}"
        if candidate.exists():
            return candidate
    # Return the .tsv path even if absent; load_gwas handles the missing-file case
    return gwas_dir / f"{slug}.tsv"


def run_pipeline(
    diseases: list[str] | None = None,
    gwas_dir: Path | None = None,
    state_path: Path | None = None,
    ncbi_api_key: str | None = None,
    lookback_hours: int | None = None,
) -> dict[str, int]:
    """
    Run a full incremental ETL cycle for each disease in ``diseases``.

    For each disease:
      1. Fetch new gene-gene interactions from PubTator (incremental window).
      2. Load GWAS summary stats (optional; falls back to empty DataFrame).
      3. Seed each interaction's prior from the current Neo4j edge weight.
      4. Apply Bayesian weight update (with GWAS boost where applicable).
      5. Upsert the updated interactions back to Neo4j.

    Returns a dict mapping disease → number of interactions written.
    """
    diseases = diseases or config.DISEASES
    gwas_dir = gwas_dir or config.GWAS_DIR
    state_path = state_path or config.PIPELINE_STATE_FILE
    ncbi_api_key = ncbi_api_key or config.NCBI_API_KEY or None
    lookback_hours = lookback_hours or config.LOOKBACK_HOURS

    graph = get_driver()
    init_schema(graph)

    results: dict[str, int] = {}

    for disease in diseases:
        logger.info("=" * 60)
        logger.info("Disease: %s", disease)
        logger.info("=" * 60)

        # ── 1. Ingest ────────────────────────────────────────────────
        interactions = fetch_incremental(
            disease, state_path=state_path, ncbi_api_key=ncbi_api_key,
            default_lookback_hours=lookback_hours,
        )
        if not interactions:
            logger.info("No new interactions found for '%s'. Skipping.", disease)
            results[disease] = 0
            continue

        logger.info("Fetched %d raw interactions.", len(interactions))

        # ── 2. GWAS data ─────────────────────────────────────────────
        gwas_df = load_gwas(_gwas_path_for(disease, gwas_dir))

        # ── 3. Seed priors from Neo4j ────────────────────────────────
        fetch_existing_weights(graph, interactions)

        # ── 4. Bayesian update ───────────────────────────────────────
        update_interaction_weights(interactions, gwas_df)

        boosted = sum(1 for r in interactions if r.get("gwas_significant"))
        logger.info("GWAS-boosted interactions: %d/%d.", boosted, len(interactions))

        # ── 5. Persist ───────────────────────────────────────────────
        count = upsert_interactions(graph, interactions)
        results[disease] = count

    logger.info("Pipeline complete. Summary: %s", results)
    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the LPGC ETL pipeline.")
    parser.add_argument(
        "--diseases",
        default=",".join(config.DISEASES),
        help="Comma-separated disease query strings (default: from LPGC_DISEASES env var).",
    )
    parser.add_argument(
        "--gwas-dir",
        type=Path,
        default=config.GWAS_DIR,
        help="Directory containing GWAS summary stat files.",
    )
    parser.add_argument(
        "--lookback",
        type=int,
        default=config.LOOKBACK_HOURS,
        help="Fallback lookback window in hours for first-time disease fetches.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stdout,
    )
    args = _parse_args()
    diseases = [d.strip() for d in args.diseases.split(",") if d.strip()]
    run_pipeline(diseases=diseases, gwas_dir=args.gwas_dir, lookback_hours=args.lookback)
