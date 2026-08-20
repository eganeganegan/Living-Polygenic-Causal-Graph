"""
Orchestrates incremental PubTator fetches per disease.

Maintains a state file (pipeline_state.json) that records the ISO-8601
timestamp of the last successful fetch for each disease. On subsequent
runs the lookback window is calculated from that timestamp rather than
a fixed number of hours, preventing both gaps and duplicate fetches.
"""

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from ingestion.pubtator_client import fetch_gene_interactions

logger = logging.getLogger(__name__)


def _load_state(state_path: Path) -> dict:
    if state_path.exists():
        with state_path.open() as f:
            return json.load(f)
    return {}


def _save_state(state_path: Path, state: dict) -> None:
    state_path.parent.mkdir(parents=True, exist_ok=True)
    with state_path.open("w") as f:
        json.dump(state, f, indent=2)


def _hours_since(iso_timestamp: str) -> float:
    last = datetime.fromisoformat(iso_timestamp)
    delta = datetime.now(timezone.utc) - last
    return delta.total_seconds() / 3600


def fetch_and_tag(
    disease: str,
    lookback_hours: int = 24,
    ncbi_api_key: str | None = None,
) -> list[dict]:
    """
    Fetch gene interactions for ``disease`` and stamp every record with
    ``disease_tag`` so downstream Neo4j writes can filter by disease.
    """
    interactions = fetch_gene_interactions(disease, lookback_hours=lookback_hours, ncbi_api_key=ncbi_api_key)
    for record in interactions:
        record["disease_tag"] = disease
    return interactions


def fetch_incremental(
    disease: str,
    state_path: Path,
    default_lookback_hours: int = 24,
    ncbi_api_key: str | None = None,
) -> list[dict]:
    """
    Fetch interactions since the last recorded run for ``disease``.

    Falls back to ``default_lookback_hours`` on the first run.
    Updates ``state_path`` only after a successful fetch.
    """
    state = _load_state(state_path)
    last_run = state.get(disease)

    if last_run:
        lookback = _hours_since(last_run)
        logger.info("Incremental fetch for '%s': %.1f hours since last run.", disease, lookback)
    else:
        lookback = float(default_lookback_hours)
        logger.info("First fetch for '%s': using default lookback of %g hours.", disease, lookback)

    interactions = fetch_and_tag(disease, lookback_hours=int(lookback) + 1, ncbi_api_key=ncbi_api_key)

    # Record run time only after successful API response
    state[disease] = datetime.now(timezone.utc).isoformat()
    _save_state(state_path, state)

    return interactions
