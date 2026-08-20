import json
import logging
import time
from collections.abc import Generator
from datetime import datetime, timedelta, timezone

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ingestion.parsers import extract_gene_relations

logger = logging.getLogger(__name__)

_BASE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
_EXPORT_URL = f"{_BASE_URL}/publications/export/biocjson"

# PubTator3 search only accepts text/page/size/sort, no date params.
# NCBI ESearch is used instead for date-filtered PMID retrieval.
_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

_ESEARCH_RETMAX = 500     # PMIDs per ESearch page (offset-based pagination)
_MAX_TOTAL_PMIDS = 2000   # hard cap, prevents overwhelming PubTator on broad queries
_BATCH_SIZE = 100         # max PMIDs per BioC export request
_INITIAL_BACKOFF = 1.0    # seconds
_MAX_BACKOFF = 64.0
_MAX_RATE_LIMIT_RETRIES = 7


def _build_session(ncbi_api_key: str | None = None) -> requests.Session:
    session = requests.Session()
    # Automatic retry only for transient server errors; 429 handled manually below
    adapter = HTTPAdapter(
        max_retries=Retry(
            total=3,
            backoff_factor=0.5,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods=["GET"],
        )
    )
    session.mount("https://", adapter)
    session.headers.update({"Accept": "application/json"})
    if ncbi_api_key:
        session.params = {"api_key": ncbi_api_key}  # type: ignore[assignment]
    return session


def _get(session: requests.Session, url: str, params: dict) -> dict:
    """GET with exponential backoff on HTTP 429."""
    backoff = _INITIAL_BACKOFF
    for attempt in range(_MAX_RATE_LIMIT_RETRIES):
        response = session.get(url, params=params, timeout=30)

        if response.status_code == 200:
            # strict=False permits control characters that NCBI responses sometimes contain
            return json.loads(response.text, strict=False)

        if response.status_code == 429:
            # Honour server-supplied wait time when present
            wait = float(response.headers.get("Retry-After", backoff))
            wait = min(max(wait, backoff), _MAX_BACKOFF)
            logger.warning(
                "Rate limited (429). Waiting %.1fs, attempt %d/%d.",
                wait, attempt + 1, _MAX_RATE_LIMIT_RETRIES,
            )
            time.sleep(wait)
            backoff = min(backoff * 2, _MAX_BACKOFF)
            continue

        response.raise_for_status()

    raise RuntimeError(
        f"Exceeded {_MAX_RATE_LIMIT_RETRIES} retries on rate-limit for {url}"
    )


def _search_pmids(
    session: requests.Session, disease_query: str, since: datetime
) -> list[str]:
    """
    Return date-filtered PMIDs via NCBI ESearch.

    PubTator 3's /search/ endpoint does not support date parameters and
    returns HTTP 400 when any are supplied. ESearch supports mindate/datetype
    and is the correct layer for date-scoped PubMed queries.
    """
    date_floor = since.strftime("%Y/%m/%d")
    pmids: list[str] = []
    retstart = 0

    while True:
        data = _get(session, _ESEARCH_URL, {
            "db": "pubmed",
            "term": disease_query,
            "mindate": date_floor,
            "datetype": "pdat",
            "retmax": _ESEARCH_RETMAX,
            "retstart": retstart,
            "retmode": "json",
            "sort": "date",
        })
        result = data.get("esearchresult", {})
        ids = result.get("idlist", [])
        if not ids:
            break
        pmids.extend(ids)
        retstart += len(ids)
        if len(pmids) >= _MAX_TOTAL_PMIDS:
            logger.warning(
                "PMID cap (%d) reached for '%s'. Increase _MAX_TOTAL_PMIDS or narrow the query.",
                _MAX_TOTAL_PMIDS, disease_query,
            )
            pmids = pmids[:_MAX_TOTAL_PMIDS]
            break
        if retstart >= int(result.get("count", 0)):
            break

    logger.info("Found %d PMIDs for '%s' since %s.", len(pmids), disease_query, date_floor)
    return pmids


def _iter_publications(
    session: requests.Session, pmids: list[str]
) -> Generator[dict, None, None]:
    """Yield individual BioC JSON publication dicts in batches of _BATCH_SIZE."""
    for offset in range(0, len(pmids), _BATCH_SIZE):
        chunk = pmids[offset : offset + _BATCH_SIZE]
        try:
            data = _get(session, _EXPORT_URL, {"pmids": ",".join(chunk)})
        except requests.exceptions.HTTPError as exc:
            # PubTator 3 returns 400/404 for batches it has no annotations for
            # (papers too new or not yet processed). Skip and continue.
            if exc.response is not None and exc.response.status_code in (400, 404):
                logger.warning(
                    "PubTator3 returned %d for batch at offset %d (%d PMIDs). Skipping.",
                    exc.response.status_code, offset, len(chunk),
                )
                continue
            raise
        for publication in data.get("PubTator3", []):
            yield publication


def fetch_gene_interactions(
    disease_query: str,
    lookback_hours: int = 24,
    ncbi_api_key: str | None = None,
) -> list[dict]:
    """
    Fetch gene-gene interactions from PubTator 3 for papers mentioning
    `disease_query` published in the last `lookback_hours` hours.

    Returns a list of dicts with keys:
        source_gene_id, source_gene_name,
        target_gene_id, target_gene_name,
        relationship_type, pmid, publication_date

    Args:
        disease_query:  Free-text disease term (e.g. "Schizophrenia").
        lookback_hours: Rolling window of publication recency.
        ncbi_api_key:   Optional NCBI API key (raises rate limit to 10 req/s).
    """
    since = datetime.now(timezone.utc) - timedelta(hours=lookback_hours)
    session = _build_session(ncbi_api_key)

    pmids = _search_pmids(session, disease_query, since)
    if not pmids:
        logger.info("No new publications found for '%s'.", disease_query)
        return []

    interactions: list[dict] = []
    for publication in _iter_publications(session, pmids):
        interactions.extend(extract_gene_relations(publication))

    logger.info(
        "Extracted %d gene-gene interactions from %d publications.",
        len(interactions), len(pmids),
    )
    return interactions
