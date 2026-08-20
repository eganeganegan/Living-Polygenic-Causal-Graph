import logging

from neo4j import Driver

from graph.neo4j_client import run_query

logger = logging.getLogger(__name__)

_DEFAULT_WEIGHT = 0.10

# ---------------------------------------------------------------------------
# Fetch current edge weights to use as Bayesian priors before updating.
# Undirected match handles edges stored in either direction.
# ---------------------------------------------------------------------------
_FETCH_WEIGHTS = """
UNWIND $pairs AS pair
MATCH (src:Gene {gene_id: pair.src_id})-[r:INTERACTS_WITH {relationship_type: pair.rel_type}]-(tgt:Gene {gene_id: pair.tgt_id})
RETURN pair.src_id     AS src_id,
       pair.tgt_id     AS tgt_id,
       pair.rel_type   AS rel_type,
       r.statistical_weight AS statistical_weight
"""

# ---------------------------------------------------------------------------
# MERGE upsert, directed from lexicographically smaller → larger gene_id
# to prevent duplicate edges in both directions for the same pair.
#
# ON MATCH semantics:
#   statistical_weight : overwritten with the already-computed posterior
#   publication_dates  : list append, deduplicated via CASE
#   disease_tags       : list append, deduplicated via CASE
#   gwas_significant   : latched, once True it stays True across updates
# ---------------------------------------------------------------------------
_UPSERT = """
UNWIND $batch AS row
MERGE (src:Gene {gene_id: row.src_id})
  ON CREATE SET src.symbol = row.src_name
MERGE (tgt:Gene {gene_id: row.tgt_id})
  ON CREATE SET tgt.symbol = row.tgt_name
MERGE (src)-[r:INTERACTS_WITH {relationship_type: row.rel_type}]->(tgt)
  ON CREATE SET
    r.statistical_weight = row.weight,
    r.publication_dates  = [row.pub_date],
    r.disease_tags       = CASE WHEN row.disease_tag IS NULL OR row.disease_tag = '' THEN [] ELSE [row.disease_tag] END,
    r.gwas_significant   = row.gwas_sig
  ON MATCH SET
    r.statistical_weight = row.weight,
    r.publication_dates  = CASE
        WHEN row.pub_date IS NULL OR row.pub_date IN r.publication_dates
        THEN r.publication_dates
        ELSE r.publication_dates + [row.pub_date]
        END,
    r.disease_tags       = CASE
        WHEN row.disease_tag IS NULL OR row.disease_tag = '' OR row.disease_tag IN r.disease_tags
        THEN r.disease_tags
        ELSE r.disease_tags + [row.disease_tag]
        END,
    r.gwas_significant   = (r.gwas_significant OR row.gwas_sig)
"""

# ---------------------------------------------------------------------------
# Dashboard read queries
# ---------------------------------------------------------------------------

# Returns distinct disease tags present on any edge in the graph
_FETCH_AVAILABLE_DISEASES = """
MATCH ()-[r:INTERACTS_WITH]->()
UNWIND r.disease_tags AS tag
RETURN DISTINCT tag AS disease_tag
ORDER BY disease_tag
"""

# Two-stage: top 50 genes by weighted degree → all edges between them
_FETCH_DISEASE_SUBGRAPH = """
MATCH (g:Gene)-[r:INTERACTS_WITH]-(:Gene)
WHERE ANY(tag IN r.disease_tags WHERE tag = $disease_tag)
WITH g, sum(r.statistical_weight) AS wdeg
ORDER BY wdeg DESC
LIMIT 50
WITH collect(g.gene_id) AS hub_ids
MATCH (src:Gene)-[r:INTERACTS_WITH]->(tgt:Gene)
WHERE src.gene_id IN hub_ids
  AND tgt.gene_id IN hub_ids
  AND ANY(tag IN r.disease_tags WHERE tag = $disease_tag)
RETURN
  src.gene_id          AS src_id,
  src.symbol           AS src_symbol,
  tgt.gene_id          AS tgt_id,
  tgt.symbol           AS tgt_symbol,
  r.statistical_weight AS weight,
  r.gwas_significant   AS gwas_sig
"""


def _normalise(record: dict) -> dict:
    """
    Canonicalise edge direction: smaller gene_id is always the source.
    This prevents (A→B) and (B→A) from coexisting as separate edges.
    """
    src_id = str(record["source_gene_id"])
    tgt_id = str(record["target_gene_id"])
    if src_id > tgt_id:
        src_id, tgt_id = tgt_id, src_id
        src_name = record.get("target_gene_name", tgt_id)
        tgt_name = record.get("source_gene_name", src_id)
    else:
        src_name = record.get("source_gene_name", src_id)
        tgt_name = record.get("target_gene_name", tgt_id)

    return {
        "src_id":     src_id,
        "src_name":   src_name,
        "tgt_id":     tgt_id,
        "tgt_name":   tgt_name,
        "rel_type":   record.get("relationship_type", "Association"),
        "weight":     float(record.get("statistical_weight", _DEFAULT_WEIGHT)),
        "pub_date":   record.get("publication_date") or "",
        "gwas_sig":   bool(record.get("gwas_significant", False)),
        "disease_tag": record.get("disease_tag") or "",
    }


def fetch_existing_weights(driver: Driver, interactions: list[dict]) -> None:
    """
    Populate ``statistical_weight`` on each interaction dict with the value
    currently stored in Neo4j so it can serve as the Bayesian prior.

    Interactions with no existing edge are left unchanged; the Bayesian
    updater will fall back to its ``default_prior`` for those.

    Mutates ``interactions`` in-place; returns None.
    """
    if not interactions:
        return

    pairs = [
        {
            "src_id":   min(str(r["source_gene_id"]), str(r["target_gene_id"])),
            "tgt_id":   max(str(r["source_gene_id"]), str(r["target_gene_id"])),
            "rel_type": r.get("relationship_type", "Association"),
        }
        for r in interactions
    ]

    db_weights: dict[tuple[str, str, str], float] = {
        (row["src_id"], row["tgt_id"], row["rel_type"]): float(row["statistical_weight"])
        for row in run_query(driver, _FETCH_WEIGHTS, pairs=pairs)
    }

    logger.info(
        "Fetched existing weights for %d/%d edges from Neo4j.",
        len(db_weights), len(interactions),
    )

    for record in interactions:
        key = (
            min(str(record["source_gene_id"]), str(record["target_gene_id"])),
            max(str(record["source_gene_id"]), str(record["target_gene_id"])),
            record.get("relationship_type", "Association"),
        )
        if key in db_weights:
            record["statistical_weight"] = db_weights[key]


def upsert_interactions(driver: Driver, interactions: list[dict]) -> int:
    """
    Write gene-gene interactions to Neo4j using MERGE semantics.

    Assumes ``statistical_weight`` in each dict is the ready-to-write
    posterior (call ``fetch_existing_weights`` then ``update_interaction_weights``
    before this function).

    Returns the number of records submitted.
    """
    if not interactions:
        return 0

    batch = [_normalise(r) for r in interactions]
    run_query(driver, _UPSERT, batch=batch)
    logger.info("Upserted %d interaction(s) to Neo4j.", len(batch))
    return len(batch)


def fetch_available_diseases(driver: Driver) -> list[str]:
    """Return all disease tags that appear on at least one edge."""
    return [row["disease_tag"] for row in run_query(driver, _FETCH_AVAILABLE_DISEASES)]


def fetch_disease_subgraph(driver: Driver, disease_tag: str) -> list[dict]:
    """Return edge rows for the top-50 hub genes in the given disease context."""
    return run_query(driver, _FETCH_DISEASE_SUBGRAPH, disease_tag=disease_tag)
