"""
Hub gene scoring: computes graph-theoretic centrality metrics on the
disease-specific subgraph fetched from Neo4j.

Uses NetworkX for centrality calculations so no APOC/GDS plugin is required.
Returns a DataFrame sorted by weighted degree (gene "strength") descending.
"""

import logging

import networkx as nx
import pandas as pd
from neo4j import Driver

from graph.queries import fetch_disease_subgraph

logger = logging.getLogger(__name__)


def _build_nx_graph(rows: list[dict]) -> tuple[nx.Graph, dict[str, str]]:
    """Return (NetworkX Graph, gene_id → symbol mapping) from Neo4j edge rows."""
    G: nx.Graph = nx.Graph()
    symbols: dict[str, str] = {}
    for row in rows:
        src, tgt = row["src_id"], row["tgt_id"]
        G.add_edge(src, tgt, weight=float(row["weight"]))
        symbols.setdefault(src, row.get("src_symbol") or src)
        symbols.setdefault(tgt, row.get("tgt_symbol") or tgt)
    return G, symbols


def compute_hub_scores(driver: Driver, disease_tag: str) -> pd.DataFrame:
    """
    Return centrality metrics for every gene in the top-50 hub subgraph
    for ``disease_tag``.

    Columns:
        gene_id, symbol, degree, weighted_degree,
        degree_centrality, betweenness_centrality, eigenvector_centrality

    Returns an empty DataFrame with those columns when no data exists.
    """
    _empty = pd.DataFrame(columns=[
        "gene_id", "symbol", "degree", "weighted_degree",
        "degree_centrality", "betweenness_centrality", "eigenvector_centrality",
    ])

    rows = fetch_disease_subgraph(driver, disease_tag)
    if not rows:
        logger.info("No subgraph data for '%s'.", disease_tag)
        return _empty

    G, symbols = _build_nx_graph(rows)

    degree_c = nx.degree_centrality(G)
    between_c = nx.betweenness_centrality(G, weight="weight", normalized=True)

    try:
        eigen_c = nx.eigenvector_centrality_numpy(G, weight="weight")
    except (nx.PowerIterationFailedConvergence, nx.NetworkXException):
        # Disconnected or degenerate graphs may fail; fall back to zeros
        logger.warning("Eigenvector centrality failed for '%s'; values set to 0.", disease_tag)
        eigen_c = {n: 0.0 for n in G.nodes}

    records = [
        {
            "gene_id": n,
            "symbol": symbols.get(n, n),
            "degree": G.degree(n),
            "weighted_degree": round(sum(d["weight"] for _, _, d in G.edges(n, data=True)), 6),
            "degree_centrality": round(degree_c[n], 6),
            "betweenness_centrality": round(between_c[n], 6),
            "eigenvector_centrality": round(eigen_c[n], 6),
        }
        for n in G.nodes
    ]

    df = pd.DataFrame(records).sort_values("weighted_degree", ascending=False).reset_index(drop=True)
    logger.info("Scored %d hub genes for '%s'.", len(df), disease_tag)
    return df
