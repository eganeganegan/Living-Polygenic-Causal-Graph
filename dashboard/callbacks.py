import logging
from collections import Counter, defaultdict

from dash import Input, Output, callback, no_update
from neo4j import Driver

from dashboard.cytoscape_styles import STYLESHEET
from graph.queries import fetch_disease_subgraph

logger = logging.getLogger(__name__)

_SIZE_MIN = 18.0
_SIZE_MAX = 72.0


def _normalise_sizes(degree: Counter) -> dict[str, float]:
    """Map raw degree counts to the [_SIZE_MIN, _SIZE_MAX] pixel range."""
    lo, hi = min(degree.values()), max(degree.values())
    span = hi - lo
    return {
        gene_id: _SIZE_MIN if span == 0 else _SIZE_MIN + (d - lo) / span * (_SIZE_MAX - _SIZE_MIN)
        for gene_id, d in degree.items()
    }


def _build_elements(rows: list[dict]) -> list[dict]:
    """Convert Neo4j edge rows into a flat Cytoscape elements list."""
    if not rows:
        return []

    degree: Counter = Counter()
    # gwas_sig is an edge property; latch it onto each endpoint node
    gwas_by_node: dict[str, bool] = defaultdict(bool)

    for row in rows:
        degree[row["src_id"]] += 1
        degree[row["tgt_id"]] += 1
        gwas_by_node[row["src_id"]] |= bool(row.get("gwas_sig"))
        gwas_by_node[row["tgt_id"]] |= bool(row.get("gwas_sig"))

    sizes = _normalise_sizes(degree)
    elements: list[dict] = []
    seen_nodes: set[str] = set()
    seen_edges: set[tuple[str, str]] = set()

    for row in rows:
        for gene_id, symbol in ((row["src_id"], row["src_symbol"]), (row["tgt_id"], row["tgt_symbol"])):
            if gene_id not in seen_nodes:
                elements.append({
                    "data": {
                        "id": gene_id,
                        "label": symbol or gene_id,
                        "size": round(sizes[gene_id], 1),
                        "degree": degree[gene_id],
                        # Integer flag used by Cytoscape selector [gwas_sig = 1]
                        "gwas_sig": int(gwas_by_node[gene_id]),
                    }
                })
                seen_nodes.add(gene_id)

        edge_key = (row["src_id"], row["tgt_id"])
        if edge_key not in seen_edges:
            elements.append({
                "data": {
                    "source": row["src_id"],
                    "target": row["tgt_id"],
                    "weight": round(float(row["weight"]), 4),
                }
            })
            seen_edges.add(edge_key)

    return elements


def register_callbacks(app, driver: Driver) -> None:

    @app.callback(
        Output("lpgc-graph", "elements"),
        Output("lpgc-graph", "stylesheet"),
        Output("status-label", "children"),
        Input("disease-select", "value"),
    )
    def update_graph(disease_tag: str | None):
        if not disease_tag:
            return [], STYLESHEET, "Select a disease to load the graph."

        try:
            rows = fetch_disease_subgraph(driver, disease_tag)
        except Exception:
            logger.exception("Neo4j query failed for disease '%s'.", disease_tag)
            return no_update, no_update, "Error querying Neo4j, check logs."

        if not rows:
            return [], STYLESHEET, f"No interaction data found for '{disease_tag}'."

        elements = _build_elements(rows)
        node_count = sum(1 for e in elements if "source" not in e["data"])
        edge_count = len(elements) - node_count
        status = f"{node_count} genes · {edge_count} interactions · {disease_tag}"
        return elements, STYLESHEET, status
