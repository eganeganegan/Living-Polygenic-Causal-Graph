from typing import Any


def _build_gene_index(publication: dict) -> dict[str, dict]:
    """Map passage-level annotation ID → gene metadata for name resolution."""
    genes: dict[str, dict] = {}
    for passage in publication.get("passages", []):
        for annotation in passage.get("annotations", []):
            infons = annotation.get("infons", {})
            if infons.get("type") != "Gene":
                continue
            ann_id = annotation.get("id")
            if ann_id:
                genes[ann_id] = {
                    "gene_id": infons.get("identifier", ""),
                    "name": infons.get("name") or annotation.get("text", ""),
                }
    return genes


def extract_gene_relations(publication: dict[str, Any]) -> list[dict]:
    """
    Parse a single BioC JSON publication dict into gene-gene interaction records.

    PubTator 3 encodes relation roles in two ways; both are handled:
      1. infons.role1 / infons.role2  (preferred, role-keyed)
      2. nodes[0] / nodes[1]          (fallback, index-keyed)

    Skips: non-Gene roles, missing identifiers, self-loops.
    """
    pmid = str(publication.get("id", ""))
    pub_date = publication.get("date", "")
    gene_index = _build_gene_index(publication)
    records: list[dict] = []

    for relation in publication.get("relations", []):
        infons = relation.get("infons", {})
        relationship_type = infons.get("type", "Association")

        # --- resolve the two roles ----------------------------------------
        role1 = infons.get("role1")
        role2 = infons.get("role2")

        if not (role1 and role2):
            # Fallback: derive roles from the nodes list + gene_index
            nodes = relation.get("nodes", [])
            if len(nodes) < 2:
                continue
            role1 = gene_index.get(nodes[0].get("refid", ""))
            role2 = gene_index.get(nodes[1].get("refid", ""))
            if not (role1 and role2):
                continue
            # Nodes-path already carries gene_id; synthesise type field
            role1 = {"type": "Gene", "identifier": role1["gene_id"], "refid": nodes[0]["refid"]}
            role2 = {"type": "Gene", "identifier": role2["gene_id"], "refid": nodes[1]["refid"]}

        if role1.get("type") != "Gene" or role2.get("type") != "Gene":
            continue

        src_id = role1.get("identifier", "").strip()
        tgt_id = role2.get("identifier", "").strip()

        if not src_id or not tgt_id or src_id == tgt_id:
            continue

        # Prefer name from infons; fall back to passage-level annotation index
        src_name = (
            role1.get("name")
            or gene_index.get(role1.get("refid", ""), {}).get("name")
            or src_id
        )
        tgt_name = (
            role2.get("name")
            or gene_index.get(role2.get("refid", ""), {}).get("name")
            or tgt_id
        )

        records.append({
            "source_gene_id": src_id,
            "source_gene_name": src_name,
            "target_gene_id": tgt_id,
            "target_gene_name": tgt_name,
            "relationship_type": relationship_type,
            "pmid": pmid,
            "publication_date": pub_date,
        })

    return records
