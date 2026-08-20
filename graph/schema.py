"""
Neo4j schema initialization: constraints and indexes for the LPGC graph.

Requires Neo4j 4.4+ (Community or Enterprise).
Relationship property indexes require Neo4j 4.3+.
"""

from neo4j import Driver

# ---------------------------------------------------------------------------
# Uniqueness constraints, also implicitly create a backing b-tree index
# ---------------------------------------------------------------------------

CONSTRAINT_GENE_ID = """
CREATE CONSTRAINT constraint_gene_id IF NOT EXISTS
FOR (g:Gene) REQUIRE g.gene_id IS UNIQUE
"""

CONSTRAINT_DISEASE_ID = """
CREATE CONSTRAINT constraint_disease_id IF NOT EXISTS
FOR (d:Disease) REQUIRE d.disease_id IS UNIQUE
"""

# ---------------------------------------------------------------------------
# Node property indexes for frequent lookup / text search fields
# ---------------------------------------------------------------------------

INDEX_GENE_SYMBOL = """
CREATE INDEX index_gene_symbol IF NOT EXISTS
FOR (g:Gene) ON (g.symbol)
"""

INDEX_GENE_TAXON = """
CREATE INDEX index_gene_taxon IF NOT EXISTS
FOR (g:Gene) ON (g.taxon_id)
"""

INDEX_DISEASE_NAME = """
CREATE INDEX index_disease_name IF NOT EXISTS
FOR (d:Disease) ON (d.name)
"""

# ---------------------------------------------------------------------------
# Relationship property indexes on INTERACTS_WITH
# Enables efficient range queries on weight thresholds and date windows
# ---------------------------------------------------------------------------

INDEX_INTERACTS_WEIGHT = """
CREATE INDEX index_interacts_weight IF NOT EXISTS
FOR ()-[r:INTERACTS_WITH]-() ON (r.statistical_weight)
"""

INDEX_INTERACTS_DATE = """
CREATE INDEX index_interacts_date IF NOT EXISTS
FOR ()-[r:INTERACTS_WITH]-() ON (r.publication_date)
"""

# ---------------------------------------------------------------------------
# Ordered list used by init_schema(), constraints must precede indexes
# ---------------------------------------------------------------------------

SCHEMA_STATEMENTS = [
    CONSTRAINT_GENE_ID,
    CONSTRAINT_DISEASE_ID,
    INDEX_GENE_SYMBOL,
    INDEX_GENE_TAXON,
    INDEX_DISEASE_NAME,
    INDEX_INTERACTS_WEIGHT,
    INDEX_INTERACTS_DATE,
]


def init_schema(driver: Driver) -> None:
    """Apply all constraints and indexes idempotently."""
    with driver.session() as session:
        for statement in SCHEMA_STATEMENTS:
            session.run(statement)
