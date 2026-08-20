"""
Load GWAS summary statistics into a normalised DataFrame for use by the
Bayesian updater.

Two input formats are supported:
  1. Simple format, any TSV/CSV with 'gene_id' and 'p_value' columns.
  2. GWAS Catalog,   the full GWAS Catalog download TSV (v1.0.3).
                       Uses 'MAPPED_GENE' and 'P-VALUE' columns.
                       Gene symbols in MAPPED_GENE are preserved as-is;
                       callers must map to Entrez IDs if needed.

Output DataFrame always has exactly two columns:
    gene_id  (str), identifier for the gene
    p_value  (float), association p-value
"""

import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# GWAS Catalog column names
_CATALOG_GENE_COL = "MAPPED_GENE"
_CATALOG_PVAL_COL = "P-VALUE"

# Separators used between multiple genes in GWAS Catalog MAPPED_GENE column
_GENE_SEP = re.compile(r"[,;\s]+")


def _explode_gene_column(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    """Split multi-gene cells (e.g. 'BRCA1, TP53') into one row per gene."""
    df = df.copy()
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df[gene_col] = df[gene_col].str.split(_GENE_SEP)
    return df.explode(gene_col)


def load_gwas(path: Path, gene_col: str = "gene_id", pval_col: str = "p_value") -> pd.DataFrame:
    """
    Load a GWAS summary stats file and return a tidy DataFrame.

    Automatically detects GWAS Catalog format when the file contains a
    'MAPPED_GENE' column; otherwise expects ``gene_col`` and ``pval_col``.

    Rows with missing or non-numeric p-values are silently dropped.
    """
    path = Path(path)
    if not path.exists():
        logger.warning("GWAS file not found: %s. Returning empty DataFrame.", path)
        return pd.DataFrame(columns=["gene_id", "p_value"])

    sep = "\t" if path.suffix in (".tsv", ".txt") else ","
    raw = pd.read_csv(path, sep=sep, low_memory=False)

    # Auto-detect GWAS Catalog format
    if _CATALOG_GENE_COL in raw.columns and _CATALOG_PVAL_COL in raw.columns:
        logger.info("Detected GWAS Catalog format in %s.", path.name)
        df = raw[[_CATALOG_GENE_COL, _CATALOG_PVAL_COL]].rename(
            columns={_CATALOG_GENE_COL: "gene_id", _CATALOG_PVAL_COL: "p_value"}
        )
        df = _explode_gene_column(df, "gene_id")
    elif gene_col in raw.columns and pval_col in raw.columns:
        df = raw[[gene_col, pval_col]].rename(columns={gene_col: "gene_id", pval_col: "p_value"})
    else:
        raise ValueError(
            f"Expected columns '{gene_col}' and '{pval_col}' (or GWAS Catalog columns) "
            f"in {path.name}. Found: {list(raw.columns)}"
        )

    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    before = len(df)
    df = df.dropna(subset=["gene_id", "p_value"])
    df = df[df["gene_id"].str.len() > 0]
    dropped = before - len(df)
    if dropped:
        logger.debug("Dropped %d rows with invalid gene_id or p_value.", dropped)

    logger.info("Loaded %d gene-pvalue rows from %s.", len(df), path.name)
    return df.reset_index(drop=True)
