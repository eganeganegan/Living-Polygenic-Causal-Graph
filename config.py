import os
from pathlib import Path

from dotenv import load_dotenv

# Load .env when present (development); silently a no-op in CI / production
load_dotenv()

ROOT = Path(__file__).parent

# --- Neo4j -------------------------------------------------------------------
NEO4J_URI: str = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER: str = os.environ.get("NEO4J_USER", "neo4j")

# Fail-fast: an unset password is almost certainly a misconfiguration
_raw_pwd = os.environ.get("NEO4J_PASSWORD")
if _raw_pwd is None:
    raise EnvironmentError(
        "NEO4J_PASSWORD is not set. "
        "Copy .env.example to .env and fill in your credentials. "
        "Set NEO4J_PASSWORD= (empty) if your instance has no password."
    )
NEO4J_PASSWORD: str = _raw_pwd

# --- NCBI / PubTator ---------------------------------------------------------
# Optional, registers at https://www.ncbi.nlm.nih.gov/account/
# Raises the API rate limit from 3 req/s to 10 req/s.
NCBI_API_KEY: str = os.environ.get("NCBI_API_KEY", "")

# --- Data directories --------------------------------------------------------
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
GWAS_DIR = DATA_DIR / "gwas"

# --- Pipeline settings -------------------------------------------------------
DISEASES: list[str] = [
    d.strip() for d in os.environ.get("LPGC_DISEASES", "Schizophrenia").split(",") if d.strip()
]
LOOKBACK_HOURS: int = int(os.environ.get("LPGC_LOOKBACK_HOURS", "24"))

# State file tracks per-disease last-run timestamps for incremental fetching
PIPELINE_STATE_FILE = RAW_DIR / "pipeline_state.json"
