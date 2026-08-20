"""
Pytest configuration and shared fixtures.

Stubs out the neo4j driver at import time so graph module tests run without a
live Neo4j installation or the neo4j package installed.
"""

import sys
from unittest.mock import MagicMock

# Insert a stub before any test module imports graph.queries or graph.neo4j_client
if "neo4j" not in sys.modules:
    neo4j_stub = MagicMock()
    # Provide a real type so isinstance / type-hint checks don't explode
    neo4j_stub.Driver = type("Driver", (), {})
    sys.modules["neo4j"] = neo4j_stub
