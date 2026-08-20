import logging

from neo4j import Driver, GraphDatabase

import config

logger = logging.getLogger(__name__)

_driver: Driver | None = None


def get_driver(
    uri: str | None = None,
    user: str | None = None,
    password: str | None = None,
) -> Driver:
    """Return a module-level cached neo4j Driver. Falls back to config values."""
    global _driver
    if _driver is None:
        _uri = uri or config.NEO4J_URI
        _user = user or config.NEO4J_USER
        # Use `is not None` so an explicit empty-string password is respected
        _pwd = password if password is not None else config.NEO4J_PASSWORD
        _driver = GraphDatabase.driver(_uri, auth=(_user, _pwd))
        logger.info("Connected to Neo4j at %s as '%s'.", _uri, _user)
    return _driver


# Backwards-compatible alias
get_graph = get_driver


def reset_driver() -> None:
    """Force a fresh connection on the next get_driver() call (useful in tests)."""
    global _driver
    if _driver:
        _driver.close()
    _driver = None


# Backwards-compatible alias
reset_graph = reset_driver


def run_query(driver: Driver, query: str, **params) -> list[dict]:
    """Execute a read or write query and return results as a list of dicts."""
    with driver.session() as session:
        return session.run(query, **params).data()
