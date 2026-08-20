import logging
import os

import dash
import dash_cytoscape as cyto

from dashboard.callbacks import register_callbacks
from dashboard.layout import build_layout
from graph.neo4j_client import get_driver
from graph.queries import fetch_available_diseases
from graph.schema import init_schema

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

# cose-bilkent gives far better layouts for dense biological networks than the default cose
cyto.load_extra_layouts()

app = dash.Dash(__name__, title="Living Polygenic Causal Graph")
server = app.server  # expose underlying Flask server for gunicorn / production


def _build_disease_options(driver) -> list[dict]:
    return [{"label": tag, "value": tag} for tag in fetch_available_diseases(driver)]


def create_app() -> dash.Dash:
    driver = get_driver()
    init_schema(driver)
    app.layout = build_layout(_build_disease_options(driver))
    register_callbacks(app, driver)
    return app


if __name__ == "__main__":
    create_app()
    app.run(
        host=os.environ.get("DASH_HOST", "127.0.0.1"),
        port=int(os.environ.get("DASH_PORT", "8050")),
        debug=os.environ.get("DASH_DEBUG", "true").lower() == "true",
    )
