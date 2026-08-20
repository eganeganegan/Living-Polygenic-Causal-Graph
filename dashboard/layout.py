import dash_cytoscape as cyto
from dash import dcc, html

_DARK_BG = "#1A1A2E"
_PANEL_BG = "#16213E"
_GRAPH_BG = "#0F0F23"
_TEXT = "#EAEAEA"


def build_layout(disease_options: list[dict]) -> html.Div:
    return html.Div(
        style={"fontFamily": "Inter, sans-serif", "backgroundColor": _DARK_BG, "minHeight": "100vh", "color": _TEXT},
        children=[
            # ---- header --------------------------------------------------
            html.Div(
                style={"backgroundColor": _PANEL_BG, "padding": "14px 24px", "borderBottom": "1px solid #2A2A4A"},
                children=[
                    html.H2("Living Polygenic Causal Graph", style={"margin": 0, "fontSize": "1.3rem"}),
                    html.P(
                        "Top 50 hub genes by weighted degree · node size = degree centrality · amber = GWAS-significant",
                        style={"margin": "4px 0 0", "fontSize": "0.78rem", "color": "#9090B0"},
                    ),
                ],
            ),
            # ---- controls ------------------------------------------------
            html.Div(
                style={"display": "flex", "alignItems": "center", "gap": "12px", "padding": "12px 24px", "backgroundColor": _PANEL_BG, "borderBottom": "1px solid #2A2A4A"},
                children=[
                    html.Label("Disease context:", style={"whiteSpace": "nowrap", "fontSize": "0.88rem"}),
                    dcc.Dropdown(
                        id="disease-select",
                        options=disease_options,
                        placeholder="Select a disease…",
                        clearable=True,
                        style={"width": "340px", "color": "#1A1A2E", "fontSize": "0.88rem"},
                    ),
                    html.Div(id="status-label", style={"fontSize": "0.78rem", "color": "#9090B0", "marginLeft": "auto"}),
                ],
            ),
            # ---- graph canvas --------------------------------------------
            cyto.Cytoscape(
                id="lpgc-graph",
                elements=[],
                layout={"name": "cose-bilkent", "animate": True, "randomize": False, "nodeDimensionsIncludeLabels": True},
                style={"width": "100%", "height": "82vh", "backgroundColor": _GRAPH_BG},
                stylesheet=[],  # populated by callback
                responsive=True,
            ),
        ],
    )
