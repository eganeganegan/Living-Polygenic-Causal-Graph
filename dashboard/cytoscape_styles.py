# Dark-mode palette: standard nodes in slate-blue, GWAS-significant in amber

STYLESHEET = [
    {
        "selector": "node",
        "style": {
            "label": "data(label)",
            "width": "data(size)",
            "height": "data(size)",
            "background-color": "#4A90D9",
            "color": "#FFFFFF",
            "font-size": 10,
            "text-halign": "center",
            "text-valign": "center",
            "text-wrap": "ellipsis",
            "text-max-width": 80,
        },
    },
    {
        # Genome-wide significant genes highlighted in amber
        "selector": "node[gwas_sig = 1]",
        "style": {
            "background-color": "#E8A838",
            "border-color": "#B87C18",
            "border-width": 2,
        },
    },
    {
        "selector": "node:selected",
        "style": {
            "border-color": "#FFFFFF",
            "border-width": 3,
            "overlay-opacity": 0.1,
        },
    },
    {
        "selector": "edge",
        "style": {
            # mapData: linearly interpolate weight [0,1] → visual range
            "width": "mapData(weight, 0.0, 1.0, 1, 8)",
            "opacity": "mapData(weight, 0.0, 1.0, 0.15, 0.85)",
            "line-color": "#6C8EBF",
            "curve-style": "bezier",
        },
    },
    {
        "selector": "edge:selected",
        "style": {
            "line-color": "#FFFFFF",
            "opacity": 1.0,
        },
    },
]
