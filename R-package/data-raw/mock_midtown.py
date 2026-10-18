"""Export the Python package's bundled mock city to CSV (run from the repo root).

Edges are written in networkx's iteration order so that the R and Python
packages see the same edge order.
"""
import sys
from pathlib import Path

import pandas as pd

from rd2dnet.datasets import load_mock_midtown

out = Path(sys.argv[1] if len(sys.argv) > 1 else "R-package/data-raw")
city = load_mock_midtown()
G = city.graph
nodes = pd.DataFrame([{"name": n, "x": d["x"], "y": d["y"], "far": d["far"]} for n, d in G.nodes(data=True)])
edges = pd.DataFrame([{"from": u, "to": v, **{k: d[k] for k in ("length", "highway", "speed_kph", "travel_time")}}
                      for u, v, d in G.edges(data=True)])
nodes.to_csv(out / "mock_midtown_nodes.csv", index=False, float_format="%.17g")
edges.to_csv(out / "mock_midtown_edges.csv", index=False, float_format="%.17g")
city.obs.to_csv(out / "mock_midtown_obs.csv", index=False, float_format="%.17g")
print(len(nodes), len(edges), len(city.obs))
