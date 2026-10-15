import networkx as nx
import numpy as np
import pytest
from shapely.geometry import LineString, box

import rdnet
from rdnet.datasets import load_mock_midtown


@pytest.fixture(scope="session")
def city():
    return load_mock_midtown()


@pytest.fixture(scope="session")
def coords(city):
    return rdnet.network_coordinates(city.network("length"), city.boundary, city.obs[["x", "y"]], n_eval=21)


@pytest.fixture
def grid():
    """5x5 unit grid (100 m spacing) with a boundary along y = 200, treated below."""
    G = nx.grid_2d_graph(5, 5)
    G = nx.convert_node_labels_to_integers(G, label_attribute="ij")
    for n, d in G.nodes(data=True):
        d["x"], d["y"] = 100.0 * d["ij"][0], 100.0 * d["ij"][1]
    for u, v, d in G.edges(data=True):
        d["length"] = 100.0
    line = LineString([(-50, 200), (450, 200)])
    bnd = rdnet.Boundary(line, treated_region=box(-1e3, -1e3, 1e3, 200))
    return G, bnd


@pytest.fixture
def highway_graph():
    """Highway A -- ramp -- highway B, with a local street hanging off the ramp.

    a0 --hw-- a1 --ramp-- r --ramp-- b1 --hw-- b0
                          |
                        local -- l1
    """
    G = nx.Graph()
    pos = {"a0": (0, 0), "a1": (100, 0), "r": (200, 0), "b1": (300, 0), "b0": (400, 0), "l1": (200, 100)}
    for n, (x, y) in pos.items():
        G.add_node(n, x=x, y=y)
    G.add_edge("a0", "a1", length=100.0, highway="motorway")
    G.add_edge("a1", "r", length=100.0, highway="motorway_link")
    G.add_edge("r", "b1", length=100.0, highway="motorway_link")
    G.add_edge("b1", "b0", length=100.0, highway="motorway")
    G.add_edge("r", "l1", length=100.0, highway="residential")
    return G
