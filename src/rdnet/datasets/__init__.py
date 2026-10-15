"""Bundled offline datasets.

``load_mock_midtown()`` reads a small synthetic Midtown-Manhattan-style street
grid shipped with the package, so examples and tests never need internet
access. ``make_mock_midtown()`` regenerates it (or variants) from scratch.

Layout (metres): 11 north-south avenues 250 m apart (``primary``), 41 east-west
streets 80 m apart numbered 40th-80th (``residential``), a west-side highway
(``motorway``) joined to the grid by ramps (``motorway_link``), and a park
north of the boundary between avenues 3 and 7 with no interior streets or
observations. The boundary is 60th Street; the treated zone lies south of it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from typing import Callable

import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, box

from ..boundary import Boundary
from ..network import NetworkGraph

AVENUE_SPACING = 250.0
STREET_SPACING = 80.0
N_AVENUES = 11
FIRST_STREET, LAST_STREET = 40, 80
BOUNDARY_STREET = 60
HIGHWAY_X = -300.0
RAMP_STREETS = (42, 50, 57, 66, 72, 79)
SPEED_KPH = {"residential": 25.0, "primary": 40.0, "motorway": 88.0, "motorway_link": 40.0}

BOUNDARY_Y = (BOUNDARY_STREET - FIRST_STREET) * STREET_SPACING
PARK = box(3 * AVENUE_SPACING, BOUNDARY_Y, 7 * AVENUE_SPACING, BOUNDARY_Y + 15 * STREET_SPACING)
BOUNDARY_LINE = LineString([(HIGHWAY_X - 100, BOUNDARY_Y), ((N_AVENUES - 1) * AVENUE_SPACING + 100, BOUNDARY_Y)])
TREATED_REGION = box(-1e5, -1e5, 1e5, BOUNDARY_Y)


def true_tau(arc) -> np.ndarray:
    """Treatment effect used to simulate the outcome, by arc length along the boundary."""
    x = np.asarray(arc, dtype=float) + BOUNDARY_LINE.coords[0][0]
    return 1.5 + 0.75 * np.sin(2 * np.pi * x / ((N_AVENUES - 1) * AVENUE_SPACING))


def mock_boundary(mask_park: bool = True) -> Boundary:
    return Boundary(BOUNDARY_LINE, treated_region=TREATED_REGION, mask=[PARK] if mask_park else None)


@dataclass
class MockCity:
    graph: nx.Graph
    obs: pd.DataFrame
    boundary: Boundary
    park: object
    true_tau: Callable

    def network(self, metric="length", **kwargs) -> NetworkGraph:
        return NetworkGraph(self.graph, metric, **kwargs)


def _grid_id(i: int, j: int) -> int:
    return i * 100 + j


def make_mock_graph() -> nx.Graph:
    G = nx.Graph(name="mock_midtown", crs="local-metres")
    streets = range(FIRST_STREET, LAST_STREET + 1)

    def in_park(i, j):
        return 3 < i < 7 and BOUNDARY_STREET < j < BOUNDARY_STREET + 15

    rng = np.random.default_rng(7)
    for i in range(N_AVENUES):
        for j in streets:
            if in_park(i, j):
                continue
            far = float(np.round(np.exp(rng.normal(1.5 - 0.02 * abs(j - 52), 0.4)), 3))
            G.add_node(_grid_id(i, j), x=i * AVENUE_SPACING, y=(j - FIRST_STREET) * STREET_SPACING, far=far)

    def add(a, b, cls):
        xa, ya = G.nodes[a]["x"], G.nodes[a]["y"]
        xb, yb = G.nodes[b]["x"], G.nodes[b]["y"]
        length = float(np.hypot(xa - xb, ya - yb))
        speed = SPEED_KPH[cls]
        G.add_edge(a, b, length=length, highway=cls, speed_kph=speed, travel_time=length / (speed / 3.6))

    for i in range(N_AVENUES):
        for j in streets:
            if _grid_id(i, j) not in G:
                continue
            if i + 1 < N_AVENUES and _grid_id(i + 1, j) in G:
                add(_grid_id(i, j), _grid_id(i + 1, j), "residential")
            if j + 1 <= LAST_STREET and _grid_id(i, j + 1) in G:
                add(_grid_id(i, j), _grid_id(i, j + 1), "primary")
    for j in streets:
        G.add_node(9000 + j, x=HIGHWAY_X, y=(j - FIRST_STREET) * STREET_SPACING, far=0.0)
        if j > FIRST_STREET:
            add(9000 + j - 1, 9000 + j, "motorway")
    for j in RAMP_STREETS:
        add(9000 + j, _grid_id(0, j), "motorway_link")
    return G


def make_mock_midtown(n_obs: int = 3000, seed: int = 0) -> MockCity:
    """Generate the mock city with ``n_obs`` simulated parcels."""
    G = make_mock_graph()
    rng = np.random.default_rng(seed)
    xmax = (N_AVENUES - 1) * AVENUE_SPACING
    ymax = (LAST_STREET - FIRST_STREET) * STREET_SPACING
    pts = []
    while len(pts) < n_obs:
        cand = rng.uniform([0, 0], [xmax, ymax], size=(n_obs, 2))
        inside = (cand[:, 0] > PARK.bounds[0]) & (cand[:, 0] < PARK.bounds[2]) & (cand[:, 1] > PARK.bounds[1]) & (
            cand[:, 1] < PARK.bounds[3]
        )
        pts.extend(cand[~inside].tolist())
    xy = np.array(pts[:n_obs])
    x, y = xy[:, 0], xy[:, 1]
    treated = y < BOUNDARY_Y
    arc = x - BOUNDARY_LINE.coords[0][0]
    far = np.round(np.exp(rng.normal(1.2, 0.5, n_obs)), 3)
    f = 0.0008 * (y - BOUNDARY_Y) + 0.6 * np.cos(x / 700) + 0.2 * np.log(far)
    outcome = f + true_tau(arc) * treated + rng.normal(0, 1, n_obs)
    obs = pd.DataFrame(
        {
            "x": np.round(x, 3),
            "y": np.round(y, 3),
            "treated": treated,
            "far": far,
            "block": (np.floor(x / AVENUE_SPACING) * 1000 + np.floor(y / STREET_SPACING)).astype(int),
            "outcome": np.round(outcome, 6),
        }
    )
    return MockCity(G, obs, mock_boundary(), PARK, true_tau)


def load_mock_midtown() -> MockCity:
    """Load the bundled mock Midtown network and observations (no internet needed)."""
    data = resources.files(__package__) / "data"
    graph_json = json.loads((data / "mock_midtown_graph.json").read_text())
    G = nx.node_link_graph(graph_json, edges="edges")
    with (data / "mock_midtown_obs.csv").open("r") as fh:
        obs = pd.read_csv(fh)
    obs["treated"] = obs["treated"].astype(bool)
    return MockCity(G, obs, mock_boundary(), PARK, true_tau)


def write_bundled_data(path=None) -> None:
    """Regenerate the files under ``rdnet/datasets/data`` (maintainers only)."""
    from pathlib import Path

    out = Path(path) if path else Path(__file__).parent / "data"
    out.mkdir(parents=True, exist_ok=True)
    city = make_mock_midtown()
    (out / "mock_midtown_graph.json").write_text(json.dumps(nx.node_link_data(city.graph, edges="edges")))
    city.obs.to_csv(out / "mock_midtown_obs.csv", index=False)
