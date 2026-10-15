import warnings

import networkx as nx
import numpy as np
import pytest
from shapely.geometry import LineString, box

import rd2dnet


def test_boundary_crossings_split_edges(grid):
    G, _ = grid
    bnd = rd2dnet.Boundary(LineString([(-50, 150), (450, 150)]), treated_region=box(-1e3, -1e3, 1e3, 150))
    net = rd2dnet.NetworkGraph(G).with_boundary_crossings(bnd)
    assert net.n_crossings_inserted == 5
    crossings = [n for n, d in net.graph.nodes(data=True) if d.get("rdnet_crossing")]
    assert len(crossings) == 5
    for c in crossings:
        assert net.graph.nodes[c]["y"] == pytest.approx(150.0)
        lengths = sorted(d["length"] for _, _, d in net.graph.edges(c, data=True))
        assert lengths == pytest.approx([50.0, 50.0])


def test_snap_to_edge_costs(grid):
    G, _ = grid
    net = rd2dnet.NetworkGraph(G)
    sn = net.snap([[30.0, 10.0]])
    assert sn.depth[0] == pytest.approx(10.0)
    assert sorted([sn.cost_u[0], sn.cost_v[0]]) == pytest.approx([30.0, 70.0])


def test_snap_respects_allowed_nodes(grid):
    G, _ = grid
    net = rd2dnet.NetworkGraph(G)
    allowed = net.coords[:, 1] >= 300
    sn = net.snap([[0.0, 0.0]], allowed=allowed, mode="node")
    assert net.coords[sn.u[0], 1] == pytest.approx(300.0)


def test_from_adjacency():
    A = np.array([[0, 2.0, 0], [2.0, 0, 3.0], [0, 3.0, 0]])
    net = rd2dnet.NetworkGraph.from_adjacency(A, [[0, 0], [1, 0], [2, 0]])
    assert net.node_distance(0, 2) == pytest.approx(5.0)


def test_from_place_offline_fallback(monkeypatch):
    import sys
    import types

    fake = types.ModuleType("osmnx")

    def boom(*a, **k):
        raise ConnectionError("no internet")

    fake.graph_from_place = boom
    monkeypatch.setitem(sys.modules, "osmnx", fake)
    with pytest.warns(UserWarning, match="bundled mock"):
        net = rd2dnet.NetworkGraph.from_place("Manhattan, New York")
    assert net.n_nodes == 450


def test_from_place_raises_without_fallback(monkeypatch):
    import sys
    import types

    fake = types.ModuleType("osmnx")
    fake.graph_from_place = lambda *a, **k: (_ for _ in ()).throw(ConnectionError("down"))
    monkeypatch.setitem(sys.modules, "osmnx", fake)
    with pytest.raises(ConnectionError):
        rd2dnet.NetworkGraph.from_place("Manhattan", offline_fallback=False)


def test_missing_coordinates_error():
    G = nx.path_graph(2)
    with pytest.raises(KeyError, match="coordinates"):
        rd2dnet.NetworkGraph(G)


@pytest.mark.online
def test_osm_download_midtown():  # pragma: no cover - needs internet
    net = rd2dnet.NetworkGraph.from_place("Midtown Manhattan, New York", offline_fallback=False)
    assert net.n_nodes > 100
