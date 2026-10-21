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


def _fake_osmnx(monkeypatch, geom=None, *, fail=None):
    """Stand-in for osmnx: geocoding returns ``geom`` (or raises ``fail``) and
    downloads return the mock city's graph, recording which loader ran."""
    import sys
    import types

    import pandas as pd
    from shapely.geometry import Point

    from rd2dnet.datasets import make_mock_graph

    calls = []
    fake = types.ModuleType("osmnx")

    def geocode_to_gdf(query, which_result=None):
        if fail is not None:
            raise fail
        calls.append(("geocode", which_result))
        g = geom if geom is not None else Point(-73.98, 40.75)
        # Polygons are given in metres, so the "projected" area is the shape's own.
        return pd.DataFrame({"geometry": [g], "display_name": ["Somewhere"], "lat": [40.75], "lon": [-73.98]})

    def loader(name):
        def load(*args, **kwargs):
            calls.append((name, args, kwargs))
            return make_mock_graph()
        return load

    fake.geocode_to_gdf = geocode_to_gdf
    fake.projection = types.SimpleNamespace(
        project_gdf=lambda gdf: pd.DataFrame({"area": [gdf.geometry.iloc[0].area]}))
    fake.graph_from_polygon = loader("polygon")
    fake.graph_from_point = loader("point")
    fake.graph_from_bbox = loader("bbox")
    monkeypatch.setitem(sys.modules, "osmnx", fake)
    return calls


def test_from_place_offline_fallback(monkeypatch):
    _fake_osmnx(monkeypatch, fail=ConnectionError("no internet"))
    with pytest.warns(UserWarning, match="bundled mock"):
        net = rd2dnet.NetworkGraph.from_place("Manhattan, New York")
    assert net.n_nodes == 450


def test_from_place_raises_without_fallback(monkeypatch):
    _fake_osmnx(monkeypatch, fail=ConnectionError("down"))
    with pytest.raises(ConnectionError):
        rd2dnet.NetworkGraph.from_place("Manhattan", offline_fallback=False)


def test_from_place_uses_the_outline_of_a_large_area(monkeypatch):
    calls = _fake_osmnx(monkeypatch, box(0, 0, 2000, 2000))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        net = rd2dnet.NetworkGraph.from_place("Manhattan, New York")
    assert net.n_nodes == 450
    assert calls[0] == ("geocode", 1)  # the top match, not the first polygon further down
    assert calls[1][0] == "polygon"


@pytest.mark.parametrize("geom", [None, box(0, 0, 100, 120)], ids=["point", "building"])
def test_from_place_uses_a_radius_around_points_and_small_areas(monkeypatch, geom):
    calls = _fake_osmnx(monkeypatch, geom)
    with pytest.warns(UserWarning, match="within 800 m"):
        rd2dnet.NetworkGraph.from_place("Midtown Manhattan, New York", dist=800.0)
    name, args, kwargs = calls[1]
    assert name == "point"
    assert args == ((40.75, -73.98),) and kwargs["dist"] == 800.0


def test_from_place_bbox(monkeypatch):
    calls = _fake_osmnx(monkeypatch)
    rd2dnet.NetworkGraph.from_place(bbox=(-74.0, 40.7, -73.9, 40.8))
    assert calls[0][0] == "bbox"


def test_missing_coordinates_error():
    G = nx.path_graph(2)
    with pytest.raises(KeyError, match="coordinates"):
        rd2dnet.NetworkGraph(G)


@pytest.mark.online
def test_osm_download_midtown():  # pragma: no cover - needs internet
    net = rd2dnet.NetworkGraph.from_place("Midtown Manhattan, New York", offline_fallback=False)
    assert net.n_nodes > 100


def test_flow_direction_measures_both_sides_of_a_river():
    # a river drawn downstream (west to east) crossing a border at x = 200
    G = nx.DiGraph()
    for k in range(5):
        G.add_node(k, x=100.0 * k, y=0.0)
    for k in range(4):
        G.add_edge(k, k + 1, length=100.0)
    G.add_node(9, x=200.0, y=300.0)  # a tributary joining downstream of the border
    G.add_edge(9, 3, length=float(np.hypot(100, 300)))
    bnd = rd2dnet.Boundary([(200, -500), (200, 500)], treated_side="right")  # treated: east (downstream)
    pts = np.array([[50.0, 0.0], [350.0, 0.0]])
    out = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G, direction="out"), bnd, pts, eval_points=[500.0])
    assert np.isnan(out.s[0])  # upstream cannot be reached from the border going downstream
    flow = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G, direction="flow"), bnd, pts, eval_points=[500.0])
    np.testing.assert_allclose(flow.s, [-150.0, 150.0])
    np.testing.assert_allclose(flow.at(0)[0], [150.0, 150.0])
