import warnings

import networkx as nx
import numpy as np
import pandas as pd
import pytest
import shapely
from shapely.geometry import LineString, box

import rd2dnet


def _segments(net):
    E = net.edge_frame()
    return pd.DataFrame({"x1": net.coords[E.u, 0], "y1": net.coords[E.u, 1],
                         "x2": net.coords[E.v, 0], "y2": net.coords[E.v, 1]})


@pytest.fixture(scope="module")
def eval_b(city):
    return city.boundary.interpolate(city.boundary.grid(9))


# --- rd2d-style entry ---------------------------------------------------------

def test_rd2dnet_matches_step_by_step(city, eval_b):
    obs, net = city.obs, city.network()
    res = rd2dnet.rd2dnet(obs["outcome"], obs[["x", "y"]], obs["treated"], eval_b, net,
                          boundary=city.boundary, h=900, seed=0)
    coords = rd2dnet.network_coordinates(net, city.boundary, obs[["x", "y"]], eval_points=eval_b)
    ref = rd2dnet.rdnet(obs["outcome"], coords, h=900, seed=0)
    pd.testing.assert_frame_equal(res.table, ref.table)
    assert res.coords is not None
    assert np.isfinite(rd2dnet.density(res.coords).pvalue)


def test_rd2dnet_infers_the_treated_side(city, eval_b):
    obs = city.obs
    X, t = obs[["x", "y"]].to_numpy(), obs["treated"].to_numpy()
    side = rd2dnet.Boundary.from_points(eval_b, X, t).treated_side
    # reversing b flips the side but describes the same design
    assert rd2dnet.Boundary.from_points(eval_b[::-1], X, t).treated_side != side
    a = rd2dnet.rd2dnet(obs["outcome"], X, t, eval_b, city.network(), h=900, seed=0)
    b = rd2dnet.rd2dnet(obs["outcome"], X, t, eval_b, city.network(), treated_side=side, h=900, seed=0)
    pd.testing.assert_frame_equal(a.table, b.table)
    np.testing.assert_allclose(a.coords.s[t], np.abs(a.coords.s[t]))


def test_rd2dnet_accepts_networkx_and_segments(city, eval_b):
    obs, net = city.obs, city.network()
    args = (obs["outcome"], obs[["x", "y"]], obs["treated"], eval_b)
    ref = rd2dnet.rd2dnet(*args, net, boundary=city.boundary, h=900, seed=0).table
    for alt in (net.graph, _segments(net)):
        got = rd2dnet.rd2dnet(*args, alt, boundary=city.boundary, h=900, seed=0).table
        np.testing.assert_allclose(got["tau_bc"], ref["tau_bc"], rtol=1e-10)


def test_rd2dnet_input_errors(city, eval_b):
    obs, net = city.obs, city.network()
    with pytest.raises(ValueError, match="0/1"):
        rd2dnet.rd2dnet(obs["outcome"], obs[["x", "y"]], np.full(len(obs), 2), eval_b, net)
    with pytest.raises(ValueError, match="same number"):
        rd2dnet.rd2dnet(obs["outcome"], obs[["x", "y"]], obs["treated"][:10], eval_b, net)
    with pytest.raises(ValueError, match="both treated and control"):
        rd2dnet.Boundary.from_points(eval_b, obs[["x", "y"]], np.ones(len(obs)))
    with pytest.raises(ValueError, match="two distinct"):
        rd2dnet.Boundary.from_points([[0, 0], [0, 0]], treated_side="left")


def test_side_inference_warns_when_t_disagrees(city, eval_b):
    rng = np.random.default_rng(0)
    t = rng.integers(0, 2, len(city.obs))
    with pytest.warns(UserWarning, match="side of the boundary"):
        rd2dnet.Boundary.from_points(eval_b, city.obs[["x", "y"]].to_numpy(), t)


# --- networks from lines and tables ----------------------------------------------

def test_from_lines_segment_table_reproduces_network(city):
    net = city.network()
    net2 = rd2dnet.NetworkGraph.from_lines(_segments(net))
    assert net2.n_nodes == net.n_nodes
    assert net2.graph.number_of_edges() == len(net.edge_frame())
    pts = city.obs[["x", "y"]].iloc[:200]
    c1 = rd2dnet.network_coordinates(net, city.boundary, pts)
    c2 = rd2dnet.network_coordinates(net2, city.boundary, pts)
    np.testing.assert_allclose(c1.s, c2.s, atol=1e-8)


def test_from_lines_keeps_attributes_for_metrics():
    seg = pd.DataFrame({"x1": [0, 100], "y1": [0, 0], "x2": [100, 200], "y2": [0, 0], "speed": [10.0, 1.0]})
    net = rd2dnet.NetworkGraph.from_lines(seg, rd2dnet.travel_time(speed_attr="speed"))
    a, c = net.nodes[0], net.nodes[-1]
    assert net.node_distance(a, c) == pytest.approx(100 / (10 / 3.6) + 100 / (1 / 3.6))


def test_from_lines_joins_nearby_ends_and_splits_crossings():
    gap = [[(0, 0), (100, 0)], [(100.4, 0), (200, 0)]]
    assert rd2dnet.NetworkGraph.from_lines(gap).n_nodes == 4
    assert rd2dnet.NetworkGraph.from_lines(gap, tol=1).n_nodes == 3
    cross = [[(0, -500), (0, 500)], [(-500, 0), (500, 0)]]
    assert rd2dnet.NetworkGraph.from_lines(cross).n_nodes == 4
    noded = rd2dnet.NetworkGraph.from_lines(cross, split_at_intersections=True)
    assert noded.n_nodes == 5
    i = noded.nodes[int(np.argmin(np.abs(noded.coords).sum(1)))]
    assert noded.graph.degree(i) == 4


def test_from_lines_curved_geometry_and_rings():
    arc = LineString([(0, 0), (50, 50), (100, 0)])
    net = rd2dnet.NetworkGraph.from_lines([arc])
    assert net.node_distance(*net.nodes) == pytest.approx(arc.length)
    ring = rd2dnet.NetworkGraph.from_lines([[(0, 0), (100, 0), (100, 100), (0, 100), (0, 0)]])
    assert ring.n_nodes == 2 and ring.graph.number_of_edges() == 2
    with pytest.warns(UserWarning, match="zero-length"):
        rd2dnet.NetworkGraph.from_lines([[(0, 0), (0, 0)], [(0, 0), (1, 0)]])


def test_from_lines_directed_river():
    # a river flowing east; distance only downstream with direction="out"
    river = [[(0, 0), (100, 0)], [(100, 0), (200, 0)]]
    net = rd2dnet.NetworkGraph.from_lines(river, directed=True, direction="out")
    up, down = net.nodes[int(np.argmin(net.coords[:, 0]))], net.nodes[int(np.argmax(net.coords[:, 0]))]
    assert net.node_distance(up, down) == pytest.approx(200)
    assert np.isinf(net.node_distance(down, up))


def test_from_lines_geodataframe_lonlat_is_projected():
    gpd = pytest.importorskip("geopandas")
    pytest.importorskip("pyproj")
    lines = gpd.GeoDataFrame({"name": ["a", "b"]}, crs="EPSG:4326",
                             geometry=[LineString([(-73.99, 40.75), (-73.98, 40.75)]),
                                       LineString([(-73.98, 40.75), (-73.98, 40.76)])])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        net = rd2dnet.NetworkGraph.from_lines(lines)
    assert "326" in net.crs  # UTM zone
    total = net.node_distance(net.nodes[0], net.nodes[-1])
    assert 1800 < total < 2000  # about 840 m + 1110 m
    assert net.graph.edges[net.nodes[0], net.nodes[1], 0]["name"] == "a"
    # points in another CRS are reprojected to the network's
    pts = gpd.GeoSeries(gpd.points_from_xy([-73.985], [40.7501]), crs="EPSG:4326")
    bnd = rd2dnet.Boundary([(net.coords[1, 0] - 1, net.coords[1, 1] - 500), (net.coords[1, 0] - 1, net.coords[1, 1] + 2000)],
                           treated_side="left")
    c = rd2dnet.network_coordinates(net, bnd, pts)
    assert np.isfinite(c.s).all()


def test_from_frames():
    nodes = pd.DataFrame({"node": ["a", "b", "c"], "lon": [0.0, 100.0, 100.0], "lat": [0.0, 0.0, 50.0]})
    edges = pd.DataFrame({"u": ["a", "b"], "v": ["b", "c"], "minutes": [3.0, np.nan]})
    net = rd2dnet.NetworkGraph.from_frames(nodes, edges, node_id="node", x="lon", y="lat", source="u", target="v")
    assert net.node_distance("a", "c") == pytest.approx(150)
    assert "minutes" not in net.graph.edges["b", "c", 0]  # NaN attributes are dropped
    with pytest.raises(KeyError, match="missing from nodes"):
        rd2dnet.NetworkGraph.from_frames(nodes, pd.DataFrame({"u": ["a"], "v": ["z"]}), node_id="node",
                                         x="lon", y="lat", source="u", target="v")
    with pytest.raises(KeyError, match="x= and y="):
        rd2dnet.NetworkGraph.from_frames(nodes, edges)


def test_from_region_outline():
    zone = box(1000, -1000, 3000, 1000)
    whole = rd2dnet.Boundary.from_region(zone)
    assert whole.length == pytest.approx(8000)
    part = rd2dnet.Boundary.from_region(zone, within=box(-1e6, -1e6, 1500, 1e6))
    assert part.length == pytest.approx(3000)
    assert part.is_treated([[2000, 0], [0, 0]]).tolist() == [True, False]


def test_lonlat_network_warns():
    G = nx.Graph()
    G.add_node(0, x=-73.98, y=40.75)
    G.add_node(1, x=-73.97, y=40.76)
    G.add_edge(0, 1, length=1200.0)
    with pytest.warns(UserWarning, match="longitude/latitude"):
        rd2dnet.NetworkGraph(G)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rd2dnet.NetworkGraph(nx.relabel_nodes(G, {0: 0, 1: 1}), crs="EPSG:32618")


def test_points_from_two_column_frame(city):
    net = city.network()
    pts = city.obs[["x", "y"]].iloc[:50].rename(columns={"x": "east", "y": "north"})
    c = rd2dnet.network_coordinates(net, city.boundary, pts)
    c0 = rd2dnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]].iloc[:50])
    np.testing.assert_array_equal(c.s, c0.s)
