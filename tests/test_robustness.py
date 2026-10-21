"""Edge cases found by trying the package on many kinds of networks and data."""

import warnings

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, box

import rd2dnet


def _path_network(n=41, spacing=250.0, y=0.0):
    G = nx.path_graph(n)
    for i in G.nodes:
        G.nodes[i]["x"], G.nodes[i]["y"] = (i - n // 2) * spacing, y
    for u, v in G.edges:
        G.edges[u, v]["length"] = spacing
    return G


def _rivers(seed=0, n=3000):
    """Three parallel rivers crossing a border at x = 0: o has no variation near each crossing."""
    rng = np.random.default_rng(seed)
    lines = []
    for yy in (-3000, 0, 3000):
        xs = np.linspace(-6000, 6000, 25)
        lines += [[(xs[i], yy), (xs[i + 1], yy)] for i in range(24)]
    net = rd2dnet.NetworkGraph.from_lines(lines)
    r = rng.integers(0, 3, n)
    x = rng.uniform(-6000, 6000, n)
    y = np.array([-3000, 0, 3000])[r] + rng.normal(0, 30, n)
    t = x > 0
    Y = 1 + 0.0002 * x - 1.5 * t + rng.normal(0, 0.5, n)
    bd = rd2dnet.Boundary([(0, -4000), (0, 4000)], treated_side="right")
    return net, bd, np.c_[x, y], t, Y


def test_few_crossings_fall_back_to_distance_only_fit():
    net, bd, X, t, Y = _rivers()
    b = np.array([[0, -3000], [0, 0], [0, 3000]])
    res = rd2dnet.rd2dnet(Y, X, t, b, net, boundary=bd, seed=0)
    assert (res.table["status"] == "ok").all()
    assert (res.table["basis"] == "s").all()
    assert abs(res.wbate()["tau_bc"] + 1.5) < 0.5
    # the 2-d fit is kept wherever it works
    from rd2dnet.datasets import load_mock_midtown

    city = load_mock_midtown()
    c = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]], n_eval=9)
    r = rd2dnet.rdnet(city.obs["outcome"], c, h=900, seed=0)
    assert set(r.table.loc[r.table.status == "ok", "basis"]) == {"s,o"}


def test_many_eval_points_on_few_crossings_warn():
    net, bd, X, t, Y = _rivers(n=600)
    with pytest.warns(UserWarning, match="distinct places"):
        rd2dnet.network_coordinates(net, bd, X, n_eval=15)


def test_single_cluster_per_side_is_flagged(city, coords):
    cl = np.where(coords.treated, "A", "B")
    r = rd2dnet.rdnet(city.obs["outcome"], coords, h=900, cluster=cl, seed=0)
    assert "few clusters" in set(r.table["status"])
    assert set(r.table["status"]) <= {"few clusters", "masked", "insufficient"}


def test_bandwidth_floor_keeps_points_estimable(city, coords):
    bw = rd2dnet.rdbwnet(city.obs["outcome"], coords)
    ok = bw["status"] == "ok"
    h = np.maximum(bw["hmin0"], bw["hmin1"]).to_numpy()
    h = np.where(ok, h, 1000.0)
    r = rd2dnet.rdnet(city.obs["outcome"], coords, h=h, seed=0)
    assert not (r.table.loc[ok.to_numpy(), "status"] == "insufficient").any()


def test_closed_ring_offsets_wrap_and_grid_does_not_repeat():
    G = nx.grid_2d_graph(21, 21)
    G = nx.convert_node_labels_to_integers(G, label_attribute="ij")
    for _, d in G.nodes(data=True):
        d["x"], d["y"] = 100.0 * d["ij"][0], 100.0 * d["ij"][1]
    for u, v, d in G.edges(data=True):
        d["length"] = 100.0
    zone = box(550, 550, 1450, 1450)
    bnd = rd2dnet.Boundary.from_region(zone)
    assert bnd.closed
    rng = np.random.default_rng(0)
    pts = rng.uniform(0, 2000, (2000, 2))
    c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G), bnd, pts, n_eval=12)
    gaps = np.diff(np.sort(c.eval_arc))
    assert gaps.min() > 0.5 * bnd.length / 12
    # offsets from the point at arc 0 change sign across the ring's start
    j = int(np.argmin(c.eval_arc))
    d, o = c.at(j)
    arc = c.crossing_arc
    L = bnd.length
    near = np.isfinite(o) & (d < 400)
    before = near & (arc > L / 2)   # just before the seam going round
    after = near & (arc < L / 2) & (arc > c.eval_node_arc[j])
    assert (o[before] <= 0).all() and (o[after] >= 0).all()


def test_density_keeps_controls_on_boundary_and_flags_repeats(city, coords):
    with pytest.warns(UserWarning, match="share a location"):
        idx = np.repeat(np.arange(len(city.obs)), 3)
        c3 = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]].iloc[idx], n_eval=5)
        rd2dnet.density(c3)
    d1 = rd2dnet.density(coords)
    d3 = rd2dnet.density(c3, unique_locations=True)
    assert d3.pvalue == pytest.approx(d1.pvalue, rel=1e-6)


def test_split_edges_divide_the_metric_attribute():
    from scipy.sparse import csr_matrix

    net = rd2dnet.NetworkGraph.from_adjacency(csr_matrix([[0, 200.0], [200.0, 0]]), [[-100, 0], [100, 0]])
    c = rd2dnet.network_coordinates(net, rd2dnet.Boundary([(0, -50), (0, 50)], treated_side="right"),
                                    [[-100, 0], [100, 0]], n_eval=1)
    np.testing.assert_allclose(c.s, [-100, 100])


def test_parallel_edges_use_the_cheaper_one():
    out = []
    for order in ((300.0, 100.0), (100.0, 300.0)):
        M = nx.MultiGraph()
        M.add_node("a", x=0.0, y=-100.0)
        M.add_node("b", x=0.0, y=100.0)
        for L in order:
            M.add_edge("a", "b", length=L)
        c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(M), rd2dnet.Boundary([(-50, 0), (50, 0)], treated_side="left"),
                                        [[0, 50], [0, -50]], eval_points=[50.0])
        out.append(c.s)
    np.testing.assert_allclose(out[0], out[1])


def test_input_validation(city, coords):
    y = city.obs["outcome"]
    with pytest.raises(ValueError, match="percentage"):
        rd2dnet.rdnet(y, coords, h=900, level=0.95)
    w = np.ones(len(y))
    w[3] = np.inf
    with pytest.raises(ValueError, match="finite"):
        rd2dnet.rdnet(y, coords, h=900, weights=w)
    with pytest.raises(ValueError, match="one value per observation"):
        rd2dnet.rdnet(y, coords, h=900, cluster=np.arange(5))
    with pytest.raises(ValueError, match="missing"):
        rd2dnet.rdnet(y, coords, h=900, cluster=np.where(np.arange(len(y)) < 5, np.nan, 1.0))
    with pytest.raises(ValueError, match="positive"):
        rd2dnet.rdnet(y, coords, h=-1)
    with pytest.raises(ValueError, match="vce"):
        rd2dnet.rdbwnet(y, coords, vce="bogus")
    with pytest.raises(ValueError, match="boolean or 0/1"):
        rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]],
                                    treated=np.where(city.obs["treated"], "yes", "no"))
    with pytest.warns(UserWarning, match="Conley"):
        rd2dnet.rdnet(y, coords, h=900, cluster=np.arange(len(y)) % 50, conley_cutoff=100, seed=0)


def test_h_pair_with_two_eval_points(city):
    c2 = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]], eval_points=[600.0, 2600.0])
    r = rd2dnet.rdnet(city.obs["outcome"], c2, h=(800, 1000), seed=0)
    assert (r.table["h0"] == 800).all() and (r.table["h1"] == 1000).all()
    r2 = rd2dnet.rdnet(city.obs["outcome"], c2, h=np.array([800, 1000]), seed=0)
    assert r2.table["h0"].tolist() == [800, 1000]


def test_status_reports_why_selection_failed():
    G = _path_network()
    bd = rd2dnet.Boundary([(0, -500), (0, 500)], treated_side="right")
    rng = np.random.default_rng(0)
    x = np.r_[rng.uniform(-5000, 0, 300), rng.uniform(0, 5000, 10)]
    c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G), bd, np.c_[x, np.zeros_like(x) + 5], eval_points=[500.0])
    r = rd2dnet.rdnet(rng.normal(size=len(x)), c)
    assert r.table["status"].iloc[0] == "insufficient"


def test_wbate_defaults_to_equal_weights(city, coords):
    r = rd2dnet.rdnet(city.obs["outcome"], coords, h=900, seed=0)
    pd.testing.assert_series_equal(r.wbate(), r.wbate(r.table["status"].eq("ok").astype(float)))


def test_travel_time_default_covers_nan_speeds():
    seg = pd.DataFrame({"x1": [0, 100], "y1": [0, 0], "x2": [100, 200], "y2": [0, 0], "speed_kph": [36.0, np.nan]})
    net = rd2dnet.NetworkGraph.from_lines(seg, rd2dnet.travel_time(default_kph=18))
    assert net.node_distance(net.nodes[0], net.nodes[-1]) == pytest.approx(10 + 20)


def test_missing_metric_attribute_fails_early():
    G = nx.Graph()
    G.add_node(0, x=0, y=0)
    G.add_node(1, x=1, y=0)
    G.add_edge(0, 1, len_m=1.0)
    with pytest.raises(KeyError, match="len_m"):
        rd2dnet.NetworkGraph(G)


def test_boundary_accepts_multiline_and_rejects_polygon():
    from shapely.geometry import MultiLineString

    b = rd2dnet.Boundary(MultiLineString([[(0, 0), (1, 0)], [(1, 0), (2, 0)]]), treated_side="left")
    assert b.length == pytest.approx(2)
    with pytest.raises(TypeError, match="from_region"):
        rd2dnet.Boundary(box(0, 0, 1, 1), treated_region=box(0, 0, 1, 1))


def test_near_miss_tolerance_warns():
    G = _path_network(n=3, spacing=100.0)
    G = nx.Graph(G)
    G.add_node("top", x=0.0, y=100.0)
    G.add_node("bottom", x=0.0, y=-100.0)
    G.add_edge("top", 1, length=100.0)
    G.add_edge("bottom", 1, length=100.0)
    G.add_node("t2", x=-100.0, y=100.0)
    G.add_edge("t2", 0, length=100.0)
    G.add_node("t3", x=100.0, y=100.0)
    G.add_edge("t3", 2, length=100.0)
    # the boundary runs 1 cm off the street through nodes 0, 1, 2
    bd = rd2dnet.Boundary([(-500, 0.01), (500, 0.01)], treated_side="left")
    with pytest.warns(UserWarning, match="increase tol"):
        rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G), bd, [[0, 50], [0, -50]], n_eval=3)


def test_inputs_that_used_to_pass_silently(coords, city):
    y = city.obs["outcome"]
    with pytest.warns(UserWarning, match="matched to the observations by position"):
        rd2dnet.rdnet(y.sample(frac=1, random_state=0), coords, h=600.0, cbands=False)
    with pytest.raises(ValueError, match="one-dimensional"):
        rd2dnet.rdnet(np.c_[y, y][: coords.n_obs // 2], coords, h=600.0)
    with pytest.raises(ValueError, match="finite and non-negative"):
        rd2dnet.density(coords, weights=-np.ones(coords.n_obs))
    with pytest.raises(ValueError, match="zero length"):
        rd2dnet.Boundary([(500, 1000), (500, 1000)], treated_side="left")
    with pytest.raises(ValueError, match="two"):
        rd2dnet.Boundary([(500, 1000)], treated_side="left")
    net = city.network("length")
    with pytest.raises(ValueError, match="finite"):
        rd2dnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]], eval_points=[np.nan, 500.0])
    with pytest.raises(ValueError, match="positive integer"):
        rd2dnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]], n_eval=0)
    with pytest.raises(ValueError, match="snap_scale"):
        rd2dnet.NetworkGraph(net.graph, "length", snap_scale=-5)


def test_collinear_covariates_warn(coords, city):
    z = city.obs["far"].to_numpy()
    with pytest.warns(UserWarning, match="collinear"):
        rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0, covs=np.c_[z, 2 * z + 1], cbands=False)


def test_clear_errors_for_malformed_inputs():
    from rd2dnet import Precomputed

    G = nx.Graph()
    G.add_node(0, x=0.0, y=np.nan)
    G.add_node(1, x=1.0, y=0.0)
    G.add_edge(0, 1, length=1.0)
    with pytest.raises(ValueError, match="finite"):
        rd2dnet.NetworkGraph(G)
    with pytest.raises(ValueError, match="arc_start, arc_end"):
        rd2dnet.Boundary([(0, 0), (1, 0)], treated_side="left", mask=[5])
    with pytest.raises(ValueError, match="one row per node"):
        rd2dnet.NetworkGraph.from_adjacency(np.ones((3, 3)) - np.eye(3), np.zeros((2, 2)))
    G2 = nx.Graph()
    G2.add_node("a", x=0.0, y=0.0)
    G2.add_node("b", x=1.0, y=0.0)
    G2.add_edge("a", "b", length=1.0)
    with pytest.raises(ValueError, match="not in the precomputed"):
        rd2dnet.NetworkGraph(G2, Precomputed(np.zeros((2, 2)), ["a", "c"])).engine


@pytest.mark.parametrize("kw, msg", [
    (dict(deriv=(0.5, 0)), "deriv must be an integer"),
    (dict(p=1.5), "p must be an integer"),
    (dict(q=2.5), "q must be an integer"),
    (dict(repp=0), "repp must be at least 1"),
    (dict(tangvec=(1.0, 0.0, 0.0)), "tangvec must be"),
    (dict(tangvec=(np.nan, 1.0)), "tangvec must be finite"),
    (dict(tangvec=(0.0, 0.0)), r"tangvec must not be \(0, 0\)"),
])
def test_rdnet_rejects_malformed_options(city, coords, kw, msg):
    with pytest.raises(ValueError, match=msg):
        rd2dnet.rdnet(city.obs["outcome"], coords, h=600, **kw)


def test_rdbwnet_rejects_fractional_deriv(city, coords):
    with pytest.raises(ValueError, match="deriv must be an integer"):
        rd2dnet.rdbwnet(city.obs["outcome"], coords, deriv=(0.5, 0))


def test_manipulation_test_rejects_bad_q_and_k_max(coords):
    with pytest.raises(ValueError, match="k_max must be at least 1"):
        rd2dnet.manipulation_test(coords, k_max=0)
    with pytest.raises(ValueError, match="q must be an integer"):
        rd2dnet.manipulation_test(coords, q=50.5)


def test_rd2dnet_with_missing_locations_says_so(city):
    X = city.obs[["x", "y"]].to_numpy().copy()
    X[3, 0] = np.nan
    b = city.boundary.interpolate(city.boundary.grid(5, start=100, end=2900))
    with pytest.raises(ValueError, match="X has missing"):
        rd2dnet.rd2dnet(city.obs["outcome"], X, city.obs["treated"], b, city.network())


def test_conley_covariance_is_made_positive_semidefinite():
    from rd2dnet.estimate import _psd

    V = np.array([[1.0, 2.0], [2.0, 1.0]])  # eigenvalues 3 and -1
    with pytest.warns(UserWarning, match="not positive semi-definite"):
        P = _psd(V)
    assert np.linalg.eigvalsh(P).min() >= -1e-12
    np.testing.assert_allclose(P, 1.5 * np.ones((2, 2)))
    W = np.array([[2.0, 1.0], [1.0, 2.0]])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _psd(W) is W
