import warnings

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, box

import rd2dnet


def test_treated_side_left_right():
    line = LineString([(0, 0), (10, 0)])
    left = rd2dnet.Boundary(line, treated_side="left")
    right = rd2dnet.Boundary(line, treated_side="right")
    pts = [[5, 1], [5, -1]]
    np.testing.assert_array_equal(left.is_treated(pts), [True, False])
    np.testing.assert_array_equal(right.is_treated(pts), [False, True])


def test_boundary_requires_one_side_spec():
    with pytest.raises(ValueError):
        rd2dnet.Boundary([(0, 0), (1, 0)])


def test_mask_geometry_and_arcs():
    line = LineString([(0, 0), (100, 0)])
    b = rd2dnet.Boundary(line, treated_side="left", mask=[box(20, 0, 40, 10), (80, 90)])
    np.testing.assert_array_equal(b.is_masked([10, 30, 50, 85]), [False, True, False, True])


def test_grid_coordinates_signs_and_values(grid):
    G, bnd = grid
    net = rd2dnet.NetworkGraph(G)
    pts = np.array([[100.0, 50.0], [100.0, 350.0], [300.0, 120.0]])
    c = rd2dnet.network_coordinates(net, bnd, pts, eval_points=[150.0])  # eval at x = 100
    np.testing.assert_array_equal(c.treated, [True, False, True])
    # Points on the avenue x=100: straight up/down the avenue to y=200.
    assert c.s[0] == pytest.approx(150.0)
    assert c.s[1] == pytest.approx(-150.0)
    # Point at (300, 120) snaps onto the x=300 avenue: 80 m to the boundary, 0 depth.
    assert c.s[2] == pytest.approx(80.0)
    d, o = c.at(0)
    assert o[0] == pytest.approx(0.0) and o[1] == pytest.approx(0.0)
    assert o[2] == pytest.approx(200.0)  # crossing at x=300 is 200 m along the boundary
    assert d[2] == pytest.approx(280.0)


def test_eval_points_between_nodes(grid):
    G, bnd = grid
    c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G), bnd, [[100.0, 50.0]], eval_points=[200.0])  # x = 150
    d, o = c.at(0)
    assert c.eval_node_arc[0] == pytest.approx(200.0)
    assert d[0] == pytest.approx(50.0 + 150.0)
    assert o[0] == pytest.approx(-50.0)


def test_metric_changes_coordinates(grid):
    G, bnd = grid
    for u, v, d in G.edges(data=True):
        d["travel_time"] = d["length"] / 10.0
    net = rd2dnet.NetworkGraph(G, "travel_time")
    c = rd2dnet.network_coordinates(net, bnd, [[100.0, 0.0]])
    assert c.s[0] == pytest.approx(20.0)  # 200 m at 10 m/s, snapping depth costs 0


def test_unreachable_points_dropped(grid):
    G, bnd = grid
    G.add_node("island", x=1000.0, y=0.0)
    G.add_node("island2", x=1000.0, y=50.0)
    G.add_edge("island", "island2", length=50.0)
    with pytest.warns(UserWarning, match="no network path"):
        c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G), bnd, [[1000.0, 10.0], [0.0, 0.0]])
    assert np.isnan(c.s[0]) and np.isfinite(c.s[1])


def test_mock_city_coordinates(coords, city):
    assert coords.n_obs == len(city.obs)
    assert np.all((coords.s >= 0) == coords.treated)
    # On a grid with Euclidean edge snapping, |s| is at least the straight-line distance.
    straight = np.abs(city.obs["y"].to_numpy() - 1600.0)
    assert np.all(np.abs(coords.s) >= straight - 1e-6)
    assert coords.eval_masked.sum() > 0


def test_pairwise_distances_symmetric(coords):
    idx = np.arange(25)
    D = coords.pairwise(idx)
    np.testing.assert_allclose(D, D.T)
    assert np.all(np.diag(D) == 0)
    assert np.all(D >= 0)


def test_point_on_same_segment_as_eval_point(grid):
    """Distance along a shared segment is direct, not routed out to an endpoint and back."""
    G, bnd = grid
    c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G), bnd, [[120.0, 210.0], [120.0, 190.0]], eval_points=[200.0])
    d, o = c.at(0)  # evaluation point at x = 150 on the boundary street
    np.testing.assert_allclose(d, [10.0 + 30.0, 10.0 + 30.0])
    np.testing.assert_allclose(o, [-30.0, -30.0])
    np.testing.assert_allclose(np.abs(c.s), [10.0, 10.0])


def _dense_hac(c, idx, Psi, cutoff):
    """Reference: dense Bartlett kernel; co-located observations are at distance 0."""
    Dm = c.pairwise(idx, cutoff=cutoff)
    Dm[c.loc[idx][:, None] == c.loc[idx][None, :]] = 0.0
    K = np.clip(1 - Dm / cutoff, 0, None)
    out = Psi.T @ K @ Psi
    return (out + out.T) / 2


def test_repeated_locations_share_coordinates(city):
    xy = city.obs[["x", "y"]].to_numpy()[:400]
    rep = np.repeat(xy, 3, axis=0)[np.random.default_rng(0).permutation(1200)]
    a = rd2dnet.network_coordinates(city.network(), city.boundary, xy, n_eval=9)
    b = rd2dnet.network_coordinates(city.network(), city.boundary, rep, n_eval=9)
    assert b.n_locations == 400 and b.n_obs == 1200
    key = {tuple(p): i for i, p in enumerate(xy)}
    m = np.array([key[tuple(p)] for p in rep])
    np.testing.assert_array_equal(b.s, a.s[m])
    for j in range(9):
        da, oa = a.at(j)
        db, ob = b.at(j)
        np.testing.assert_array_equal(db, da[m])
        np.testing.assert_array_equal(ob, oa[m])


def test_sparse_hac_matches_dense_with_repeats(city):
    xy = np.repeat(city.obs[["x", "y"]].to_numpy()[:250], 2, axis=0)
    c = rd2dnet.network_coordinates(city.network(), city.boundary, xy, n_eval=5)
    idx = np.arange(0, 500, 1)
    Psi = np.random.default_rng(1).normal(size=(500, 3))
    for cutoff in (120.0, 400.0):
        ref = _dense_hac(c, idx, Psi, cutoff)
        np.testing.assert_allclose(c.hac_quadratic(idx, Psi, cutoff), ref, rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(c.hac_quadratic(idx, Psi, cutoff, max_cells=3 * 500), ref, rtol=1e-10, atol=1e-10)


def test_sparse_hac_directed_network_and_euclidean(grid):
    import networkx as nx

    G, bnd = grid
    D = nx.DiGraph(G)
    D.remove_edges_from([(u, v) for u, v in list(D.edges) if u > v and (u + v) % 3 == 0])  # some one-way streets
    rng = np.random.default_rng(2)
    pts = rng.uniform(0, 400, (120, 2))
    c = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(D, direction="out"), bnd, pts, n_eval=3)
    idx = np.flatnonzero(np.isfinite(c.s))
    Psi = rng.normal(size=(len(idx), 2))
    np.testing.assert_allclose(c.hac_quadratic(idx, Psi, 150.0), _dense_hac(c, idx, Psi, 150.0), rtol=1e-10, atol=1e-10)
    diff = c.xy[idx][:, None, :] - c.xy[idx][None, :, :]
    K = np.clip(1 - np.sqrt((diff**2).sum(-1)) / 150.0, 0, None)
    np.testing.assert_allclose(c.hac_quadratic(idx, Psi, 150.0, distance="euclidean"), Psi.T @ K @ Psi, rtol=1e-10)


def test_threaded_geometry_matches_serial(city, monkeypatch):
    from rd2dnet import _parallel

    xy = city.obs[["x", "y"]].to_numpy()
    monkeypatch.setenv("RD2DNET_NUM_THREADS", "1")
    a = rd2dnet.network_coordinates(city.network(), city.boundary, xy, n_eval=7)
    monkeypatch.setattr(_parallel, "MIN_CHUNK", 50)
    monkeypatch.setenv("RD2DNET_NUM_THREADS", "8")
    assert _parallel.n_threads(len(xy)) == 8
    b = rd2dnet.network_coordinates(city.network(), city.boundary, xy, n_eval=7)
    for attr in ("s", "treated", "crossing_arc"):
        np.testing.assert_array_equal(getattr(a, attr), getattr(b, attr))
    np.testing.assert_array_equal(a.snap.cost_u, b.snap.cost_u)


def test_grid_on_a_closed_ring_does_not_repeat_the_start():
    ring = rd2dnet.Boundary.from_region(box(0, 0, 100, 100))
    g = ring.grid(4)
    np.testing.assert_allclose(g, [0, 100, 200, 300])
    assert ring.grid(spacing=100).tolist() == [0, 100, 200, 300]


def test_side_inference_ignores_points_beyond_the_ends_of_b():
    rng = np.random.default_rng(0)
    X = rng.uniform(-1000, 1000, (2000, 2))
    t = X[:, 1] < 0
    # b covers only the middle fifth of the border; the line is drawn through it
    b = np.array([[-200.0, 0.0], [200.0, 0.0]])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        bnd = rd2dnet.Boundary.from_points(b, X, t)
    assert bnd.is_treated(np.array([[0.0, -10.0]]))[0]


def test_lonlat_eval_point_coordinates_are_projected():
    pytest.importorskip("pyproj")
    from rd2dnet import project_lonlat

    lon0, lat0 = -71.2, 42.37
    lon = lon0 + np.linspace(-0.01, 0.01, 11)
    lat = lat0 + np.linspace(-0.01, 0.01, 11)
    G = nx.Graph()
    crs = "EPSG:32619"
    for i, a in enumerate(lon):
        for j, c in enumerate(lat):
            x, y = project_lonlat([a], [c], crs)[0]
            G.add_node((i, j), x=x, y=y)
    for i in range(11):
        for j in range(11):
            for di, dj in ((1, 0), (0, 1)):
                if i + di < 11 and j + dj < 11:
                    p, q = G.nodes[(i, j)], G.nodes[(i + di, j + dj)]
                    G.add_edge((i, j), (i + di, j + dj), length=float(np.hypot(p["x"] - q["x"], p["y"] - q["y"])))
    net = rd2dnet.NetworkGraph(G, crs=crs)
    line = project_lonlat([lon0 - 0.02, lon0 + 0.02], [lat0 + 0.0005, lat0 + 0.0005], crs)
    bnd = rd2dnet.Boundary(line, treated_side="left")
    pts = np.column_stack([lon0 + np.linspace(-0.009, 0.009, 50), lat0 + np.tile([-0.004, 0.004], 25)])
    ep = np.array([[lon0 - 0.005, lat0 + 0.0005], [lon0 + 0.005, lat0 + 0.0005]])
    c = rd2dnet.network_coordinates(net, bnd, pts, lonlat=True, eval_points=ep)
    want = bnd.project(project_lonlat(ep[:, 0], ep[:, 1], crs))
    np.testing.assert_allclose(c.eval_arc, want)


def test_from_lines_projects_plain_lonlat_segments():
    pytest.importorskip("pyproj")
    seg = pd.DataFrame({"x1": [-71.20, -71.19], "y1": [42.37, 42.37], "x2": [-71.19, -71.18], "y2": [42.37, 42.37]})
    net = rd2dnet.NetworkGraph.from_lines(seg, crs="EPSG:4326")
    assert "326" in net.crs  # a UTM zone
    assert 1500 < net.node_distance(net.nodes[0], net.nodes[-1]) < 1800  # about 1.6 km, not degrees


def test_treatment_that_contradicts_geography_warns(city):
    t = city.obs["treated"].to_numpy().astype(bool).copy()
    t[:30] = ~t[:30]
    with pytest.warns(UserWarning, match="30 observations"):
        rd2dnet.network_coordinates(city.network("length"), city.boundary, city.obs[["x", "y"]], treated=t, n_eval=5)
