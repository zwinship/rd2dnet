import numpy as np
import pytest
from shapely.geometry import LineString, box

import rdnet


def test_treated_side_left_right():
    line = LineString([(0, 0), (10, 0)])
    left = rdnet.Boundary(line, treated_side="left")
    right = rdnet.Boundary(line, treated_side="right")
    pts = [[5, 1], [5, -1]]
    np.testing.assert_array_equal(left.is_treated(pts), [True, False])
    np.testing.assert_array_equal(right.is_treated(pts), [False, True])


def test_boundary_requires_one_side_spec():
    with pytest.raises(ValueError):
        rdnet.Boundary([(0, 0), (1, 0)])


def test_mask_geometry_and_arcs():
    line = LineString([(0, 0), (100, 0)])
    b = rdnet.Boundary(line, treated_side="left", mask=[box(20, 0, 40, 10), (80, 90)])
    np.testing.assert_array_equal(b.is_masked([10, 30, 50, 85]), [False, True, False, True])


def test_grid_coordinates_signs_and_values(grid):
    G, bnd = grid
    net = rdnet.NetworkGraph(G)
    pts = np.array([[100.0, 50.0], [100.0, 350.0], [300.0, 120.0]])
    c = rdnet.network_coordinates(net, bnd, pts, eval_points=[150.0])  # eval at x = 100
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
    c = rdnet.network_coordinates(rdnet.NetworkGraph(G), bnd, [[100.0, 50.0]], eval_points=[200.0])  # x = 150
    d, o = c.at(0)
    assert c.eval_node_arc[0] == pytest.approx(200.0)
    assert d[0] == pytest.approx(50.0 + 150.0)
    assert o[0] == pytest.approx(-50.0)


def test_metric_changes_coordinates(grid):
    G, bnd = grid
    for u, v, d in G.edges(data=True):
        d["travel_time"] = d["length"] / 10.0
    net = rdnet.NetworkGraph(G, "travel_time")
    c = rdnet.network_coordinates(net, bnd, [[100.0, 0.0]])
    assert c.s[0] == pytest.approx(20.0)  # 200 m at 10 m/s, snapping depth costs 0


def test_unreachable_points_dropped(grid):
    G, bnd = grid
    G.add_node("island", x=1000.0, y=0.0)
    G.add_node("island2", x=1000.0, y=50.0)
    G.add_edge("island", "island2", length=50.0)
    with pytest.warns(UserWarning, match="cannot reach"):
        c = rdnet.network_coordinates(rdnet.NetworkGraph(G), bnd, [[1000.0, 10.0], [0.0, 0.0]])
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
    c = rdnet.network_coordinates(rdnet.NetworkGraph(G), bnd, [[120.0, 210.0], [120.0, 190.0]], eval_points=[200.0])
    d, o = c.at(0)  # evaluation point at x = 150 on the boundary street
    np.testing.assert_allclose(d, [10.0 + 30.0, 10.0 + 30.0])
    np.testing.assert_allclose(o, [-30.0, -30.0])
    np.testing.assert_allclose(np.abs(c.s), [10.0, 10.0])
