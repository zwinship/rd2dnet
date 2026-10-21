import numpy as np
import pytest

import rd2dnet
from rd2dnet.datasets import make_mock_midtown


def test_density_no_sorting(coords):
    r = rd2dnet.density(coords)
    assert r.normalize == "euclidean"
    assert np.isfinite(r.diff_bc) and r.se > 0
    assert r.pvalue > 0.01


def test_density_network_normalization(coords):
    r = rd2dnet.density(coords, normalize="network")
    assert r.normalize == "network" and np.isfinite(r.diff_bc)
    assert r.pvalue > 0.01
    assert abs(r.theta) < 0.5  # uniform parcels per metre of street on both sides


def test_jackknife_matches_simulation():
    """The jackknife SE tracks the Monte Carlo spread of the density estimate."""
    from rd2dnet.density import _cjm
    from rd2dnet.kernels import triangular

    rng = np.random.default_rng(0)
    est, ses = [], []
    for _ in range(200):
        x = rng.uniform(-1, 1, 2000)
        b, V, _, _ = _cjm(x, np.ones_like(x), 0.4, 0.4, 2, triangular)
        est.append(b[3])
        ses.append(np.sqrt(V[3, 3]))
    assert np.mean(est) == pytest.approx(0.5, abs=0.02)
    assert np.mean(ses) == pytest.approx(np.std(est), rel=0.2)


def test_density_detects_bunching():
    city = make_mock_midtown(n_obs=3000, seed=5)
    obs = city.obs
    rng = np.random.default_rng(0)
    extra = obs.sample(1500, random_state=1)[["x", "y"]].to_numpy().copy()
    extra[:, 1] = 1600 - rng.uniform(1, 60, len(extra))  # pile up just inside the treated zone
    pts = np.vstack([obs[["x", "y"]].to_numpy(), extra])
    c = rd2dnet.network_coordinates(city.network(), city.boundary, pts, n_eval=5)
    for norm in ("euclidean", "network", "none"):
        r = rd2dnet.density(c, normalize=norm)
        assert r.theta > 0 and r.pvalue < 0.01


def test_density_excludes_masked_stretches(coords):
    r = rd2dnet.density(coords, normalize="network")
    masked = coords.boundary.is_masked(coords.crossing_arc)
    assert masked.any()
    assert r.n_left + r.n_right == int((~masked & np.isfinite(coords.s)).sum())
    # straight-line version: masked by where a location projects onto the boundary
    b = coords.boundary
    arc = b.project(coords.xy)
    keep = ~b.is_masked(arc) & (arc > 1e-9 * b.length) & (arc < b.length * (1 - 1e-9))
    e = rd2dnet.density(coords)
    assert e.n_left + e.n_right == int(keep.sum())


def test_density_euclidean_is_valid_on_a_round_zone():
    # parcels uniform in the plane around a round zone on a street grid: no sorting
    import networkx as nx
    from shapely.geometry import Point

    G = nx.grid_2d_graph(31, 31)
    G = nx.convert_node_labels_to_integers(G, label_attribute="ij")
    for _, d in G.nodes(data=True):
        d["x"], d["y"] = 100.0 * d["ij"][0], 100.0 * d["ij"][1]
    for u, v, d in G.edges(data=True):
        d["length"] = 100.0
    net = rd2dnet.NetworkGraph(G)
    zone = rd2dnet.Boundary.from_region(Point(1500, 1500).buffer(800, 64))
    pvals = []
    for seed in range(5):
        pts = np.random.default_rng(seed).uniform(0, 3000, (12000, 2))
        c = rd2dnet.network_coordinates(net, zone, pts, n_eval=8)
        pvals.append(rd2dnet.density(c).pvalue)
    assert np.mean(np.array(pvals) < 0.05) <= 0.4  # the street-length version rejects ~every time


def test_density_bwselect(coords):
    a = rd2dnet.density(coords, bwselect="common")
    b = rd2dnet.density(coords, bwselect="each")
    assert a.h_left == a.h_right
    assert b.h_left != b.h_right
    with pytest.raises(ValueError):
        rd2dnet.density(coords, bwselect="bogus")


def test_density_local_window(coords):
    j = int(np.flatnonzero(~coords.eval_masked)[2])
    r = rd2dnet.density(coords, eval_index=j, window=600.0, h=400.0)
    assert r.n_left + r.n_right < coords.n_obs
    with pytest.raises(ValueError):
        rd2dnet.density(coords, eval_index=j)
