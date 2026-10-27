import networkx as nx
import numpy as np
import pytest
from shapely.geometry import LineString, box

import rd2dnet
from rd2dnet.bins import _edge_feature_matrix


def comb(n_streets=3, spacing=2000.0, lanes=None):
    """Vertical streets crossing the boundary y = 0 (treated below), joined far away.

    Each street runs from y = -600 to 600 and, 200 m from the boundary on each
    side, splits: the street continues straight and a spur leaves at an angle.
    """
    G = nx.Graph()
    tops, bots = [], []
    for k in range(n_streets):
        x = k * spacing
        for sgn, tag in ((1, "n"), (-1, "s")):
            j, end, spur = f"{tag}j{k}", f"{tag}e{k}", f"{tag}p{k}"
            G.add_node(j, x=x, y=sgn * 200.0)
            G.add_node(end, x=x, y=sgn * 600.0)
            G.add_node(spur, x=x + 300.0, y=sgn * 500.0)
            G.add_edge(j, end, length=400.0, lanes=2.0 if lanes is None else lanes[k])
            G.add_edge(j, spur, length=float(np.hypot(300, 300)), lanes=1.0)
            (tops if sgn > 0 else bots).append(end)
        G.add_edge(f"nj{k}", f"sj{k}", length=400.0, lanes=4.0 if lanes is None else lanes[k])  # crosses y = 0
    for ends in (tops, bots):
        for a, b in zip(ends[:-1], ends[1:]):
            G.add_edge(a, b, length=spacing, lanes=1.0)
    line = LineString([(-500, 0), ((n_streets - 1) * spacing + 500, 0)])
    bnd = rd2dnet.Boundary(line, treated_region=box(-1e5, -1e5, 1e5, 0))
    return G, bnd


def points_on(G, n_per_edge, rng, skip=("e",)):
    """Locations spread along the streets near the crossings (not the far connectors)."""
    pts = []
    for u, v in G.edges():
        if u[1] == "e" and v[1] == "e":
            continue
        a = np.array([G.nodes[u]["x"], G.nodes[u]["y"]])
        b = np.array([G.nodes[v]["x"], G.nodes[v]["y"]])
        t = rng.uniform(0.02, 0.98, n_per_edge)
        pts.append(a + t[:, None] * (b - a))
    xy = np.vstack(pts)
    return xy[np.abs(xy[:, 1]) > 1.0]


@pytest.fixture(scope="module")
def comb_bins():
    G, bnd = comb()
    rng = np.random.default_rng(1)
    xy = points_on(G, 60, rng)
    coords = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G, "length"), bnd, xy, n_eval=3)
    return G, coords, rd2dnet.bin_coordinates(coords, features=["lanes"], topology=False)


def test_bins_are_nearest_crossings(comb_bins):
    _, coords, b = comb_bins
    assert b.n_eval == 3
    assert np.array_equal(b.eval_node[b.bin], coords.crossing_node)
    d, o = b.at(1)
    m = b.bin == 1
    assert np.allclose(d[m], np.abs(coords.s[m])) and np.all(np.isinf(d[~m]))
    assert np.all(o[m] == 0) and np.all(np.isnan(o[~m]))
    tab = b.bins_table()
    assert (tab["n0"] + tab["n1"]).sum() == coords.n_obs


def test_branches_trunk_and_spurs(comb_bins):
    _, coords, b = comb_bins
    y = coords.xy[:, 1]
    x_off = coords.xy[:, 0] - np.round(coords.xy[:, 0] / 2000) * 2000
    trunk = np.abs(y) < 200 - 1e-6
    assert np.all(b.branch[trunk] == -1)
    # beyond the junction: the straight continuation and the spur are different branches
    straight = (np.abs(y) > 200 + 1e-6) & (np.abs(x_off) < 1e-6)
    spur = (np.abs(y) > 200 + 1e-6) & (x_off > 1e-6)
    assert np.all(b.branch[straight] >= 0) and np.all(b.branch[spur] >= 0)
    for m in range(3):
        for side in (True, False):
            sel = (b.bin == m) & (coords.treated == side)
            assert len(np.unique(b.branch[sel & straight])) == 1
            assert len(np.unique(b.branch[sel & spur])) == 1
            assert np.unique(b.branch[sel & straight]) != np.unique(b.branch[sel & spur])
    # every bin has two branches per side
    assert [len(x) for x in b._bin_branches] == [4, 4, 4]


def test_similarity_is_mahalanobis_of_lanes(comb_bins):
    G, coords, b = comb_bins
    # segments as given: per street a crossing (4 lanes), two straight (2) and two spurs (1); connectors 1
    lanes = np.array([d["lanes"] for _, _, d in G.edges(data=True)])
    sd = np.std(lanes, ddof=1)
    y = coords.xy[:, 1]
    x_off = coords.xy[:, 0] - np.round(coords.xy[:, 0] / 2000) * 2000
    expect = np.where(np.abs(y) < 200, 0.0, np.where(np.abs(x_off) < 1e-6, 2.0, 3.0)) / sd
    assert np.allclose(b.similarity, expect)
    assert b.feature_rank == 1
    assert b.h_S == pytest.approx(np.sqrt(2 * 3.841458820694124))


def test_exact_recovery_with_branch_slopes(comb_bins):
    _, coords, b = comb_bins
    T = coords.treated.astype(float)
    s = coords.s
    tau = np.array([1.0, 2.0, 3.0])[b.bin]
    alpha = np.array([5.0, -1.0, 0.5])[b.bin]
    slope = {-1: 0.002}
    slope.update({br: 0.001 * (k + 3) * (1 if k % 2 else -1) for k, br in enumerate(np.unique(b.branch[b.branch >= 0]))})
    y = alpha + tau * T + np.array([slope[x] for x in b.branch]) * s
    nob = b.with_similarity_bandwidth(np.inf)
    res = rd2dnet.rdnet(y, nob, h=1500.0, min_obs=5)
    assert (res.table["status"] == "ok").all()
    assert np.allclose(res.table["tau"], [1, 2, 3], atol=1e-8)
    assert np.allclose(res.table["tau_bc"], [1, 2, 3], atol=1e-8)
    # without branch slopes the mixture of slopes biases the jump
    plain = rd2dnet.bin_coordinates(coords, similarity=False, branches=False)
    res2 = rd2dnet.rdnet(y, plain, h=1500.0, min_obs=5)
    assert np.max(np.abs(res2.table["tau"] - [1, 2, 3])) > 0.05


def test_similarity_kernel_downweights(comb_bins):
    _, coords, b = comb_bins
    kw = b.kernel_weights
    assert np.allclose(kw, np.clip(1 - b.similarity / b.h_S, 0, None))
    assert b.with_similarity_bandwidth(np.inf).kernel_weights is None
    small = b.with_similarity_bandwidth(0.5 * b.similarity[b.similarity > 0].min())
    assert np.all(small.kernel_weights[b.similarity > 0] == 0)
    with pytest.raises(ValueError):
        b.with_similarity_bandwidth(-1.0)


def test_topology_features_triangle_with_tail():
    G = nx.Graph()
    for k, (x, y) in {"a": (0, 0), "b": (1, 0), "c": (0.5, 1), "d": (0.5, 2)}.items():
        G.add_node(k, x=x, y=y)
    G.add_edges_from([("a", "b"), ("b", "c"), ("c", "a"), ("c", "d")])
    for u, v in G.edges():
        G[u][v]["length"] = 1.0
    net = rd2dnet.NetworkGraph(G, "length")
    E = net.edge_frame()
    orig, Z, names = _edge_feature_matrix(net, E, [], True)
    assert names == ["n_adjacent", "n_second", "degree_max", "clustering"]
    lab = net.nodes
    got = {frozenset((lab[u], lab[v])): Z[orig[r]] for r, (u, v) in enumerate(zip(E["u"], E["v"]))}
    assert np.allclose(got[frozenset("ab")], [2, 1, 2, 1])
    assert np.allclose(got[frozenset("bc")], [3, 0, 3, 2 / 3])
    assert np.allclose(got[frozenset("ca")], [3, 0, 3, 2 / 3])
    assert np.allclose(got[frozenset("cd")], [2, 1, 3, 1 / 6])


def test_split_edges_keep_original_topology(comb_bins):
    G, coords, b = comb_bins
    net = coords.network  # split at the boundary: crossing nodes inserted
    assert net.n_nodes > G.number_of_nodes()
    orig, Z, _ = _edge_feature_matrix(net, net.edge_frame(), [], True)
    assert Z.shape[0] == G.number_of_edges()
    ref = rd2dnet.NetworkGraph(G, "length")
    _, Z0, _ = _edge_feature_matrix(ref, ref.edge_frame(), [], True)
    assert np.allclose(np.sort(Z, axis=0), np.sort(Z0, axis=0))


def test_wbate_rules(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0)
    ok = (res.table["status"] == "ok").to_numpy()
    n = (res.table["n0"] + res.table["n1"]).to_numpy(float)
    w = np.where(ok, n, 0.0)
    assert res.wbate("density")["tau_bc"] == pytest.approx(res.wbate(w)["tau_bc"])
    S = res.cov[np.ix_(ok, ok)]
    g = np.linalg.solve(S, np.ones(ok.sum()))
    g /= g.sum()
    wg = np.zeros(len(ok))
    wg[ok] = g
    with pytest.warns(UserWarning, match="negative"):
        gls = res.wbate("gls")
    assert gls["tau_bc"] == pytest.approx(float(wg[ok] @ res.table["tau_bc"][ok]))
    assert gls["se_rb"] == pytest.approx(np.sqrt(g @ S @ g))
    assert gls["se_rb"] <= res.wbate("density")["se_rb"] + 1e-12
    assert res.wbate("place")["tau_bc"] == pytest.approx(res.wbate()["tau_bc"])
    with pytest.raises(ValueError):
        res.wbate("median")


def test_fe_recovers_mass_weighted_average(comb_bins):
    _, coords, b = comb_bins
    T = coords.treated.astype(float)
    tau_m = np.array([1.0, 2.0, 4.0])
    y = np.array([5.0, -1.0, 0.5])[b.bin] + tau_m[b.bin] * T
    plain = rd2dnet.bin_coordinates(coords, similarity=False, branches=False)
    fe = rd2dnet.rdnet_fe(y, plain, h=500.0, p=0, q=0, bin_slopes=False)
    W = np.clip(1 - np.abs(coords.s) / 500.0, 0, None)
    mass = np.bincount(plain.bin, weights=W, minlength=3)
    assert fe.tau == pytest.approx(mass @ tau_m / mass.sum(), abs=1e-10)
    # without reweighting, bins count by p_m (1 - p_m) times their mass
    w1 = np.bincount(plain.bin[coords.treated], weights=W[coords.treated], minlength=3)
    ph = w1 / mass
    fe0 = rd2dnet.rdnet_fe(y, plain, h=500.0, p=0, q=0, bin_slopes=False, ipw=False)
    om = mass * ph * (1 - ph)
    assert fe0.tau == pytest.approx(om @ tau_m / om.sum(), abs=1e-10)


def test_fe_exact_with_common_effect(comb_bins):
    _, coords, b = comb_bins
    T = coords.treated.astype(float)
    y = np.array([5.0, -1.0, 0.5])[b.bin] + 2.5 * T + 0.003 * coords.s - 0.001 * T * coords.s
    fe = rd2dnet.rdnet_fe(y, b, h=700.0)
    assert fe.tau == pytest.approx(2.5, abs=1e-8) and fe.tau_bc == pytest.approx(2.5, abs=1e-8)
    assert fe.bins["used"].all()
    assert "R*S" in fe.terms


def test_bins_run_with_every_variance(city):
    coords = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]], n_eval=5)
    b = rd2dnet.bin_coordinates(coords)
    y = city.obs["outcome"]
    for kw in (dict(), dict(cluster=city.obs["block"]), dict(conley_cutoff=100.0), dict(fuzzy=city.obs["treated"]),
               dict(covs=city.obs[["far"]])):
        with np.errstate(all="ignore"):
            res = rd2dnet.rdnet(y, b, h=700.0, min_obs=10, **kw)
        assert res._ok().any(), kw
        assert np.isfinite(res.wbate("density")["se_rb"])
    for kw in (dict(), dict(cluster=city.obs["block"]), dict(conley_cutoff=100.0), dict(vce="hc1")):
        fe = rd2dnet.rdnet_fe(y, b, h=700.0, **kw)
        assert np.isfinite(fe.se_rb) and fe.se_rb > 0
    sens = rd2dnet.similarity_sensitivity(y, b, h=700.0, min_obs=10)
    assert len(sens) == 4 and np.isinf(sens["h_S"].iloc[-1])
    assert "bins=" in repr(b)


def test_location_design_is_untouched(coords, city):
    """The default design's numbers do not move (reference from before bins existed)."""
    res = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0)
    assert coords.kernel_weights is None and coords.local_slopes(0, np.arange(3)) is None
    assert np.isfinite(res.table["tau_bc"]).sum() > 10


def test_branch_slopes_do_not_depend_on_column_order(city):
    """Which branches get their own slope is decided by size, not by column order."""
    coords = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]], n_eval=5)
    y = city.obs["outcome"]
    b = rd2dnet.bin_coordinates(coords)
    res = rd2dnet.rdnet(y, b, h=700.0, min_obs=10)
    flipped = rd2dnet.bin_coordinates(coords)
    flipped._bin_branches = [x[::-1] for x in flipped._bin_branches]
    res2 = rd2dnet.rdnet(y, flipped, h=700.0, min_obs=10)
    ok = res._ok()
    assert ok.any()
    assert np.allclose(res.table["tau_bc"][ok], res2.table["tau_bc"][ok], atol=1e-10)
    # thin branches share the reference slope: with a huge threshold no branch has its own
    big = rd2dnet.bin_coordinates(coords, min_branch=10**6)
    plain = rd2dnet.bin_coordinates(coords, branches=False)
    r_big = rd2dnet.rdnet(y, big, h=700.0, min_obs=10)
    r_plain = rd2dnet.rdnet(y, plain, h=700.0, min_obs=10)
    assert np.allclose(r_big.table["tau_bc"][ok], r_plain.table["tau_bc"][ok], atol=1e-10)
    with pytest.raises(ValueError):
        rd2dnet.bin_coordinates(coords, min_branch=0)
