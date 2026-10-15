import numpy as np
import pandas as pd
import pytest

import rd2dnet
from rd2dnet._wls import poly_basis, poly_exponents, wls
from rd2dnet.datasets import load_mock_midtown, make_mock_midtown
from rd2dnet.kernels import get_kernel


def test_poly_basis_terms():
    assert poly_exponents(1) == [(0, 0), (1, 0), (0, 1)]
    assert len(poly_exponents(2)) == 6
    X = poly_basis(np.array([2.0]), np.array([3.0]), 2)
    np.testing.assert_allclose(X, [[1, 2, 3, 4, 6, 9]])


def test_wls_recovers_linear_function():
    rng = np.random.default_rng(0)
    s, o = rng.normal(size=200), rng.normal(size=200)
    y = 1 + 2 * s - 3 * o
    f = wls(poly_basis(s, o, 1), y, rng.uniform(0.1, 1, 200))
    np.testing.assert_allclose(f.beta, [1, 2, -3], atol=1e-10)


@pytest.mark.parametrize("name", ["triangular", "epanechnikov", "uniform", "gaussian"])
def test_kernels(name):
    k = get_kernel(name)
    assert k(0.0) > 0
    assert k(np.inf) == 0
    if name != "gaussian":
        assert k(1.5) == 0


def test_bundled_data_matches_generator():
    a, b = load_mock_midtown(), make_mock_midtown()
    pd.testing.assert_frame_equal(a.obs, b.obs, check_exact=False, atol=1e-6)
    assert a.graph.number_of_edges() == b.graph.number_of_edges()


def test_rdnet_runs_and_masks(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, seed=0)
    t = res.table
    assert len(t) == coords.n_eval
    assert (t.loc[coords.eval_masked, "status"] == "masked").all()
    ok = t["status"] == "ok"
    assert ok.sum() >= 8
    assert (t.loc[ok, "se"] > 0).all()
    assert (t.loc[ok, "cb_lo"] <= t.loc[ok, "ci_lo"] + 1e-12).all()
    assert "rd2dnet" in res.summary()


def test_rdnet_recovers_large_effect():
    """With little noise and a constant effect, estimates sit near the truth."""
    city = make_mock_midtown(n_obs=6000, seed=3)
    obs = city.obs
    y = 0.0008 * (obs["y"] - 1600) + 3.0 * obs["treated"] + np.random.default_rng(1).normal(0, 0.1, len(obs))
    coords = rd2dnet.network_coordinates(city.network(), city.boundary, obs[["x", "y"]], n_eval=15)
    res = rd2dnet.rdnet(y.to_numpy(), coords, h=500.0, cbands=False)
    ok = res.table["status"] == "ok"
    assert ok.sum() >= 5
    np.testing.assert_allclose(res.table.loc[ok, "tau"], 3.0, atol=0.15)


def test_pooled_option_runs(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, pooled=True, cbands=False)
    assert (res.table["status"] == "ok").sum() >= 8


@pytest.mark.parametrize("vce", ["hc0", "hc1", "hc2", "hc3"])
def test_vce_variants(coords, city, vce):
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, vce=vce, cbands=False)
    ok = res.table["status"] == "ok"
    assert np.isfinite(res.table.loc[ok, "se"]).all()


def test_hc_ordering(coords, city):
    y = city.obs["outcome"].to_numpy()
    se = {v: rd2dnet.rdnet(y, coords, h=600.0, vce=v, cbands=False).table["se"] for v in ("hc0", "hc1", "hc3")}
    ok = se["hc0"].notna()
    assert (se["hc1"][ok] >= se["hc0"][ok]).all()
    assert (se["hc3"][ok] >= se["hc0"][ok]).all()


def test_cluster_and_conley(coords, city):
    y = city.obs["outcome"].to_numpy()
    rc = rd2dnet.rdnet(y, coords, h=600.0, cluster=city.obs["block"].to_numpy(), seed=0)
    assert rc.vce == "cluster" and (rc.table["status"] == "ok").sum() >= 8
    rn = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, cbands=False)
    re = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, conley_distance="euclidean", cbands=False)
    assert rn.vce == "conley"
    ok = rn.table["status"] == "ok"
    assert np.isfinite(rn.table.loc[ok, "se"]).all() and np.isfinite(re.table.loc[ok, "se"]).all()


def test_intensity_weights_change_estimates(coords, city):
    y = city.obs["outcome"].to_numpy()
    a = rd2dnet.rdnet(y, coords, h=600.0, cbands=False).table["tau"]
    b = rd2dnet.rdnet(y, coords, h=600.0, weights=city.obs["far"].to_numpy(), cbands=False).table["tau"]
    ok = a.notna() & b.notna()
    assert not np.allclose(a[ok], b[ok])


def test_bandwidth_selection(coords, city):
    y = city.obs["outcome"].to_numpy()
    bw = rd2dnet.rdbwnet(y, coords)
    ok = bw["status"] == "ok"
    assert ok.sum() >= 8
    assert (bw.loc[ok, "h0"] == bw.loc[ok, "h1"]).all()
    assert (bw.loc[ok, "h0"] >= np.maximum(bw.loc[ok, "hmin0"], bw.loc[ok, "hmin1"]) - 1e-9).all()
    assert np.isfinite(bw.loc[ok, ["V0", "V1", "B0", "B1", "Reg10", "Reg20", "bn0", "bn1", "dn"]]).all().all()


@pytest.mark.parametrize("sel", ["msetwo", "imserd", "imsetwo"])
def test_bandwidth_variants(coords, city, sel):
    bw = rd2dnet.rdbwnet(city.obs["outcome"].to_numpy(), coords, bwselect=sel)
    ok = bw["status"] == "ok"
    assert (bw.loc[ok, ["h0", "h1"]] > 0).all().all()
    if sel == "msetwo":
        assert not np.allclose(bw.loc[ok, "h0"], bw.loc[ok, "h1"])
    if sel.startswith("i"):
        h = bw.loc[ok, ["h0", "h1"]]
        inner = h[(h["h0"] > bw.loc[ok, "hmin0"]) & (h["h1"] > bw.loc[ok, "hmin1"])]
        assert inner["h0"].nunique() <= 2


def test_side_specific_bandwidths_in_estimation(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, bwselect="msetwo", cbands=False)
    ok = res.table["status"] == "ok"
    assert ok.sum() >= 6
    assert res.table.loc[ok, "h"].isna().any()  # h0 != h1 somewhere
    res2 = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=(400.0, 700.0), cbands=False)
    assert (res2.table.loc[res2.table["status"] == "ok", "h1"] == 700.0).all()


def test_bandwidth_with_cluster_and_fuzzy(coords, city):
    y = city.obs["outcome"].to_numpy()
    rng = np.random.default_rng(0)
    D = np.where(coords.treated, rng.uniform(size=len(y)) < 0.8, rng.uniform(size=len(y)) < 0.2).astype(float)
    bw = rd2dnet.rdbwnet(y, coords, cluster=city.obs["block"].to_numpy(), fuzzy=D)
    assert (bw.loc[bw["status"] == "ok", "h0"] > 0).all()


def test_fuzzy_design_recovers_late():
    city = make_mock_midtown(n_obs=6000, seed=11)
    obs = city.obs
    rng = np.random.default_rng(2)
    T = obs["treated"].to_numpy()
    D = np.where(T, rng.uniform(size=len(obs)) < 0.75, rng.uniform(size=len(obs)) < 0.25).astype(float)
    y = 0.0008 * (obs["y"] - 1600) + 2.0 * D + rng.normal(0, 0.1, len(obs))
    coords = rd2dnet.network_coordinates(city.network(), city.boundary, obs[["x", "y"]], n_eval=9)
    res = rd2dnet.rdnet(y.to_numpy(), coords, fuzzy=D, h=600.0, seed=0)
    t = res.table[res.table["status"] == "ok"]
    assert res.fuzzy and len(t) >= 4
    assert (np.abs(t["fs"] - 0.5) <= 4 * t["fs_se"]).all()
    assert (np.abs(t["tau"] - 2.0) <= 4 * t["se"]).all()
    assert ((t["ci_lo"] <= 2.0) & (2.0 <= t["ci_hi"])).mean() >= 0.75
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y.to_numpy(), coords, fuzzy=D[:-1])


def test_conley_uniform_band_uses_hac(coords, city):
    y = city.obs["outcome"].to_numpy()
    hc = rd2dnet.rdnet(y, coords, h=600.0, vce="hc0", seed=0)
    cl = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=1e-9, conley_resid_scale=1.0, seed=0)
    # A vanishing cutoff makes the HAC kernel diagonal: identical to HC0.
    np.testing.assert_allclose(cl.table["se_rb"], hc.table["se_rb"], rtol=1e-8)
    np.testing.assert_allclose(cl.cb_cov, hc.cb_cov, rtol=1e-8)
    wide = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=300.0, seed=0)
    assert wide.cb_crit is not None and np.isfinite(wide.table["cb_lo"]).any()


def test_hac_quadratic_blocking_invariant(coords):
    idx = np.arange(300)
    Psi = np.random.default_rng(0).normal(size=(300, 3))
    a = coords.hac_quadratic(idx, Psi, 200.0)
    b = coords.hac_quadratic(idx, Psi, 200.0, max_cells=7 * 300)
    np.testing.assert_allclose(a, b)
    K = np.clip(1 - coords.pairwise(idx, cutoff=200.0) / 200.0, 0, None)
    np.testing.assert_allclose(a, Psi.T @ K @ Psi, atol=1e-8)


def test_input_validation(coords, city):
    y = city.obs["outcome"].to_numpy()
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y[:-1], coords)
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y, coords, vce="cluster")
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y, coords, weights=-np.ones_like(y))


def test_travel_time_metric_end_to_end(city):
    net = city.network(rd2dnet.travel_time())
    c = rd2dnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]], n_eval=11)
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), c, cbands=False)
    assert (res.table["status"] == "ok").sum() >= 3
    # bandwidths are now in seconds, much smaller than metres
    assert res.table["h"].dropna().max() < 300


def test_plot(coords, city):
    import matplotlib

    matplotlib.use("Agg")
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, seed=0)
    ax = res.plot(true_tau=city.true_tau(res.table["node_arc"]))
    assert ax.get_xlabel()


@pytest.fixture(scope="module")
def big_coords():
    city = make_mock_midtown(n_obs=8000, seed=21)
    c = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]], eval_points=[600, 900, 2300, 2600])
    return city, c


def test_covariate_adjustment_matches_rd2d_plugin_and_shrinks_se(big_coords):
    from rd2dnet.kernels import triangular

    city, c = big_coords
    rng = np.random.default_rng(0)
    X = rng.normal(size=c.n_obs)
    y = 0.0005 * c.s + 2.0 * c.treated + 3.0 * X + rng.normal(0, 1, c.n_obs)
    plain = rd2dnet.rdnet(y, c, h=500.0, cbands=False).table
    adj = rd2dnet.rdnet(y, c, h=500.0, covs=X, cbands=False).table
    assert (adj["se"] < 0.6 * plain["se"]).all()
    # rd2d: common gamma from partialled-out side fits, then side-specific fits on y - X gamma
    d, o = c.at(0)
    w = triangular(d / 500.0)
    use = (w > 0) & np.isfinite(d)
    sides = [np.flatnonzero(use & c.treated), np.flatnonzero(use & ~c.treated)]
    A = B = 0.0
    for i in sides:
        R = poly_basis(c.s[i] / 500, o[i] / 500, 1)
        G = np.linalg.inv((R.T * w[i]) @ R)
        uz = (R.T * w[i]) @ X[i]
        A += (X[i] * w[i]) @ X[i] - uz @ G @ uz
        B += (X[i] * w[i]) @ y[i] - uz @ G @ ((R.T * w[i]) @ y[i])
    ya = y - (B / A) * X
    t1, t0 = (wls(poly_basis(c.s[i] / 500, o[i] / 500, 1), ya[i], w[i]).beta[0] for i in sides)
    assert adj["tau"][0] == pytest.approx(t1 - t0, rel=1e-8)


def test_constant_covariate_is_dropped(coords, city):
    y = city.obs["outcome"].to_numpy()
    a = rd2dnet.rdnet(y, coords, h=600.0, cbands=False).table
    b = rd2dnet.rdnet(y, coords, h=600.0, covs=np.ones(len(y)), cbands=False).table
    np.testing.assert_allclose(a["tau"], b["tau"])


def test_derivative_targets_exact_without_noise(big_coords):
    _, c = big_coords
    s, T = c.s, c.treated.astype(float)
    y = 0.001 * s + 0.003 * s * T + 1e-9 * np.random.default_rng(0).normal(size=c.n_obs)
    k = rd2dnet.rdnet(y, c, h=500.0, p=2, deriv=(1, 0), cbands=False).table
    np.testing.assert_allclose(k["tau"], 0.003, rtol=1e-4)
    np.testing.assert_allclose(k["tau_bc"], 0.003, rtol=1e-4)
    tv = rd2dnet.rdnet(y, c, h=500.0, p=2, tangvec=(1.0, 0.0), cbands=False).table
    np.testing.assert_allclose(tv["tau"], k["tau"])
    # a jump that grows along the boundary: d tau / d o = 0.002. Parcels that reach
    # the boundary at the highway node (via a ramp, not along 60th St) are left out.
    arc = c.crossing_arc
    y2 = np.where(arc >= 399.0, (0.5 + 0.002 * arc) * T, np.nan)
    g = rd2dnet.rdnet(y2 + 1e-9 * np.random.default_rng(1).normal(size=c.n_obs), c, h=600.0, p=1, deriv=(0, 1), cbands=False).table
    np.testing.assert_allclose(g["tau"], 0.002, rtol=1e-3)


def test_derivative_input_validation(coords, city):
    y = city.obs["outcome"].to_numpy()
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y, coords, p=1, deriv=(2, 0))
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y, coords, deriv=(1, 0), tangvec=(1, 0))
    with pytest.raises(ValueError):
        rd2dnet.rdnet(y, coords, pooled=True, deriv=(1, 0))


@pytest.mark.parametrize("sel,base", [("cerrd", "mserd"), ("certwo", "msetwo")])
def test_cer_bandwidths_shrink_mse(coords, city, sel, base):
    from rd2dnet.bandwidth import cer_factor

    y = city.obs["outcome"].to_numpy()
    m = rd2dnet.rdbwnet(y, coords, bwselect=base, bwcheck=1)
    c = rd2dnet.rdbwnet(y, coords, bwselect=sel, bwcheck=1)
    ok = (m["status"] == "ok") & (c["h0"] > c["hmin0"] + 1e-9)
    assert ok.any()
    ratio = (c.loc[ok, "h0"] / m.loc[ok, "h0"]).to_numpy()
    assert np.all(ratio < 1)
    if sel == "cerrd":
        np.testing.assert_allclose(ratio, cer_factor(len(np.unique(coords.xy, axis=0)), 1), rtol=1e-8)


def test_bandwidth_with_derivative_and_covariates(big_coords):
    city, c = big_coords
    y = city.obs["outcome"].to_numpy()
    X = np.random.default_rng(0).normal(size=c.n_obs)
    for kw in ({"p": 2, "deriv": (1, 0)}, {"p": 2, "tangvec": (0.0, 1.0)}, {"covs": X}):
        bw = rd2dnet.rdbwnet(y, c, **kw)
        assert (bw.loc[bw["status"] == "ok", "h0"] > 0).all()
    res = rd2dnet.rdnet(y, c, covs=X, bwselect="certwo", conley_cutoff=200.0, seed=0)
    assert (res.table["status"] == "ok").sum() >= 3


def test_conley_default_uses_wider_auxiliary_fit(coords, city):
    y = city.obs["outcome"].to_numpy()
    a = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, cbands=False).table
    b = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, conley_resid_scale=3.0, cbands=False).table
    c = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, conley_resid_order=3, cbands=False).table
    np.testing.assert_allclose(a["se_rb"], b["se_rb"])
    assert not np.allclose(a["se_rb"].dropna(), c["se_rb"].dropna())
