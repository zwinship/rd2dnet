import numpy as np
import pandas as pd
import pytest

import rd2dnet
from shapely.geometry import box
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
    np.testing.assert_allclose(cl.cov, hc.cov, rtol=1e-8)
    wide = rd2dnet.rdnet(y, coords, h=600.0, conley_cutoff=300.0, seed=0)
    assert wide.crit is not None and np.isfinite(wide.table["cb_lo"]).any()


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


@pytest.mark.parametrize("vce", ["hc1", "cluster", "conley"])
def test_joint_covariance_is_full_and_matches_pointwise(coords, city, vce):
    y = city.obs["outcome"].to_numpy()
    kw = {"cluster": city.obs["block"].to_numpy()} if vce == "cluster" else {}
    kw.update(conley_cutoff=150.0) if vce == "conley" else None
    res = rd2dnet.rdnet(y, coords, h=600.0, vce=vce, cbands=False, seed=0, **kw)
    ok = res.table["status"].eq("ok").to_numpy()
    assert res.cov.shape == (coords.n_eval, coords.n_eval) and not ok.all()
    assert np.isnan(res.cov[~ok]).all() and np.isnan(res.cov[:, ~ok]).all()
    np.testing.assert_allclose(np.sqrt(np.diag(res.cov)[ok]), res.table.loc[ok, "se_rb"], rtol=1e-10)
    assert "cb_lo" not in res.table and res.crit > 1.96  # stored without bands, for LBATE


def test_wbate_one_hot_and_uniform(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, seed=0)
    t, ok = res.table, res.table["status"].eq("ok").to_numpy()
    j = int(np.flatnonzero(ok)[2])
    w = np.zeros(coords.n_eval)
    w[j] = 5.0  # normalised to one
    one = res.wbate(w)
    for a, b in [("tau", "tau"), ("tau_bc", "tau_bc"), ("se_rb", "se_rb"), ("ci_lo", "ci_lo"),
                 ("ci_hi", "ci_hi"), ("pvalue", "pvalue")]:
        np.testing.assert_allclose(one[a], t.loc[j, b], rtol=1e-10)
    avg = res.wbate(ok.astype(float))
    np.testing.assert_allclose(avg["tau_bc"], t.loc[ok, "tau_bc"].mean())
    Sigma = res.cov[np.ix_(ok, ok)]
    np.testing.assert_allclose(avg["se_rb"], np.sqrt(Sigma.sum()) / ok.sum())
    # positively correlated neighbours: averaging helps less than independence would suggest
    assert avg["se_rb"] > np.sqrt((t.loc[ok, "se_rb"] ** 2).sum()) / ok.sum()
    with pytest.raises(ValueError, match="not estimated"):
        res.wbate(np.ones(coords.n_eval))
    for bad in (np.ones(3), np.zeros(coords.n_eval), np.r_[np.nan, ok[1:]]):
        with pytest.raises(ValueError):
            res.wbate(bad)


def test_lbate(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, seed=0)
    t = res.table[res.table["status"] == "ok"]
    lb = res.lbate()
    assert lb["tau_bc"] == t["tau_bc"].max() and lb["tau"] == t["tau"].max()
    assert lb["ci_lo"] <= lb["tau_bc"] <= lb["ci_hi"]
    np.testing.assert_allclose([lb["ci_lo"], lb["ci_hi"]], [t["cb_lo"].max(), t["cb_hi"].max()])
    # same interval without bands in the table (same seed, same critical value)
    nb = rd2dnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, seed=0, cbands=False).lbate()
    np.testing.assert_allclose(nb[["ci_lo", "ci_hi"]], lb[["ci_lo", "ci_hi"]])
    text = res.summary(wbate=res.table["status"].eq("ok").astype(float), lbate=True)
    assert "WBATE" in text and "LBATE" in text


def test_covariate_adjustment_with_side_bandwidths_scales_weights_by_h2(big_coords):
    # rd2d divides each side's kernel weights by h_side^2; with h0 != h1 this sets the
    # relative weight of the two sides in the common covariate coefficient
    from rd2dnet.kernels import triangular

    city, c = big_coords
    rng = np.random.default_rng(1)
    X = rng.normal(size=c.n_obs) + 0.002 * c.s
    y = 0.0005 * c.s + 2.0 * c.treated + 3.0 * X + rng.normal(0, 1, c.n_obs)
    h0, h1 = 400.0, 700.0
    adj = rd2dnet.rdnet(y, c, h=(h0, h1), covs=X, cbands=False).table
    d, o = c.at(0)
    sides = []
    for t, h in ((True, h1), (False, h0)):
        w = triangular(d / h) / h**2
        i = np.flatnonzero((w > 0) & np.isfinite(d) & (c.treated == t))
        sides.append((i, w, h))
    A = B = 0.0
    for i, w, h in sides:
        R = poly_basis(c.s[i] / h, o[i] / h, 1)
        G = np.linalg.inv((R.T * w[i]) @ R)
        uz = (R.T * w[i]) @ X[i]
        A += (X[i] * w[i]) @ X[i] - uz @ G @ uz
        B += (X[i] * w[i]) @ y[i] - uz @ G @ ((R.T * w[i]) @ y[i])
    ya = y - (B / A) * X
    t1, t0 = (wls(poly_basis(c.s[i] / h, o[i] / h, 1), ya[i], w[i]).beta[0] for i, w, h in sides)
    assert adj["tau"][0] == pytest.approx(t1 - t0, rel=1e-8)


def test_bandwidth_floor_counts_locations_not_repeated_measurements(city):
    # a panel: every location observed three times; the bwcheck-th nearest location is
    # the same as without repeats, so the floor (hmin) must not shrink
    net = city.network("length")
    xy = city.obs[["x", "y"]].to_numpy()
    c1 = rd2dnet.network_coordinates(net, city.boundary, xy, n_eval=9)
    c3 = rd2dnet.network_coordinates(net, city.boundary, np.repeat(xy, 3, axis=0), n_eval=9)
    y = city.obs["outcome"].to_numpy()
    b1 = rd2dnet.rdbwnet(y, c1)
    b3 = rd2dnet.rdbwnet(np.repeat(y, 3), c3)
    ok = (b1["status"] == "ok").to_numpy()
    np.testing.assert_allclose(b3["hmin0"][ok], b1["hmin0"][ok], rtol=1e-12)
    np.testing.assert_allclose(b3["hmin1"][ok], b1["hmin1"][ok], rtol=1e-12)


def test_one_point_band_uses_the_normal_quantile(coords, city):
    res = rd2dnet.rdnet(city.obs["outcome"], rd2dnet.network_coordinates(
        city.network("length"), city.boundary, city.obs[["x", "y"]], eval_points=[520.0]), h=600.0)
    assert res.crit == pytest.approx(1.959963984540054, rel=1e-12)


def test_min_obs_counts_locations_not_repeated_measurements(city):
    net = city.network("length")
    xy = city.obs[["x", "y"]].to_numpy()
    y = city.obs["outcome"].to_numpy()
    c1 = rd2dnet.network_coordinates(net, city.boundary, xy, n_eval=9)
    c3 = rd2dnet.network_coordinates(net, city.boundary, np.repeat(xy, 3, axis=0), n_eval=9)
    s1 = rd2dnet.rdnet(y, c1, h=180.0, cbands=False).table["status"]
    s3 = rd2dnet.rdnet(np.repeat(y, 3), c3, h=180.0, cbands=False).table["status"]
    assert (s1 == "insufficient").any()
    assert list(s3) == list(s1)


def test_plot_accepts_a_scalar_or_callable_truth(coords, city):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    res = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0, cbands=False)
    res.plot(true_tau=2.0)
    res.plot(true_tau=lambda arc: 1.5 + 0 * arc)


def test_few_clusters_within_the_bandwidth_warn(coords, city):
    blocks = (city.obs["x"] // 700).to_numpy()  # a handful of wide blocks
    with pytest.warns(UserWarning, match="few clusters"):
        rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0, cluster=blocks, cluster_vce="cr1", seed=0)


def test_selector_treats_badly_conditioned_local_designs_as_singular():
    # plug-in constants from a Gram matrix with condition number ~1e10 can flip sign with rounding
    from rd2dnet.bandwidth import _inv_gram

    A = np.diag([1.0, 1e-10])
    with pytest.raises(np.linalg.LinAlgError):
        _inv_gram(A)
    np.testing.assert_allclose(_inv_gram(np.diag([1.0, 1e-6])), np.diag([1.0, 1e6]))


def test_identifiable_basis_drops_aliased_and_leverage_one_columns():
    from rd2dnet.estimate import _identifiable

    rng = np.random.default_rng(0)
    s = rng.uniform(0, 1, 40)
    o = np.r_[np.zeros(20), np.ones(20)]  # two crossings: o identified, o^2 aliased with o
    X = np.column_stack([np.ones(40), s, o, s * s, s * o, o * o])
    assert _identifiable(X, np.ones(40)).tolist() == [True, True, True, True, True, False]
    o3 = np.r_[np.zeros(20), np.ones(19), 2.0]  # a third crossing with one observation
    X3 = np.column_stack([np.ones(40), s, o3, s * s, s * o3, o3 * o3])
    keep = _identifiable(X3, np.ones(40))
    assert not keep[5]  # o^2 would be determined by that single observation (leverage 1)


def test_singular_bias_correction_fit_does_not_change_the_main_estimate(city):
    # noise-free quadratic design; at the boundary's end only 4 crossings are reachable on
    # the treated side, so the order-3 fit cannot identify every o term
    c = rd2dnet.network_coordinates(city.network("length"), city.boundary, city.obs[["x", "y"]], n_eval=15)
    T, s, arc = c.treated, c.s, c.crossing_arc
    y = 1 + 2.5 * T + 0.004 * s + 3e-6 * s**2 - 2e-6 * T * s**2 - 0.0015 * arc + 1e-6 * arc**2
    t = rd2dnet.rdnet(y, c, h=500.0, p=2, q=3, cbands=False).table
    ok = t["status"] == "ok"
    np.testing.assert_allclose(t.loc[ok, "tau"], 2.5, atol=1e-8)
    np.testing.assert_allclose(t.loc[ok, "tau_bc"], 2.5, atol=1e-6)


def test_lbate_kinds_and_place_weighted_wbate(coords, city):
    r = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0, vce="hc1")
    t = r.table[r.table["status"] == "ok"]
    lo, hi = t["tau_bc"] - r.crit * t["se_rb"], t["tau_bc"] + r.crit * t["se_rb"]
    mn, ab = r.lbate("min"), r.lbate("absmax")
    assert mn["tau_bc"] == pytest.approx(t["tau_bc"].min()) and mn["ci_hi"] == pytest.approx(hi.min())
    assert ab["ci_lo"] == pytest.approx(np.maximum.reduce([0 * lo, lo, -hi]).max())
    assert ab["ci_hi"] == pytest.approx(np.maximum(lo.abs(), hi.abs()).max())
    with pytest.raises(ValueError, match="kind"):
        r.lbate("median")
    # with a wide bandwidth the two evaluation points on the isolated highway node are
    # estimated; they repeat one estimate and share one place's weight
    r2 = rd2dnet.rdnet(city.obs["outcome"], coords, h=2000.0, vce="hc1")
    ok = (r2.table["status"] == "ok").to_numpy()
    w = ok.astype(float)
    w[1] = 0.0  # the second copy of the highway place
    assert r2.wbate()["tau_bc"] == pytest.approx(r2.wbate(w)["tau_bc"])


def test_default_seed_makes_the_band_reproducible(coords, city):
    a = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0)
    b = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0)
    assert a.crit == b.crit


def test_points_on_the_line_are_treated_with_either_side_rule():
    line = [(0.0, 0.0), (100.0, 0.0)]
    on = np.array([[50.0, 0.0]])
    for side in ("left", "right"):
        assert rd2dnet.Boundary(line, treated_side=side).is_treated(on)[0]
    assert rd2dnet.Boundary(line, treated_region=box(-10, 0, 110, 50)).is_treated(on)[0]


def test_pooled_fit_with_side_bandwidths_is_exact_without_noise(coords):
    y = 2.0 * coords.treated + 0.004 * coords.s
    t = rd2dnet.rdnet(y, coords, h=(400.0, 700.0), pooled=True, cbands=False).table
    ok = t["status"] == "ok"
    np.testing.assert_allclose(t.loc[ok, "tau"], 2.0, atol=1e-9)


def test_covariate_units_do_not_change_the_fit(coords, city):
    z = city.obs["far"].to_numpy() + np.random.default_rng(0).normal(size=coords.n_obs)
    y = city.obs["outcome"]
    a = rd2dnet.rdnet(y, coords, h=600.0, covs=z, cbands=False).table
    b = rd2dnet.rdnet(y, coords, h=600.0, covs=z * 1e7, cbands=False).table
    assert (a["status"] == b["status"]).all() and (b["status"] == "ok").sum() > 5
    ok = a["status"] == "ok"
    np.testing.assert_allclose(b.loc[ok, "tau"], a.loc[ok, "tau"], rtol=1e-8)


def test_conley_cutoff_above_the_bandwidth_warns(coords, city):
    with pytest.warns(UserWarning, match="exceeds the bandwidth"):
        rd2dnet.rdnet(city.obs["outcome"], coords, h=400.0, conley_cutoff=1500.0, cbands=False)


def test_fuzzy_points_with_a_first_stage_near_zero_are_flagged(city, coords):
    rng = np.random.default_rng(0)
    d = rng.binomial(1, 0.5, coords.n_obs).astype(float)  # take-up unrelated to the boundary
    r = rd2dnet.rdnet(city.obs["outcome"], coords, h=600, fuzzy=d, cbands=False)
    t = r.table
    weak = t["status"] == "weak first stage"
    assert weak.any()
    assert (np.abs(t.loc[weak, "fs"]) < 2 * t.loc[weak, "fs_se"]).all()
    assert np.isfinite(t.loc[weak, "tau_bc"]).all()  # estimates are kept
    strong = t["status"] == "ok"
    assert (np.abs(t.loc[strong, "fs"]) >= 2 * t.loc[strong, "fs_se"]).all()


def test_rdbwnet_reproduces_the_bandwidths_rdnet_selects_with_cluster_and_conley(city, coords):
    y = city.obs["outcome"]
    cl = (city.obs["x"] // 500).to_numpy()
    r = rd2dnet.rdnet(y, coords, cluster=cl, cbands=False)
    bw = rd2dnet.rdbwnet(y, coords, cluster=cl, vce="cluster", scaleregul=3)
    np.testing.assert_allclose(r.table["h0"], bw["h0"], rtol=1e-12)
    r = rd2dnet.rdnet(y, coords, conley_cutoff=150.0, cbands=False)
    bw = rd2dnet.rdbwnet(y, coords, vce="conley", scaleregul=3)
    np.testing.assert_allclose(r.table["h0"], bw["h0"], rtol=1e-12)


def _cr2_dense(X, w, c, y, cl):
    """Bell-McCaffrey CR2 and its degrees of freedom from the dense n x n matrices."""
    n = len(y)
    Ginv = np.linalg.inv((X.T * w) @ X)
    B = Ginv @ (X.T * w)
    M = np.eye(n) - X @ B
    a = c @ B
    e = M @ y
    R = np.empty(n)
    for g in np.unique(cl):
        m = np.flatnonzero(cl == g)
        ev, U = np.linalg.eigh((M @ M.T)[np.ix_(m, m)])
        ev = np.where(ev > 1e-10 * ev.max(), 1 / np.sqrt(np.clip(ev, 1e-300, None)), 0.0)
        R[m] = (U * ev) @ (U.T @ a[m])
    labels = np.unique(cl)
    Rb = np.zeros((n, len(labels)))
    for j, g in enumerate(labels):
        Rb[cl == g, j] = R[cl == g]
    S = Rb.T @ e
    P = Rb.T @ M @ M.T @ Rb
    return float(np.sqrt(S @ S)), float(np.trace(P) ** 2 / np.sum(P * P))


def test_cr2_matches_the_dense_formulas():
    from rd2dnet._wls import wls

    rng = np.random.default_rng(4)
    n = 120
    X = np.c_[np.ones(n), rng.normal(size=(n, 2))]
    w = rng.uniform(0.2, 2.0, n)
    y = rng.normal(size=n)
    cl = rng.integers(0, 9, n)
    c = np.array([1.0, 0.3, 0.0])
    f = wls(X, y, w)
    R, df = f.cr2(c, cl)
    _, codes = np.unique(cl, return_inverse=True)
    se = np.sqrt(np.sum(np.bincount(codes, weights=R * f.resid) ** 2))
    se_ref, df_ref = _cr2_dense(X, w, c, y, cl)
    np.testing.assert_allclose([se, df], [se_ref, df_ref], rtol=1e-10)


def test_cluster_errors_default_to_cr2_with_t_intervals(city, coords):
    from scipy import stats

    blocks = (city.obs["x"] // 300).to_numpy() * 100 + (city.obs["y"] // 300).to_numpy()
    y = city.obs["outcome"]
    r2 = rd2dnet.rdnet(y, coords, h=600.0, cluster=blocks, seed=0)
    r1 = rd2dnet.rdnet(y, coords, h=600.0, cluster=blocks, cluster_vce="cr1", seed=0)
    t2, t1 = r2.table, r1.table
    ok = t2["status"] == "ok"
    assert "df" in t2 and "df" not in t1
    np.testing.assert_allclose(t2.loc[ok, "tau_bc"], t1.loc[ok, "tau_bc"], rtol=1e-12)
    q = stats.t.ppf(0.975, t2.loc[ok, "df"])
    np.testing.assert_allclose(t2.loc[ok, "ci_hi"] - t2.loc[ok, "tau_bc"], q * t2.loc[ok, "se_rb"], rtol=1e-10)
    np.testing.assert_allclose(np.sqrt(np.diag(r2.cov)[ok.to_numpy()]), t2.loc[ok, "se_rb"], rtol=1e-10)
    assert (t2.loc[ok, "se_rb"] > t1.loc[ok, "se_rb"]).mean() > 0.5  # CR1 is biased down
    assert np.isfinite(r2.wbate()["ci_lo"]) and np.isfinite(r2.lbate()["ci_hi"])
    with pytest.raises(ValueError, match="cluster_vce"):
        rd2dnet.rdnet(y, coords, h=600.0, cluster=blocks, cluster_vce="cr3")


def test_wbate_ignores_points_whose_first_stage_is_exactly_zero(city, coords):
    d = np.where(city.obs["x"].to_numpy() < 1200, 0.0, (city.obs["treated"].to_numpy() > 0).astype(float))
    r = rd2dnet.rdnet(city.obs["outcome"], coords, h=600.0, fuzzy=d, cbands=False)
    assert r.table["tau_bc"].isna().any()
    out = r.wbate()
    assert np.isfinite(out["tau_bc"]) and np.isfinite(out["se_rb"])
