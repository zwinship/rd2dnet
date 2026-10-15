import numpy as np
import pandas as pd
import pytest

import rdnet
from rdnet._wls import poly_basis, poly_exponents, wls
from rdnet.datasets import load_mock_midtown, make_mock_midtown
from rdnet.kernels import get_kernel


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
    res = rdnet.rdnet(city.obs["outcome"].to_numpy(), coords, seed=0)
    t = res.table
    assert len(t) == coords.n_eval
    assert (t.loc[coords.eval_masked, "status"] == "masked").all()
    ok = t["status"] == "ok"
    assert ok.sum() >= 8
    assert (t.loc[ok, "se"] > 0).all()
    assert (t.loc[ok, "cb_lo"] <= t.loc[ok, "ci_lo"] + 1e-12).all()
    assert "rdnet" in res.summary()


def test_rdnet_recovers_large_effect():
    """With little noise and a constant effect, estimates sit near the truth."""
    city = make_mock_midtown(n_obs=6000, seed=3)
    obs = city.obs
    y = 0.0008 * (obs["y"] - 1600) + 3.0 * obs["treated"] + np.random.default_rng(1).normal(0, 0.1, len(obs))
    coords = rdnet.network_coordinates(city.network(), city.boundary, obs[["x", "y"]], n_eval=15)
    res = rdnet.rdnet(y.to_numpy(), coords, h=500.0, cbands=False)
    ok = res.table["status"] == "ok"
    assert ok.sum() >= 5
    np.testing.assert_allclose(res.table.loc[ok, "tau"], 3.0, atol=0.15)


def test_pooled_option_runs(coords, city):
    res = rdnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, pooled=True, cbands=False)
    assert (res.table["status"] == "ok").sum() >= 8


@pytest.mark.parametrize("vce", ["hc0", "hc1", "hc2", "hc3"])
def test_vce_variants(coords, city, vce):
    res = rdnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, vce=vce, cbands=False)
    ok = res.table["status"] == "ok"
    assert np.isfinite(res.table.loc[ok, "se"]).all()


def test_hc_ordering(coords, city):
    y = city.obs["outcome"].to_numpy()
    se = {v: rdnet.rdnet(y, coords, h=600.0, vce=v, cbands=False).table["se"] for v in ("hc0", "hc1", "hc3")}
    ok = se["hc0"].notna()
    assert (se["hc1"][ok] >= se["hc0"][ok]).all()
    assert (se["hc3"][ok] >= se["hc0"][ok]).all()


def test_cluster_and_conley(coords, city):
    y = city.obs["outcome"].to_numpy()
    rc = rdnet.rdnet(y, coords, h=600.0, cluster=city.obs["block"].to_numpy(), seed=0)
    assert rc.vce == "cluster" and (rc.table["status"] == "ok").sum() >= 8
    rn = rdnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, cbands=False)
    re = rdnet.rdnet(y, coords, h=600.0, conley_cutoff=150.0, conley_distance="euclidean", cbands=False)
    assert rn.vce == "conley"
    ok = rn.table["status"] == "ok"
    assert np.isfinite(rn.table.loc[ok, "se"]).all() and np.isfinite(re.table.loc[ok, "se"]).all()


def test_intensity_weights_change_estimates(coords, city):
    y = city.obs["outcome"].to_numpy()
    a = rdnet.rdnet(y, coords, h=600.0, cbands=False).table["tau"]
    b = rdnet.rdnet(y, coords, h=600.0, weights=city.obs["far"].to_numpy(), cbands=False).table["tau"]
    ok = a.notna() & b.notna()
    assert not np.allclose(a[ok], b[ok])


def test_bandwidth_selection(coords, city):
    y = city.obs["outcome"].to_numpy()
    bw = rdnet.rdbwnet(y, coords)
    ok = bw["status"] == "ok"
    assert ok.sum() >= 8
    assert (bw.loc[ok, "h0"] == bw.loc[ok, "h1"]).all()
    assert (bw.loc[ok, "h0"] >= np.maximum(bw.loc[ok, "hmin0"], bw.loc[ok, "hmin1"]) - 1e-9).all()
    assert np.isfinite(bw.loc[ok, ["V0", "V1", "B0", "B1", "Reg10", "Reg20", "bn0", "bn1", "dn"]]).all().all()


@pytest.mark.parametrize("sel", ["msetwo", "imserd", "imsetwo"])
def test_bandwidth_variants(coords, city, sel):
    bw = rdnet.rdbwnet(city.obs["outcome"].to_numpy(), coords, bwselect=sel)
    ok = bw["status"] == "ok"
    assert (bw.loc[ok, ["h0", "h1"]] > 0).all().all()
    if sel == "msetwo":
        assert not np.allclose(bw.loc[ok, "h0"], bw.loc[ok, "h1"])
    if sel.startswith("i"):
        h = bw.loc[ok, ["h0", "h1"]]
        inner = h[(h["h0"] > bw.loc[ok, "hmin0"]) & (h["h1"] > bw.loc[ok, "hmin1"])]
        assert inner["h0"].nunique() <= 2


def test_side_specific_bandwidths_in_estimation(coords, city):
    res = rdnet.rdnet(city.obs["outcome"].to_numpy(), coords, bwselect="msetwo", cbands=False)
    ok = res.table["status"] == "ok"
    assert ok.sum() >= 6
    assert res.table.loc[ok, "h"].isna().any()  # h0 != h1 somewhere
    res2 = rdnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=(400.0, 700.0), cbands=False)
    assert (res2.table.loc[res2.table["status"] == "ok", "h1"] == 700.0).all()


def test_bandwidth_with_cluster_and_fuzzy(coords, city):
    y = city.obs["outcome"].to_numpy()
    rng = np.random.default_rng(0)
    D = np.where(coords.treated, rng.uniform(size=len(y)) < 0.8, rng.uniform(size=len(y)) < 0.2).astype(float)
    bw = rdnet.rdbwnet(y, coords, cluster=city.obs["block"].to_numpy(), fuzzy=D)
    assert (bw.loc[bw["status"] == "ok", "h0"] > 0).all()


def test_fuzzy_design_recovers_late():
    city = make_mock_midtown(n_obs=6000, seed=11)
    obs = city.obs
    rng = np.random.default_rng(2)
    T = obs["treated"].to_numpy()
    D = np.where(T, rng.uniform(size=len(obs)) < 0.75, rng.uniform(size=len(obs)) < 0.25).astype(float)
    y = 0.0008 * (obs["y"] - 1600) + 2.0 * D + rng.normal(0, 0.1, len(obs))
    coords = rdnet.network_coordinates(city.network(), city.boundary, obs[["x", "y"]], n_eval=9)
    res = rdnet.rdnet(y.to_numpy(), coords, fuzzy=D, h=600.0, seed=0)
    t = res.table[res.table["status"] == "ok"]
    assert res.fuzzy and len(t) >= 4
    assert (np.abs(t["fs"] - 0.5) <= 4 * t["fs_se"]).all()
    assert (np.abs(t["tau"] - 2.0) <= 4 * t["se"]).all()
    assert ((t["ci_lo"] <= 2.0) & (2.0 <= t["ci_hi"])).mean() >= 0.75
    with pytest.raises(ValueError):
        rdnet.rdnet(y.to_numpy(), coords, fuzzy=D[:-1])


def test_conley_uniform_band_uses_hac(coords, city):
    y = city.obs["outcome"].to_numpy()
    hc = rdnet.rdnet(y, coords, h=600.0, vce="hc0", seed=0)
    cl = rdnet.rdnet(y, coords, h=600.0, conley_cutoff=1e-9, conley_resid_scale=1.0, seed=0)
    # A vanishing cutoff makes the HAC kernel diagonal: identical to HC0.
    np.testing.assert_allclose(cl.table["se_rb"], hc.table["se_rb"], rtol=1e-8)
    np.testing.assert_allclose(cl.cb_cov, hc.cb_cov, rtol=1e-8)
    wide = rdnet.rdnet(y, coords, h=600.0, conley_cutoff=300.0, seed=0)
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
        rdnet.rdnet(y[:-1], coords)
    with pytest.raises(ValueError):
        rdnet.rdnet(y, coords, vce="cluster")
    with pytest.raises(ValueError):
        rdnet.rdnet(y, coords, weights=-np.ones_like(y))


def test_travel_time_metric_end_to_end(city):
    net = city.network(rdnet.travel_time())
    c = rdnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]], n_eval=11)
    res = rdnet.rdnet(city.obs["outcome"].to_numpy(), c, cbands=False)
    assert (res.table["status"] == "ok").sum() >= 3
    # bandwidths are now in seconds, much smaller than metres
    assert res.table["h"].dropna().max() < 300


def test_plot(coords, city):
    import matplotlib

    matplotlib.use("Agg")
    res = rdnet.rdnet(city.obs["outcome"].to_numpy(), coords, h=600.0, seed=0)
    ax = res.plot(true_tau=city.true_tau(res.table["node_arc"]))
    assert ax.get_xlabel()
