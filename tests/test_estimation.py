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
    assert (bw.loc[ok, "h"] >= bw.loc[ok, "h_min"]).all()
    assert (bw.loc[ok, "h"] > 0).all()
    common = rdnet.rdbwnet(y, coords, bwselect="imserd")
    assert common.loc[common["status"] == "ok", "h_mse"].nunique() == 1


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
