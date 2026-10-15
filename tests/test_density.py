import numpy as np
import pytest

import rdnet
from rdnet.datasets import make_mock_midtown


def test_density_no_sorting(coords):
    r = rdnet.density(coords)
    assert np.isfinite(r.theta) and r.se > 0
    assert r.pvalue > 0.01
    assert r.n_left + r.n_right == coords.n_obs


def test_density_network_normalization(coords):
    r = rdnet.density(coords, normalize="network")
    assert r.normalize == "network" and np.isfinite(r.theta)


def test_density_detects_bunching():
    city = make_mock_midtown(n_obs=3000, seed=5)
    obs = city.obs
    rng = np.random.default_rng(0)
    extra = obs.sample(1500, random_state=1)[["x", "y"]].to_numpy().copy()
    extra[:, 1] = 1600 - rng.uniform(1, 60, len(extra))  # pile up just inside the treated zone
    pts = np.vstack([obs[["x", "y"]].to_numpy(), extra])
    c = rdnet.network_coordinates(city.network(), city.boundary, pts, n_eval=5)
    r = rdnet.density(c)
    assert r.theta > 0 and r.pvalue < 0.01


def test_density_local_window(coords):
    j = int(np.flatnonzero(~coords.eval_masked)[2])
    r = rdnet.density(coords, eval_index=j, window=600.0, h=400.0)
    assert r.n_left + r.n_right < coords.n_obs
    with pytest.raises(ValueError):
        rdnet.density(coords, eval_index=j)
