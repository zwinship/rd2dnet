import numpy as np
import pytest

import rd2dnet
from rd2dnet.manipulation import _kmeans_1d_partitions, _lcg_labels


def test_reference_values(coords):
    # also checked by the R and Stata ports (same clusters and the same portable draws)
    r = rd2dnet.manipulation_test(coords, q=[50, 200], reps=2000)
    p = r.path
    assert p["statistic"].tolist() == pytest.approx([1.1666666666666667, 1.1666666666666667])
    assert p["pvalue"].tolist() == pytest.approx([0.2375, 0.689])
    assert p["n_treated"].tolist() == [23, 100]
    assert r.q == 50 and r.statistic == p["statistic"][0]


def test_portable_draws():
    assert _lcg_labels(0, 2, 10).tolist() == [[1, 1, 0, 1, 0, 1, 1, 0, 0, 0], [1, 0, 0, 1, 1, 0, 0, 1, 1, 1]]
    lab = _lcg_labels(7, 4000, 50)
    assert abs(lab.mean() - 0.5) < 0.01


def test_kmeans_dp_matches_exhaustive_search():
    import itertools

    rng = np.random.default_rng(1)
    for _ in range(100):
        pos = np.sort(rng.integers(0, 12, rng.integers(3, 9)).astype(float))
        vals = np.unique(pos)
        for k, f in enumerate(_kmeans_1d_partitions(pos, len(vals)), start=1):
            def sse(firsts):
                ends = list(firsts[1:]) + [len(pos)]
                return sum(((pos[a:e] - pos[a:e].mean()) ** 2).sum() for a, e in zip(firsts, ends))

            cands = [i for i in range(1, len(pos)) if pos[i] != pos[i - 1]]
            best = min(sse([0, *c]) for c in itertools.combinations(cands, k - 1))
            assert sse(list(f)) == pytest.approx(best, abs=1e-9)


def test_detects_localized_sorting(city):
    # control locations missing near one stretch of the boundary
    xy = city.obs[["x", "y"]].to_numpy()
    near = (np.abs(xy[:, 1] - 1600) < 150) & (xy[:, 1] > 1600) & (xy[:, 0] < 700)
    keep = ~(near & (np.random.default_rng(0).uniform(size=len(xy)) < 0.9))
    c = rd2dnet.network_coordinates(city.network("length"), city.boundary, xy[keep], n_eval=5)
    assert rd2dnet.manipulation_test(c, q=200, reps=2000).pvalue < 0.05


def test_input_checks(coords):
    with pytest.raises(ValueError, match="q must be"):
        rd2dnet.manipulation_test(coords, q=1)
    with pytest.raises(ValueError, match="level"):
        rd2dnet.manipulation_test(coords, level=0.95)
