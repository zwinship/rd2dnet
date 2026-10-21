"""Manipulation (sorting) test at the boundary, after Bugni, Crippa and Restrepo (2026).

Without sorting across the boundary, the locations closest to it are about
equally likely to lie on either side, anywhere along the boundary. The test
takes the ``q`` locations nearest the boundary (straight-line distance in the
plane: sorting is about where units locate, not about the network), projects
them onto the boundary, groups the projections into ``k`` contiguous clusters
for every ``k`` from 1 to ``k_max``, and measures the largest standardised
imbalance between the treated share in a cluster and 1/2:

    S = max_k max_a sqrt(n_ka) * |share treated in cluster a of partition k - 1/2|.

Its null distribution is simulated by relabelling the ``q`` locations as
independent Bernoulli(1/2) draws with the clusters held fixed. ``q`` is a fixed
number of locations, so the test is valid as the sample grows for any
boundary that is smooth except at finitely many corners; the paper recommends
looking at a range of ``q``.

Unlike a density test of the signed distance to the boundary, it needs no
assumption about how the distance to the boundary is distributed, which on a
curved boundary (a zone outline) or along a street network differs between
the two sides even without sorting.

Implementation choices (so that Python, R and Stata give identical answers):

- Projections are positions along the boundary (arc length), and the clusters
  are the exact optimal 1-D k-means partitions of those positions (dynamic
  programming over distinct positions; ties go to the earliest split). On a
  closed ring the positions are cut at their largest gap.
- The relabelling draws come from the Park-Miller generator
  ``x <- 16807 x mod (2^31 - 1)``, seeded by ``seed``; a draw is treated when
  ``x <= 2^30 - 1``.
- Locations count once (repeated measurements share a location). Locations on
  the boundary line, in masked stretches, or projecting onto an end point of
  an open boundary are left out.

Reference: Bugni, F. A., Crippa, F. and Restrepo, D. (2026). Manipulation
Testing in Boundary Discontinuity Designs. arXiv:2609.00350.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .bandwidth import _as_int
from .coordinates import NetworkCoordinates

_M = 2147483647  # 2^31 - 1
_A = 16807
_HALF = 1073741823  # draws 1..2^30-1 are "treated": exactly half of 1..2^31-2


@dataclass
class ManipulationTestResult:
    """Output of :func:`manipulation_test`.

    ``statistic`` is the largest standardised imbalance ``S``; ``pvalue`` the
    share of relabelled draws with a statistic at least as large; ``crit`` the
    simulated critical value at ``level``. ``distance`` is the distance to the
    boundary of the farthest selected location, ``n_treated`` the number of
    the ``q`` selected locations on the treated side, ``k_max`` the largest
    number of clusters. ``path`` holds one row per ``q`` when several were
    requested (the first row's values fill the scalar fields).
    """

    statistic: float
    pvalue: float
    crit: float
    q: int
    n_treated: int
    distance: float
    k_max: int
    reps: int
    level: float
    path: pd.DataFrame = field(repr=False, default=None)

    def __str__(self) -> str:
        return (
            f"rd2dnet manipulation test (Bugni, Crippa and Restrepo): q={self.q} locations within "
            f"{self.distance:.4g} of the boundary, {self.n_treated} treated; S={self.statistic:.3f}, "
            f"p={self.pvalue:.4f} (critical value {self.crit:.3f} at {self.level:g}%, {self.reps} draws)"
        )


def _lcg_labels(seed: int, reps: int, q: int) -> np.ndarray:
    """``reps x q`` Bernoulli(1/2) labels from the portable Park-Miller stream."""
    x0 = int(seed) % (_M - 1) + 1
    # state after n steps: x0 * A^n mod M; draw r, position j uses step n = r * q + j + 1
    pow_j = np.empty(q, dtype=np.int64)
    p = 1
    for j in range(q):
        p = p * _A % _M
        pow_j[j] = p  # A^(j+1)
    aq = int(pow_j[-1])  # A^q
    starts = np.empty(reps, dtype=np.int64)
    s = x0
    for r in range(reps):
        starts[r] = s  # x0 * A^(r q)
        s = s * aq % _M
    x = (starts[:, None] * pow_j[None, :]) % _M
    return (x <= _HALF).astype(np.int8)


def _kmeans_1d_partitions(pos: np.ndarray, k_max: int) -> list[np.ndarray]:
    """Optimal 1-D k-means partitions of sorted positions ``pos`` for k = 1..k_max.

    Clusters are contiguous runs of distinct positions; returns, for each k, the
    index (into ``pos``) of the first element of every cluster.
    """
    vals, start, counts = np.unique(pos, return_index=True, return_counts=True)
    m = len(vals)
    k_max = min(k_max, m)
    w = counts.astype(float)
    s1 = np.r_[0.0, np.cumsum(w * vals)]
    s2 = np.r_[0.0, np.cumsum(w * vals * vals)]
    sw = np.r_[0.0, np.cumsum(w)]

    def cost(i, j):  # distinct values i..j (0-based, inclusive); i may be an array
        n = sw[j + 1] - sw[i]
        a = s1[j + 1] - s1[i]
        return np.maximum((s2[j + 1] - s2[i]) - a * a / n, 0.0)

    D = np.full((k_max + 1, m), np.inf)
    arg = np.zeros((k_max + 1, m), dtype=np.int64)
    D[1] = cost(np.zeros(m, dtype=np.int64), np.arange(m))
    for k in range(2, k_max + 1):
        for j in range(k - 1, m):
            i = np.arange(k - 1, j + 1)  # the last cluster is values i..j
            c = D[k - 1, i - 1] + cost(i, j)
            best = c.min()
            bi = int(i[np.flatnonzero(c <= best + 1e-12 * max(1.0, abs(best)))[0]])
            D[k, j], arg[k, j] = c[bi - (k - 1)], bi
    out = []
    for k in range(1, k_max + 1):
        firsts = []
        j, kk = m - 1, k
        while kk >= 1:
            i = arg[kk, j] if kk > 1 else 0
            firsts.append(i)
            j, kk = i - 1, kk - 1
        out.append(start[np.array(firsts[::-1], dtype=np.int64)])
    return out


def _statistic(labels: np.ndarray, firsts: list[np.ndarray], q: int) -> np.ndarray:
    """``S`` for each row of ``labels`` (draws x q, in position order)."""
    P = np.concatenate([np.zeros((labels.shape[0], 1)), np.cumsum(labels, axis=1)], axis=1)
    S = np.zeros(labels.shape[0])
    for f in firsts:
        ends = np.r_[f[1:], q]
        n = (ends - f).astype(float)
        t = P[:, ends] - P[:, f]
        dev = np.sqrt(n) * np.abs(t / n - 0.5)
        S = np.maximum(S, dev.max(axis=1))
    return S


def manipulation_test(
    coords: NetworkCoordinates,
    *,
    q=None,
    k_max: int | None = None,
    reps: int = 10_000,
    level: float = 95.0,
    seed: int = 0,
) -> ManipulationTestResult:
    """Test for sorting across the boundary (Bugni, Crippa and Restrepo, 2026).

    Parameters
    ----------
    coords
        Network coordinates; the test uses the locations, the boundary and the
        treatment indicator (straight-line distances, not the network).
    q
        Number of locations nearest the boundary to use, or a list of values
        (the result's ``path`` then has one row per value). Defaults to 200
        (more power than 50 or 100 in our Monte Carlo, with correct size), or a
        tenth of the usable locations when there are fewer than 2,000. Report
        the path over several values, e.g. ``q=[50, 100, 200]``.
    k_max
        Largest number of clusters along the boundary (default: ``q``).
    reps
        Relabelling draws for the null distribution.
    level
        Confidence level (percent) for the reported critical value.
    seed
        Seed of the relabelling draws.

    Returns
    -------
    ManipulationTestResult

    Notes
    -----
    Designed for units located in the plane (parcels, households, firms). It is
    not valid for observations that lie on the network itself (stations along a
    river, points on streets): where a line crosses the boundary at an angle,
    position along the boundary sorts the nearest points by side, and in Monte
    Carlo runs the test rejected up to 94% of the time without sorting. It has
    power against sorting that offsets along the boundary (one side denser on
    one stretch, the other on another), which a density test of the distance
    to the boundary cannot detect; against sorting in one direction everywhere,
    ``density()`` is more powerful.
    """
    from .density import _warn_if_on_network

    _warn_if_on_network(coords, stacklevel=3)
    if not (1 < level < 100):
        raise ValueError(f"level is a percentage, e.g. 95 (got {level!r})")
    if int(reps) < 100:
        raise ValueError("reps must be at least 100")
    b = coords.boundary
    first = coords._first
    xy = coords.xy[first]
    T = np.asarray(coords.treated)[first].astype(bool)
    dist = b.distance(xy)
    arc = b.project(xy)
    L = b.length
    ok = np.isfinite(dist) & (dist > b.tol) & ~b.is_masked(arc)
    if not b.closed:
        ok &= (arc > 1e-9 * L) & (arc < L * (1 - 1e-9))
    idx = np.flatnonzero(ok)
    n_use = len(idx)
    qs = [q] if q is None or np.isscalar(q) else list(q)
    if qs == [None]:
        qs = [200 if n_use >= 2000 else max(n_use // 10, 1)]
    rows = []
    if k_max is not None:
        k_max = _as_int(k_max, "k_max", 1)
    for qq in qs:
        qq = _as_int(qq, "q")
        if qq < 2 or qq > n_use:
            what = "q" if q is not None else f"the default q ({qq}, a tenth of the usable locations)"
            raise ValueError(f"{what} must be between 2 and the number of usable locations ({n_use})")
        if qq < 20:
            warnings.warn(f"q = {qq} locations give a test with little power", stacklevel=2)
        if qq > 0.4 * n_use:
            warnings.warn(f"q = {qq} is a large share of the {n_use} usable locations; the test's approximation "
                          "assumes the selected locations are close to the boundary", stacklevel=2)
        order = np.lexsort((idx, dist[idx]))  # nearest first; ties by location order
        sel = idx[order[:qq]]
        pos = arc[sel]
        if b.closed:  # cut the ring at the largest gap between neighbouring positions
            srt = np.sort(np.unique(pos))
            gaps = np.diff(np.r_[srt, srt[0] + L])
            g = int(np.argmax(gaps))
            cut = srt[(g + 1) % len(srt)]
            pos = (pos - cut) % L
        po = np.lexsort((sel, pos))  # position order; ties by location order
        pos, lab = pos[po], T[sel][po].astype(np.int8)
        km = qq if k_max is None else int(min(k_max, qq))
        firsts = _kmeans_1d_partitions(pos, km)
        S = float(_statistic(lab[None, :].astype(float), firsts, qq)[0])
        null = _statistic(_lcg_labels(seed, int(reps), qq).astype(float), firsts, qq)
        pval = float(np.mean(null >= S - 1e-12 * max(1.0, S)))
        crit = float(np.sort(null)[int(np.ceil(level / 100 * len(null))) - 1])
        rows.append(dict(q=qq, statistic=S, pvalue=pval, crit=crit, n_treated=int(lab.sum()),
                         distance=float(dist[sel].max()), k_max=len(firsts)))
    path = pd.DataFrame(rows)
    r0 = rows[0]
    return ManipulationTestResult(statistic=r0["statistic"], pvalue=r0["pvalue"], crit=r0["crit"], q=r0["q"],
                                  n_treated=r0["n_treated"], distance=r0["distance"], k_max=r0["k_max"],
                                  reps=int(reps), level=float(level), path=path)
