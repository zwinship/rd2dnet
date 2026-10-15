"""Network-adapted McCrary-style density (sorting) test at the boundary.

Observations are binned by signed network distance ``s`` to the boundary and a
local linear fit of bin heights on each side estimates the density at ``s = 0``.
The test statistic is the log difference ``theta = log f+ - log f-``.

With ``normalize="network"`` bin counts are divided by the street length lying
in each distance band. On a network the amount of street at distance ``s``
varies with ``s`` and differs by side (a highway on one side, a dense grid on
the other), so raw counts can jump at the boundary even without sorting. The
normalized test compares observations per unit of street length instead.

Standard errors treat bin counts as Poisson.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .coordinates import NetworkCoordinates
from .kernels import get_kernel


@dataclass
class DensityResult:
    theta: float
    se: float
    z: float
    pvalue: float
    f_left: float
    f_right: float
    h: float
    binsize: float
    n_left: int
    n_right: int
    normalize: str
    bins: pd.DataFrame

    def __str__(self) -> str:
        return (
            f"rdnet density test ({self.normalize}): theta={self.theta:.4f} se={self.se:.4f} "
            f"z={self.z:.3f} p={self.pvalue:.4f}  [n_left={self.n_left}, n_right={self.n_right}, "
            f"h={self.h:.1f}, bin={self.binsize:.2f}]"
        )


def _local_linear_at_zero(mid, height, var, h, kfun):
    w = kfun(np.abs(mid) / h)
    keep = w > 0
    mid, height, var, w = mid[keep], height[keep], var[keep], w[keep]
    if len(mid) < 3:
        raise ValueError("too few non-empty bins within the bandwidth; increase h or the bin count")
    X = np.column_stack([np.ones_like(mid), mid])
    G = (X.T * w) @ X
    ell = np.linalg.solve(G, X.T * w)[0]
    return float(ell @ height), float(np.sum(ell**2 * var))


def density(
    coords: NetworkCoordinates,
    *,
    eval_index: int | None = None,
    window: float | None = None,
    h: float | None = None,
    binsize: float | None = None,
    normalize: str = "none",
    kernel="triangular",
) -> DensityResult:
    """Test for a discontinuity in the density of observations at the boundary.

    Parameters
    ----------
    eval_index, window
        Restrict to observations (and streets) whose boundary crossing lies
        within ``window`` (arc length) of evaluation point ``eval_index``.
        By default the whole boundary is pooled.
    h, binsize
        Bandwidth and bin width in metric units. Defaults follow McCrary (2008):
        ``binsize = 2 sd(s) n^(-1/2)`` and ``h = 2 sd(s) n^(-1/5)`` (at least 5 bins).
    normalize
        ``"none"`` (observations per unit of distance) or ``"network"``
        (observations per unit of street length).
    """
    if normalize not in ("none", "network"):
        raise ValueError("normalize must be 'none' or 'network'")
    kfun = get_kernel(kernel)
    s = np.abs(coords.s)
    T = coords.treated
    sel = np.isfinite(s)
    center = None
    if eval_index is not None:
        if window is None:
            raise ValueError("give a window (arc length) together with eval_index")
        center = coords.eval_node_arc[eval_index]
        sel &= np.abs(coords.crossing_arc - center) <= window
    s, T = s[sel], T[sel]
    n = len(s)
    if n < 20:
        raise ValueError("too few observations for a density test")
    signed = np.where(T, s, -s)
    sd = np.std(signed)
    binsize = 2 * sd * n ** -0.5 if binsize is None else float(binsize)
    h = max(2 * sd * n ** -0.2, 5 * binsize) if h is None else float(h)
    nb = int(np.ceil(h / binsize))
    edges = np.arange(nb + 1) * binsize

    exposure = {}
    if normalize == "network":
        net = coords.network
        E = net.edge_frame()
        nd = np.where(coords.node_on_boundary, 0.0, coords.node_dist)
        side = coords.node_treated
        u, v = E["u"].to_numpy(), E["v"].to_numpy()
        both_on = coords.node_on_boundary[u] & coords.node_on_boundary[v]
        e_side = np.where(coords.node_on_boundary[u], side[v], side[u])
        mid = (nd[u] + nd[v]) / 2
        keep = ~both_on & np.isfinite(mid)
        if center is not None:
            src = coords.node_source
            arc_u = coords.boundary.project(net.coords[np.maximum(src[u], 0)])
            arc_v = coords.boundary.project(net.coords[np.maximum(src[v], 0)])
            keep &= np.abs((arc_u + arc_v) / 2 - center) <= window
        for flag in (True, False):
            m = keep & (e_side == flag)
            exposure[flag] = np.histogram(mid[m], bins=edges, weights=E["length"].to_numpy()[m])[0]

    mids = (edges[:-1] + edges[1:]) / 2
    est = {}
    frames = []
    for flag in (True, False):
        counts = np.histogram(s[T == flag], bins=edges)[0].astype(float)
        if normalize == "network":
            ex = exposure[flag]
            ok = ex > 0
            height = np.where(ok, counts / np.where(ok, ex, 1), np.nan)
            var = np.where(ok, counts / np.where(ok, ex, 1) ** 2, np.nan)
        else:
            ok = np.ones_like(counts, bool)
            height = counts / (n * binsize)
            var = counts / (n * binsize) ** 2
        est[flag] = _local_linear_at_zero(mids[ok], height[ok], var[ok], h, kfun)
        frames.append(pd.DataFrame({"side": "treated" if flag else "control",
                                    "s": mids if flag else -mids, "count": counts, "height": height}))

    (f1, v1), (f0, v0) = est[True], est[False]
    if f1 <= 0 or f0 <= 0:
        raise ValueError("estimated density at the boundary is not positive; adjust h or binsize")
    theta = np.log(f1) - np.log(f0)
    se = float(np.sqrt(v1 / f1**2 + v0 / f0**2))
    zstat = theta / se
    return DensityResult(
        float(theta), se, float(zstat), float(2 * stats.norm.sf(abs(zstat))), f0, f1, h, binsize,
        int((~T).sum()), int(T.sum()), normalize,
        pd.concat(frames, ignore_index=True).sort_values("s", ignore_index=True),
    )
