"""Network-adapted density (sorting) test at the boundary.

The estimator and jackknife variance are ported from rddensity (Cattaneo,
Jansson and Ma; GPL-3.0; ``rddensity/funs.py``) and adapted to networks.

Uses the local polynomial density estimator of Cattaneo, Jansson and Ma
(2020; the ``rddensity`` package) on the signed network distance ``s``: the
empirical CDF is regressed on a polynomial in ``s`` separately on each side
within a bandwidth, and the density at the boundary is the slope at ``s = 0``.
Standard errors use rddensity's jackknife variance; inference uses the
bias-corrected order ``q = p + 1`` fit.

With ``normalize="network"`` each observation's distance is re-expressed as
cumulative street length: ``G(a)`` is the total length of street on that side
lying within network distance ``a`` of the boundary, and the test is run on
``sign * G(|s|)``. The density in that coordinate is observations per metre
of street, so jumps in street supply at the boundary (a highway on one side,
a dense grid on the other, the gap before the first parallel street) do not
masquerade as sorting. Bandwidths are then in street-length units. Street
geometry is treated as fixed.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy import stats

from .bandwidth import _check_weights
from .coordinates import NetworkCoordinates
from .kernels import get_kernel


def _street_closeness(coords: NetworkCoordinates, n_ref: int = 45) -> float:
    """Median snap distance of the locations relative to points spread evenly over their extent.

    About 1 when units are spread over the plane (parcels, households), near 0
    when they sit on or next to the streets (stations on a river, GPS fixes of
    vehicles). Scale-free: it compares with the same network's own spacing.
    """
    depth = np.asarray(coords.snap.depth, dtype=float)[coords._first]
    xy = coords.xy[coords._first]
    ok = np.isfinite(depth) & np.all(np.isfinite(xy), axis=1)
    if ok.sum() < 20:
        return np.nan
    lo, hi = xy[ok].min(axis=0), xy[ok].max(axis=0)
    g = (np.arange(n_ref) + 0.5) / n_ref  # a regular 45 x 45 lattice over the locations' extent
    gx, gy = np.meshgrid(lo[0] + g * (hi[0] - lo[0]), lo[1] + g * (hi[1] - lo[1]))
    ref_pts = np.column_stack([gx.ravel(), gy.ravel()])
    try:
        ref = np.asarray(coords.network.snap(ref_pts).depth, dtype=float)
    except ValueError:
        return np.nan
    med = float(np.median(ref))
    if not med > 0:
        return np.nan
    return float(np.median(depth[ok])) / med


def _warn_if_on_network(coords: NetworkCoordinates, stacklevel: int = 2) -> None:
    """Warn when most locations lie on or next to the network itself (no sorting test is valid then)."""
    r = _street_closeness(coords)
    if np.isfinite(r) and r < 0.6:
        where = "on the network itself" if r < 1e-6 else "on or next to the network"
        warnings.warn(f"most locations lie {where} (stations on a river, GPS points on streets): their median "
                      f"distance to the nearest street is {r:.2f} times that of points spread over the same area. "
                      "Sorting tests assume units spread over the plane; in Monte Carlo runs none of them held its "
                      "size for points on or near streets (GPS noise up to 15 m)", stacklevel=stacklevel)


@dataclass
class DensityResult:
    """Output of :func:`density`.

    The test statistic is rddensity's jump ``f+ - f-`` in the density at the
    boundary (per unit of ``s``, or per metre of street with
    ``normalize="network"``). ``diff`` is the order-``p`` estimate; ``diff_bc``,
    ``se``, ``z`` and ``pvalue`` come from the bias-corrected order-``q`` fit.
    ``theta`` is the descriptive log ratio ``log(f+/f-)`` of the order-``p``
    estimates (``nan`` if either is not positive).
    """

    diff: float
    diff_bc: float
    se: float
    z: float
    pvalue: float
    theta: float
    f_left: float
    f_right: float
    h_left: float
    h_right: float
    n_left: int
    n_right: int
    n_left_h: int
    n_right_h: int
    p: int
    q: int
    normalize: str

    def __str__(self) -> str:
        return (
            f"rd2dnet density test ({self.normalize}): jump={self.diff:.3g} (bc {self.diff_bc:.3g}) "
            f"se={self.se:.3g} z={self.z:.3f} p={self.pvalue:.4f} log-ratio={self.theta:.3f} "
            f"[n={self.n_left}/{self.n_right}, in window {self.n_left_h}/{self.n_right_h}, "
            f"h={self.h_left:.1f}/{self.h_right:.1f}]"
        )


def _cjm(x, om, hl, hr, p, kfun):
    """Local polynomial CDF fit; returns raw coefficients, their jackknife covariance, window counts.

    Columns are ordered ``[1_L, 1_R, u_L, u_R, u_L^2, u_R^2, ...]``; the
    density on each side is the coefficient on ``u`` (indices 2 and 3).
    """
    order = np.argsort(x, kind="stable")
    x, om = x[order], om[order]
    denom = om.sum() - om.mean()
    F = (np.cumsum(om) - om) / denom
    win = (x >= -hl) & (x <= hr)
    xw, Fw, omw = x[win], F[win], om[win]
    left = xw < 0
    u = np.where(left, xw / hl, xw / hr)
    K = np.where(left, kfun(np.abs(u)) / hl, kfun(np.abs(u)) / hr)
    keep = K > 0
    xw, Fw, omw, left, u, K = xw[keep], Fw[keep], omw[keep], left[keep], u[keep], K[keep]
    nl, nr = int(left.sum()), int((~left).sum())
    if min(nl, nr) <= p + 1:
        raise ValueError("too few observations within the bandwidth on one side")
    Xp = np.zeros((len(xw), 2 * p + 2))
    Hp = np.zeros(2 * p + 2)
    for k in range(p + 1):
        Xp[:, 2 * k] = np.where(left, u**k, 0.0)
        Xp[:, 2 * k + 1] = np.where(~left, u**k, 0.0)
        Hp[2 * k], Hp[2 * k + 1] = hl**k, hr**k
    XpW = Xp * K[:, None]
    Sinv = np.linalg.inv(XpW.T @ Xp)
    b = (Sinv @ (XpW.T @ Fw)) / Hp
    tail = np.cumsum(XpW[::-1], axis=0)[::-1]  # sum over k with x_k >= x_i
    L = tail * (omw / denom)[:, None]
    Hinv = np.diag(1.0 / Hp)
    V = Hinv @ Sinv @ (L.T @ L) @ Sinv @ Hinv
    return b, V, nl, nr


def _side_bandwidth(x, om, right, p, kfun, scaleregul, n_min):
    """MSE-optimal bandwidth for one side's density (plug-in, rd2d-style constants)."""
    side = x >= 0 if right else x < 0
    a = np.abs(x[side])
    if a.size < n_min:
        raise ValueError("too few observations on one side for bandwidth selection")
    sd = np.std(x)
    lo = float(np.partition(a, n_min - 1)[n_min - 1])
    hi = float(a.max())
    h0 = min(max(2.0 * sd * len(x) ** (-1.0 / (2 * p + 1)), lo), hi)
    bn = min(max(float(np.median(a)), lo), hi)
    sgn = 3 if right else 2  # column of the density coefficient
    # Variance constant at the pilot.
    _, V0, _, _ = _cjm(x, om, h0, h0, p, kfun)
    Cv = V0[sgn, sgn] * h0
    # Bias constant: degree p+1 coefficient from an order p+1 fit at bn, projected on the slope.
    b1, V1, _, _ = _cjm(x, om, bn, bn, p + 1, kfun)
    gi = 2 * (p + 1) + (sgn - 2)
    gamma, var_gamma = b1[gi], V1[gi, gi]
    # Projection of u^(p+1) onto the slope at h0, in scaled units.
    xs = x[(x >= -h0) & (x <= h0)]
    left = xs < 0
    u = xs / h0
    K = kfun(np.abs(u)) / h0
    on = ~left if sgn == 3 else left
    us, Ks = u[on], K[on]
    X = np.column_stack([us**k for k in range(p + 1)])
    proj = np.linalg.solve((X.T * Ks) @ X, (X.T * Ks) @ us ** (p + 1))[1]
    Cb2 = (proj * gamma) ** 2 + scaleregul * proj**2 * var_gamma
    h = (Cv / (2 * p * max(Cb2, 1e-300))) ** (1.0 / (2 * p + 1))
    return float(min(max(h, lo), hi))


def _arc_gap(boundary, arc, center):
    """Distance along the boundary between arc positions (the shorter way round a closed ring)."""
    gap = np.abs(np.asarray(arc, dtype=float) - center)
    if boundary.closed:
        gap = np.minimum(gap, boundary.length - gap)
    return gap


def _street_length_transform(coords: NetworkCoordinates, center, window):
    """Return ``G_left, G_right``: cumulative street length within distance ``a`` of the boundary."""
    net = coords.network
    E = net.edge_frame()
    u, v = E["u"].to_numpy(), E["v"].to_numpy()
    on = coords.node_on_boundary
    nd = np.where(on, 0.0, coords.node_dist)
    side = np.where(on[u], coords.node_treated[v], coords.node_treated[u])
    cost = E["cost"].to_numpy()
    length = E["length"].to_numpy()
    keep = ~(on[u] & on[v]) & np.isfinite(nd[u]) & np.isfinite(nd[v])
    src = coords.node_source
    arc = (coords.boundary.project(net.coords[np.maximum(src[u], 0)])
           + coords.boundary.project(net.coords[np.maximum(src[v], 0)])) / 2
    keep &= ~coords.boundary.is_masked(arc)
    if center is not None:
        keep &= _arc_gap(coords.boundary, arc, center) <= window
    exact = np.isfinite(cost)
    med = np.nanmedian(cost[keep]) if np.isfinite(cost[keep]).any() else np.nanmedian(np.abs(nd[u] - nd[v])[keep])
    out = []
    for flag in (False, True):
        ks = np.flatnonzero(keep & (side == flag))
        if ks.size == 0:
            raise ValueError("no streets on one side of the boundary in this window")
        a, b, ex = nd[u[ks]], nd[v[ks]], exact[ks]
        c = np.where(ex, cost[ks], np.abs(b - a))
        # Sample each street at m points (more on long ones) and spread its length over them.
        m = np.maximum(4, np.ceil(16 * c / max(med, 1e-12)).astype(np.int64))
        e = np.repeat(np.arange(ks.size), m)
        tt = (np.arange(m.sum()) - np.repeat(np.cumsum(m) - m, m) + 0.5) / m[e]
        xs = np.where(ex[e], np.minimum(a[e] + tt * c[e], b[e] + (1 - tt) * c[e]), a[e] + tt * (b[e] - a[e]))
        ws = (length[ks] / m)[e]
        order = np.argsort(xs, kind="stable")
        grid, cum = np.r_[0.0, xs[order]], np.r_[0.0, np.cumsum(ws[order])]
        out.append(lambda a, g=grid, c=cum: np.interp(a, g, c))
    return out


def density(
    coords: NetworkCoordinates,
    *,
    eval_index: int | None = None,
    window: float | None = None,
    h=None,
    p: int = 2,
    q: int | None = None,
    normalize: str = "euclidean",
    kernel="triangular",
    weights=None,
    bwselect: str = "common",
    scaleregul: float = 1.0,
    unique_locations: bool = False,
) -> DensityResult:
    """Test for a discontinuity in the density of observations at the boundary.

    Parameters
    ----------
    eval_index, window
        Restrict to observations (and streets) whose boundary crossing lies
        within ``window`` (arc length) of evaluation point ``eval_index``.
        By default the whole boundary is pooled. Masked stretches of the
        boundary are always excluded.
    h
        Bandwidth, scalar or ``(h_control, h_treated)``, in metric units (street-length
        units with ``normalize="network"``). Selected by an MSE plug-in rule per
        side when omitted.
    p, q
        Local polynomial order of the CDF fit (density is its slope) and of the
        bias-corrected fit used for inference (default ``p + 1``).
    normalize
        The running variable. ``"euclidean"`` (default): the signed
        straight-line distance from each location to the boundary line; sorting
        is about where units locate, and this is the test that held its size
        (0.03-0.06 at 5%) on straight, round, wavy and irregular boundaries in
        our Monte Carlo, with the most power against one-sided sorting.
        ``"none"``: the signed network distance ``s`` (as rddensity would on
        ``s``); its density jumps at curved or irregular boundaries without any
        sorting (rejection 0.08-0.33). ``"network"``: observations per metre of
        street; valid only when observations lie on the streets of a straight
        boundary, and it rejects 34-100% of the time at round or wavy zone
        outlines. No density test holds its size everywhere when the
        observations lie on or next to the network (points on streets): see
        :func:`manipulation_test` and the README.
    weights
        Optional non-negative observation weights (e.g. floor area).
    bwselect
        ``"common"`` (default): one bandwidth, the smaller of the two side-wise
        MSE-optimal choices, for both sides. ``"each"``: separate bandwidths.
        On street grids the density of ``s`` is jagged, and unequal bandwidths
        let unequal smoothing bias pass for a discontinuity.
    scaleregul
        Weight on the variance of the bias estimate in the bandwidth rule.
    unique_locations
        Count each location once (the first observation there). Use this with
        repeated measurements (panels, repeat sales, sensor readings): the
        sorting question is about where units are, and repeated rows at the
        same place otherwise look like a density spike.

    Returns
    -------
    DensityResult
        The density jump at the boundary (order ``p`` and bias-corrected), its
        jackknife standard error, z statistic and p-value, the side-wise
        densities, bandwidths and sample sizes.
    """
    if normalize not in ("none", "network", "euclidean"):
        raise ValueError("normalize must be 'euclidean', 'none' or 'network'")
    _warn_if_on_network(coords, stacklevel=3)
    q = p + 1 if q is None else int(q)
    if p < 1 or q < p:
        raise ValueError("need 1 <= p <= q")
    kfun = get_kernel(kernel)
    T = coords.treated
    om = np.ones(coords.n_obs) if weights is None else _check_weights(weights, coords.n_obs)
    b = coords.boundary
    if normalize == "euclidean":
        # straight-line distance to the boundary line, where it is reached in the plane
        s = b.distance(coords.xy)
        arc = b.project(coords.xy)
        reach = ~b.is_masked(arc)
        if not b.closed:  # beyond the ends of an open boundary the side is not defined
            reach &= (arc > 1e-9 * b.length) & (arc < b.length * (1 - 1e-9))
    else:
        s = np.abs(coords.s)
        arc = coords.crossing_arc
        reach = ~b.is_masked(np.nan_to_num(arc, nan=-1.0))
    # Masked stretches (parks, rivers) are void by construction: leave them out.
    sel = np.isfinite(s) & (om > 0) & reach
    first = np.zeros(coords.n_obs, dtype=bool)
    first[np.unique(coords.loc, return_index=True)[1]] = True
    if unique_locations:
        sel &= first
    elif (sel & ~first).any():
        warnings.warn("several observations share a location (repeated measurements?); the density test counts "
                      "each one. Use unique_locations=True to test the density of locations", stacklevel=2)
    center = None
    if eval_index is not None:
        if window is None:
            raise ValueError("give a window (arc length) together with eval_index")
        center = coords.eval_node_arc[eval_index]
        sel &= _arc_gap(coords.boundary, arc, center) <= window
    if normalize == "network":
        G_left, G_right = _street_length_transform(coords, center, window)
        s = np.where(T, G_right(s), G_left(s))
    # controls exactly on the boundary (s = 0) stay on the control side
    x = np.where(T, s, np.where(s == 0, -np.finfo(float).tiny, -s))[sel]
    om = om[sel]
    if len(x) < 20 or (x >= 0).sum() < 10 or (x < 0).sum() < 10:
        raise ValueError("too few observations for a density test")

    if h is None:
        n_min = max(10, 2 * (q + 2))
        hl = _side_bandwidth(x, om, False, p, kfun, scaleregul, n_min)
        hr = _side_bandwidth(x, om, True, p, kfun, scaleregul, n_min)
        if bwselect == "common":
            hl = hr = min(hl, hr)
        elif bwselect != "each":
            raise ValueError("bwselect must be 'common' or 'each'")
    else:
        hh = np.broadcast_to(np.asarray(h, dtype=float), (2,))
        hl, hr = float(hh[0]), float(hh[1])

    bp, _, nl, nr = _cjm(x, om, hl, hr, p, kfun)
    bq, Vq, _, _ = _cjm(x, om, hl, hr, q, kfun)
    fl, fr = bp[2], bp[3]
    fl_q, fr_q = bq[2], bq[3]
    diff = fr - fl
    diff_bc = fr_q - fl_q
    Vlr = Vq[np.ix_([2, 3], [2, 3])]
    se = float(np.sqrt(max(Vlr[0, 0] + Vlr[1, 1] - 2 * Vlr[0, 1], 0.0)))
    zstat = diff_bc / se if se > 0 else np.nan
    theta = np.log(fr / fl) if fl > 0 and fr > 0 else np.nan
    return DensityResult(
        float(diff), float(diff_bc), se, float(zstat), float(2 * stats.norm.sf(abs(zstat))), float(theta),
        float(fl), float(fr), hl, hr,
        int((x < 0).sum()), int((x >= 0).sum()), nl, nr, p, q, normalize,
    )
