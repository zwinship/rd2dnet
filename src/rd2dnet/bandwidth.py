"""Data-driven bandwidth selection for network boundary discontinuity estimation.

Ported from rd2d (Cattaneo, Titiunik and Yu; GPL-3.0; ``rd2d/location.py``,
``rdbw2d``) and adapted to network coordinates.

A port of rd2d's direct plug-in (``rdbw2d(method="dpi")``) to network
coordinates. Kernel weights are ``K(d_i(b)/h) / h**2`` with ``d`` the network
distance (rd2d's radial kernel) and the local basis is in ``(s/h, o/h)``.
At each evaluation point and on each side:

1. A rule-of-thumb pilot ``dn`` from a bivariate normal reference for
   ``(s, o(b))``, clipped to ``[bwcheck-th nearest, farthest]`` network distance.
2. The variance constant ``V`` from the order-``p`` fit at ``dn``.
3. A DPI pilot ``bn`` for the leading bias: the order ``p+1`` coefficients are
   estimated at ``bn``, itself chosen by the same MSE rule with its bias taken
   from an order ``p+2`` fit at the median distance.
4. The bias constant ``B``, its variance ``Reg1`` (regularisation, weighted by
   ``scaleregul``) and the next-order bias ``Reg2`` (weighted by
   ``scalebiascrct``) from an order ``p+2`` fit at the median distance.

Then, with bias ``~ h**(p+1)`` and variance ``~ h**-2``::

    h = (2 V / ((2p+2) * ((B + Reg2)**2 + scaleregul * Reg1))) ** (1 / (2p+4))

``mserd`` uses one bandwidth for both sides (constants summed / differenced),
``msetwo`` one per side, and ``imserd``/``imsetwo`` average the constants over
all evaluation points. Fuzzy designs use the linearised outcome
``y / tau_fs - tau_itt / tau_fs**2 * D`` as in rd2d.
"""

from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from ._parallel import n_threads
from ._wls import poly_basis, poly_exponents
from .coordinates import NetworkCoordinates
from .kernels import get_kernel

BWSELECT = ("mserd", "msetwo", "imserd", "imsetwo", "cerrd", "certwo", "icerrd", "icertwo")

# (mu_2(K)^2, R(K)) style constants used by rd2d's rule of thumb.
_ROT_CONST = {
    "triangular": (3.0 / 20.0, 3.0 / (2.0 * np.pi)),
    "epanechnikov": (1.0 / 6.0, 4.0 / (3.0 * np.pi)),
    "uniform": (1.0 / 4.0, 1.0 / np.pi),
    "gaussian": (1.0, 1.0 / (4.0 * np.pi)),
}
_ALIASES = {"tri": "triangular", "epa": "epanechnikov", "uni": "uniform", "gau": "gaussian"}


def default_min_obs(q: int) -> int:
    return 20 + len(poly_exponents(q))


def _n_terms(p: int, dim: int = 2) -> int:
    return p + 1 if dim == 1 else (p + 1) * (p + 2) // 2


def _inv_gram(A: np.ndarray) -> np.ndarray:
    """Inverse of a local Gram matrix; raises when the local design is degenerate."""
    if not np.all(np.isfinite(A)) or np.linalg.cond(A) > 1e12:
        raise np.linalg.LinAlgError("local design matrix is singular")
    L = np.linalg.cholesky(A)
    Li = np.linalg.inv(L)
    return Li.T @ Li


def _rot(s: np.ndarray, o: np.ndarray, kname: str, n_eff: int, dim: int = 2) -> float:
    """rd2d's bivariate normal-reference rule of thumb (a univariate one with ``dim=1``)."""
    if dim == 1:
        return float(1.06 * np.std(s, ddof=1) * n_eff ** (-1.0 / 5.0))
    mu2, l2 = _ROT_CONST.get(_ALIASES.get(kname, kname), _ROT_CONST["triangular"])
    cov = np.cov(np.vstack([s, o]), ddof=1)
    inv = np.linalg.pinv(cov, rcond=1e-20)
    det = np.linalg.det(cov)
    sqrt_det = np.sqrt(det) if np.isfinite(det) and det > 0 else np.sqrt(max(np.prod(np.linalg.eigvalsh(cov)), 1e-300))
    trace_const = 1.0 / (2.0**4 * np.pi * sqrt_det) * (2.0 * np.trace(inv @ inv) + np.trace(inv) ** 2)
    return float(((2.0 * l2) / (n_eff * mu2 * trace_const)) ** (1.0 / 6.0))


class _Side:
    """Local data on one side of the boundary at one evaluation point."""

    def __init__(self, s, o, d, y, om, kfun, cluster=None, dim=2):
        self.s, self.o, self.d, self.y, self.om, self.kfun = s, o, d, y, om, kfun
        self.dim = dim
        self.cluster = cluster
        self._basis: dict[float, tuple] = {}

    def design(self, h: float, order: int, y=None):
        # Kernel weights and the basis depend only on h; the order-``order``
        # basis is a column prefix of any higher-order one, so cache per h.
        hit = self._basis.get(h)
        if hit is None or hit[3] < order:
            w = self.kfun(self.d / h) / h**2 * self.om
            keep = np.flatnonzero(w > 0)
            hit = (keep, w[keep], poly_basis(self.s[keep] / h, self.o[keep] / h, order, self.dim), order)
            self._basis[h] = hit
        keep, w, R, _ = hit
        yy = (self.y if y is None else y)[keep]
        cl = None if self.cluster is None else self.cluster[keep]
        return R[:, :_n_terms(order, self.dim)], w, yy, cl

    @staticmethod
    def _sigma(WR, resd, h, cl):
        scores = WR * resd[:, None]
        if cl is None:
            return scores.T @ scores * h**2
        _, codes = np.unique(cl, return_inverse=True)
        summed = np.zeros((codes.max() + 1, scores.shape[1]))
        np.add.at(summed, codes, scores)
        n, k, g = len(cl), scores.shape[1], summed.shape[0]
        adj = ((n - 1) / (n - k)) * (g / (g - 1)) if g > 1 and n > k else 1.0
        return summed.T @ summed * h**2 * adj

    @staticmethod
    def _adjust(resd, R, w, invG, vce):
        n, k = R.shape
        if vce == "hc1" and n > k:
            return resd * np.sqrt(n / (n - k))
        if vce in ("hc2", "hc3"):
            hii = w * np.einsum("ij,jk,ik->i", R, invG, R)
            return resd / (np.sqrt(np.maximum(1 - hii, 1e-8)) if vce == "hc2" else np.maximum(1 - hii, 1e-8))
        return resd

    def lm(self, h: float, order: int, vce: str = "hc1", varr: bool = False, y=None):
        """Raw-scale coefficients (and their scaled covariance constant)."""
        R, w, yy, cl = self.design(h, order, y)
        invG = _inv_gram((R.T * w) @ R)
        beta_s = invG @ ((R.T * w) @ yy)
        degs = np.array([a + c for a, c in poly_exponents(order, self.dim)])
        beta = beta_s / h**degs
        cov = None
        if varr:
            resd = self._adjust(yy - R @ beta_s, R, w, invG, vce)
            cov = invG @ self._sigma(R * w[:, None], resd, h, cl) @ invG
        return beta, cov

    def coeff(self, dn: float, p: int, vec: np.ndarray) -> np.ndarray:
        """Projection of the degree p+1 terms onto the target ``vec``: rd2d's ``vec_q``."""
        R, w, _, _ = self.design(dn, p + 1)
        k = _n_terms(p, self.dim)
        Rp, S = R[:, :k], R[:, k:]
        invG = _inv_gram((Rp.T * w) @ Rp)
        tail = vec @ invG @ ((Rp.T * w) @ S)
        return np.r_[np.zeros(k), tail]

    def consts(self, p: int, vec: np.ndarray, dn: float, bn1: float, bn2: float | None, vce: str, y=None) -> dict:
        R, w, yy, cl = self.design(dn, p + (2 if bn2 is not None else 1), y)
        kp, kp1 = _n_terms(p, self.dim), _n_terms(p + 1, self.dim)
        Rp, S, Tm = R[:, :kp], R[:, kp:kp1], R[:, kp1:]
        WR = Rp * w[:, None]
        invG = _inv_gram(WR.T @ Rp)
        vec_q = np.r_[np.zeros(kp), vec @ invG @ (WR.T @ S)]
        beta_s = invG @ (WR.T @ yy)
        resd = self._adjust(yy - Rp @ beta_s, Rp, w, invG, vce)
        V = float(vec @ invG @ self._sigma(WR, resd, dn, cl) @ invG @ vec)

        beta1, cov1 = self.lm(bn1, p + 1, vce, varr=True, y=y)
        B = float(vec_q @ beta1)
        reg1 = float(vec_q @ cov1 @ vec_q / bn1 ** (2 + 2 * (p + 1)))
        reg2 = np.nan
        if bn2 is not None:
            vec_t = np.r_[np.zeros(kp1), vec @ invG @ (WR.T @ Tm)]
            beta2, _ = self.lm(bn2, p + 2, vce, y=y)
            reg2 = float(dn * (vec_t @ beta2))
        return {"V": V, "B": B, "Reg1": reg1, "Reg2": reg2}


def cer_factor(n: int, p: int) -> float:
    """rd2d's shrinkage from MSE-optimal to coverage-error-optimal bandwidth."""
    n = max(int(n), 1)
    return n ** (1.0 / (2.0 * p + 4.0) - 1.0 / (p + 4.0))


def _covariate_gamma(sides, Zs, dns, p, outcomes) -> np.ndarray:
    """Covariate coefficients common to both sides, partialling out side-specific polynomials."""
    A = 0.0
    b = 0.0
    for sd, Zt, dn, yt in zip(sides, Zs, dns, outcomes):
        w = sd.kfun(sd.d / dn) / dn**2 * sd.om
        k = w > 0
        R = poly_basis(sd.s[k] / dn, sd.o[k] / dn, p, sd.dim)
        W, Zk, yk = w[k], Zt[k], yt[k]
        invG = np.linalg.pinv((R.T * W) @ R, rcond=1e-12)
        uz = (R.T * W) @ Zk
        A = A + (Zk.T * W) @ Zk - uz.T @ invG @ uz
        b = b + (Zk.T * W) @ yk - uz.T @ invG @ ((R.T * W) @ yk)
    return np.linalg.lstsq(A, b, rcond=None)[0]


def _check_obs(v, n: int, name: str) -> np.ndarray:
    v = np.asarray(v, dtype=float).ravel()
    if v.shape != (n,):
        raise ValueError(f"{name} must have one value per observation ({n}), got {v.size}")
    return v


def _check_weights(weights, n: int) -> np.ndarray:
    w = _check_obs(weights, n, "weights")
    if not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("weights must be finite and non-negative")
    return w


def _check_cluster(cluster, n: int) -> np.ndarray:
    cl = np.asarray(cluster).ravel()
    if cl.shape != (n,):
        raise ValueError(f"cluster must have one value per observation ({n}), got {cl.size}")
    if pd.isna(cl).any():
        raise ValueError("cluster has missing values")
    return cl


def _per_location(d: np.ndarray, loc: np.ndarray) -> np.ndarray:
    """One distance per distinct location."""
    _, first = np.unique(loc, return_index=True)
    return d[first]


def _kth(d: np.ndarray, k: int) -> float:
    """Distance to the k-th nearest observation, nudged up so that observation gets positive kernel weight."""
    k = int(min(max(k, 1), d.size))
    return float(np.partition(d, k - 1)[k - 1]) * (1 + 1e-9)


def rdbwnet(
    y,
    coords: NetworkCoordinates,
    *,
    p: int = 1,
    deriv=(0, 0),
    tangvec=None,
    covs=None,
    kernel="triangular",
    bwselect: str = "mserd",
    vce: str = "hc1",
    cluster=None,
    fuzzy=None,
    weights=None,
    bwcheck: int | None = None,
    scaleregul: float = 1.0,
    scalebiascrct: float = 1.0,
) -> pd.DataFrame:
    """Select bandwidths (in metric units) at each evaluation point.

    ``deriv``/``tangvec`` target a derivative jump and ``covs`` are partialled
    out (with coefficients common to both sides, estimated at the pilot
    bandwidth), as in rd2d. The ``cer*`` variants shrink the MSE-optimal
    bandwidth by rd2d's coverage-error rate factor ``n^(1/(2p+4) - 1/(p+4))``.

    Parameters
    ----------
    y
        Outcome, one value per observation in ``coords``.
    coords
        Output of :func:`rd2dnet.network_coordinates`.
    p
        Polynomial order of the local fit.
    deriv, tangvec, covs, kernel, cluster, fuzzy, weights
        As in :func:`rd2dnet.rdnet`.
    bwselect
        ``"mserd"`` (one MSE-optimal bandwidth for both sides), ``"msetwo"``
        (one per side), ``"imserd"``/``"imsetwo"`` (one bandwidth for all
        evaluation points), or the coverage-error-optimal ``"cerrd"``,
        ``"certwo"``, ``"icerrd"``, ``"icertwo"``.
    vce
        ``"hc0"``-``"hc3"`` variance used inside the selector (others fall back to ``"hc1"``).
    bwcheck
        Every bandwidth is at least the distance to the ``bwcheck``-th nearest
        observation on each side (default ``20 +`` the number of order-``p + 1`` terms).
    scaleregul, scalebiascrct
        Weights on the bias-variance regularisation and the higher-order bias term.

    Returns
    -------
    pandas.DataFrame
        One row per evaluation point with the control/treated bandwidths
        ``h0``/``h1`` (and ``h`` when they are equal), the pilots ``dn``,
        ``bn0``, ``bn1``, the MSE constants (``V``, ``B``, ``Reg1``, ``Reg2`` per
        side), the bounds ``hmin``/``hmax``, local sample sizes and a ``status``.
    """
    bwselect = bwselect.lower()
    if bwselect not in BWSELECT:
        raise ValueError(f"bwselect must be one of {BWSELECT}")
    if vce not in ("hc0", "hc1", "hc2", "hc3", "cluster", "conley"):
        raise ValueError("vce must be one of 'hc0', 'hc1', 'hc2', 'hc3' (or 'cluster'/'conley', which select with 'hc1')")
    if vce not in ("hc0", "hc1", "hc2", "hc3"):
        vce = "hc1"
    y = _check_obs(y, coords.n_obs, "y")
    kfun = get_kernel(kernel)
    kname = kernel.lower() if isinstance(kernel, str) else "triangular"
    om = np.ones_like(y) if weights is None else _check_weights(weights, coords.n_obs)
    D = None if fuzzy is None else _check_obs(fuzzy, coords.n_obs, "fuzzy")
    cl = None if cluster is None else _check_cluster(cluster, coords.n_obs)
    Z = None
    if covs is not None:
        Z = np.asarray(covs, dtype=float)
        Z = Z[:, None] if Z.ndim == 1 else Z
    bwcheck = default_min_obs(p + 1) if bwcheck is None else int(bwcheck)
    T, s = coords.treated, coords.s
    # Effective sample sizes count distinct locations (mass points count once), as in rd2d.
    n_eff = len(np.unique(coords.loc[np.isfinite(s)]))
    n_eff_side = [len(np.unique(coords.loc[np.isfinite(s) & (T == t)])) for t in (False, True)]
    deriv = tuple(int(v) for v in deriv)
    if len(deriv) != 2 or min(deriv) < 0 or sum(deriv) > p:
        raise ValueError("deriv must be (a, c) with a, c >= 0 and a + c <= p")
    tang = None if tangvec is None else np.broadcast_to(np.asarray(tangvec, dtype=float), (coords.n_eval, 2))
    v_ord = 1 if tang is not None else sum(deriv)
    expo_b = 1.0 / (2 * p + 6)

    def point(j):
        rec = point_dim(j, 2)
        if rec["status"] == "singular" and deriv[1] == 0 and tang is None:
            # Few boundary crossings near this point (a river, rail line or highway
            # crossing the border): the offset o barely varies, so fit in s alone.
            rec1 = point_dim(j, 1)
            if rec1["status"] == "ok":
                return rec1
        return rec

    def point_dim(j, dim):
        rec = {"eval": j, "arc": coords.eval_arc[j], "status": "ok", "basis": "s,o" if dim == 2 else "s"}
        if coords.eval_masked[j]:
            rec["status"] = "masked"
            return rec
        exps_d = poly_exponents(p, dim)
        vec = np.zeros(len(exps_d))
        vec[exps_d.index(deriv)] = math.factorial(deriv[0]) * math.factorial(deriv[1])
        if tang is not None:
            vec[:] = 0.0
            vec[exps_d.index((1, 0))], vec[exps_d.index((0, 1))] = tang[j]
        d, o = coords.at(j)
        ok = np.isfinite(d) & np.isfinite(o) & np.isfinite(s) & np.isfinite(y) & (om > 0)
        if D is not None:
            ok &= np.isfinite(D)
        if Z is not None:
            ok &= np.all(np.isfinite(Z), axis=1)
        masks = [ok & ~T, ok & T]
        # bwcheck counts distinct locations (repeated measurements count once), as
        # rdrobust and rd2d count unique observations
        dloc = [_per_location(d[m], coords.loc[m]) for m in masks]
        if min(len(x) for x in dloc) < bwcheck:
            rec["status"] = "insufficient"
            return rec
        sides = [_Side(s[m], o[m], d[m], y[m], om[m], kfun, None if cl is None else cl[m], dim) for m in masks]
        dn = _rot(s[ok], o[ok], kname, n_eff, dim)
        bounds = [(_kth(dl, bwcheck), float(sd.d.max())) for dl, sd in zip(dloc, sides)]
        dns = [min(max(dn, lo), hi) for lo, hi in bounds]
        thr = [float(np.median(sd.d)) for sd in sides]
        try:
            Dside = [None if D is None else D[m] for m in masks]
            if Z is not None:
                zc = np.ptp(Z[ok], axis=0) > 1e-12
                if zc.any():
                    Zs = [Z[m][:, zc] for m in masks]
                    g_y = _covariate_gamma(sides, Zs, dns, p, [sd.y for sd in sides])
                    for t, sd in enumerate(sides):
                        sd.y = sd.y - Zs[t] @ g_y
                    if D is not None:
                        g_d = _covariate_gamma(sides, Zs, dns, p, Dside)
                        Dside = [Dside[t] - Zs[t] @ g_d for t in (0, 1)]
            ybw = [None, None]
            if D is not None:
                itt = fs = 0.0
                for t, sd in enumerate(sides):
                    sgn = 1.0 if t == 1 else -1.0
                    itt += sgn * float(vec @ sd.lm(dns[t], p)[0])
                    fs += sgn * float(vec @ sd.lm(dns[t], p, y=Dside[t])[0])
                if np.isfinite(fs) and abs(fs) > 1e-8:
                    ybw = [sides[t].y / fs - itt / fs**2 * Dside[t] for t in (0, 1)]
            hc, bn = [], []
            for t, sd in enumerate(sides):
                vq = sd.coeff(dns[t], p, vec)
                bc = sd.consts(p + 1, vq, dns[t], thr[t], None, vce, y=ybw[t])
                if not (bc["V"] > 0 and bc["Reg1"] >= 0 and np.isfinite(bc["B"])):
                    raise np.linalg.LinAlgError("invalid bias-pilot constants")
                b_t = ((2 + 2 * (p + 1)) * bc["V"] / (2 * (bc["B"] ** 2 + scaleregul * bc["Reg1"]))) ** expo_b
                b_t = min(max(b_t, bounds[t][0]), bounds[t][1])
                bn.append(b_t)
                hc.append(sd.consts(p, vec, dns[t], b_t, thr[t], vce, y=ybw[t]))
                if not (hc[-1]["V"] > 0 and hc[-1]["Reg1"] >= 0 and np.isfinite(hc[-1]["B"])):
                    raise np.linalg.LinAlgError("invalid MSE constants")
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            rec["status"] = "singular"
            return rec
        rec.update(dn=dn, bn0=bn[0], bn1=bn[1], n0=int(masks[0].sum()), n1=int(masks[1].sum()),
                   hmin0=bounds[0][0], hmin1=bounds[1][0], hmax0=bounds[0][1], hmax1=bounds[1][1])
        for t in (0, 1):
            for k, v in hc[t].items():
                rec[f"{k}{t}"] = v
        return rec

    # Evaluation points are independent; numpy releases the GIL on the large
    # array operations, so big samples run points in parallel threads. Beyond
    # about four threads memory bandwidth, not compute, is the limit (and each
    # thread holds O(n) temporaries), hence the cap.
    k = min(n_threads(coords.n_obs), coords.n_eval, 4)
    if k > 1:
        coords.node_rows(range(coords.n_eval))
        with ThreadPoolExecutor(k) as ex:
            rows = list(ex.map(point, range(coords.n_eval)))
    else:
        rows = [point(j) for j in range(coords.n_eval)]

    out = pd.DataFrame(rows)
    good = out["status"] == "ok"
    for col in ("h0", "h1"):
        out[col] = np.nan
    if not good.any():
        out["h"] = np.nan
        return out
    g = out.loc[good]
    Bc = [g[f"B{t}"] + scalebiascrct * g[f"Reg2{t}"].fillna(0.0) for t in (0, 1)]
    expo = 1.0 / (2 * p + 4)
    base = {"cerrd": "mserd", "certwo": "msetwo", "icerrd": "imserd", "icertwo": "imsetwo"}.get(bwselect, bwselect)
    cer = bwselect != base
    num_c, den_c = 2 + 2 * v_ord, 2 * p + 2 - 2 * v_ord
    if base in ("mserd", "imserd"):
        num = num_c * (g["V0"] + g["V1"])
        den = den_c * ((Bc[0] - Bc[1]) ** 2 + scaleregul * (g["Reg10"] + g["Reg11"]))
        h = (num.mean() / den.mean()) ** expo if base == "imserd" else (num / den) ** expo
        if cer:
            h = h * cer_factor(n_eff, p)
        h = np.clip(h, np.maximum(g["hmin0"], g["hmin1"]), np.maximum(g["hmax0"], g["hmax1"]))
        out.loc[good, "h0"] = out.loc[good, "h1"] = h
    else:
        for t in (0, 1):
            num = num_c * g[f"V{t}"]
            den = den_c * (Bc[t] ** 2 + scaleregul * g[f"Reg1{t}"])
            h = (num.mean() / den.mean()) ** expo if base == "imsetwo" else (num / den) ** expo
            if cer:
                h = h * cer_factor(n_eff_side[t], p)
            out.loc[good, f"h{t}"] = np.clip(h, g[f"hmin{t}"], g[f"hmax{t}"])
    out["h"] = np.where(np.isclose(out["h0"], out["h1"]), out["h0"], np.nan)
    return out
