"""Data-driven bandwidth selection for network boundary discontinuity estimation.

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

import numpy as np
import pandas as pd

from ._wls import poly_basis, poly_exponents
from .coordinates import NetworkCoordinates
from .kernels import get_kernel

BWSELECT = ("mserd", "msetwo", "imserd", "imsetwo")

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


def _n_terms(p: int) -> int:
    return (p + 1) * (p + 2) // 2


def _inv_gram(A: np.ndarray) -> np.ndarray:
    """Inverse of a local Gram matrix; raises when the local design is degenerate."""
    if not np.all(np.isfinite(A)) or np.linalg.cond(A) > 1e12:
        raise np.linalg.LinAlgError("local design matrix is singular")
    L = np.linalg.cholesky(A)
    Li = np.linalg.inv(L)
    return Li.T @ Li


def _rot(s: np.ndarray, o: np.ndarray, kname: str, n_eff: int) -> float:
    """rd2d's bivariate normal-reference rule of thumb."""
    mu2, l2 = _ROT_CONST.get(_ALIASES.get(kname, kname), _ROT_CONST["triangular"])
    cov = np.cov(np.vstack([s, o]), ddof=1)
    inv = np.linalg.pinv(cov, rcond=1e-20)
    det = np.linalg.det(cov)
    sqrt_det = np.sqrt(det) if np.isfinite(det) and det > 0 else np.sqrt(max(np.prod(np.linalg.eigvalsh(cov)), 1e-300))
    trace_const = 1.0 / (2.0**4 * np.pi * sqrt_det) * (2.0 * np.trace(inv @ inv) + np.trace(inv) ** 2)
    return float(((2.0 * l2) / (n_eff * mu2 * trace_const)) ** (1.0 / 6.0))


class _Side:
    """Local data on one side of the boundary at one evaluation point."""

    def __init__(self, s, o, d, y, om, kfun, cluster=None):
        self.s, self.o, self.d, self.y, self.om, self.kfun = s, o, d, y, om, kfun
        self.cluster = cluster

    def design(self, h: float, order: int, y=None):
        w = self.kfun(self.d / h) / h**2 * self.om
        keep = w > 0
        R = poly_basis(self.s[keep] / h, self.o[keep] / h, order)
        yy = (self.y if y is None else y)[keep]
        cl = None if self.cluster is None else self.cluster[keep]
        return R, w[keep], yy, cl

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
            return resd / (np.sqrt(np.maximum(1 - hii, 1e-12)) if vce == "hc2" else np.maximum(1 - hii, 1e-12))
        return resd

    def lm(self, h: float, order: int, vce: str = "hc1", varr: bool = False, y=None):
        """Raw-scale coefficients (and their scaled covariance constant)."""
        R, w, yy, cl = self.design(h, order, y)
        invG = _inv_gram((R.T * w) @ R)
        beta_s = invG @ ((R.T * w) @ yy)
        degs = np.array([a + c for a, c in poly_exponents(order)])
        beta = beta_s / h**degs
        cov = None
        if varr:
            resd = self._adjust(yy - R @ beta_s, R, w, invG, vce)
            cov = invG @ self._sigma(R * w[:, None], resd, h, cl) @ invG
        return beta, cov

    def coeff(self, dn: float, p: int) -> np.ndarray:
        """Projection of the degree p+1 terms onto the intercept: rd2d's ``vec_q``."""
        R, w, _, _ = self.design(dn, p + 1)
        k = _n_terms(p)
        Rp, S = R[:, :k], R[:, k:]
        invG = _inv_gram((Rp.T * w) @ Rp)
        tail = invG[0] @ ((Rp.T * w) @ S)
        return np.r_[np.zeros(k), tail]

    def consts(self, p: int, vec: np.ndarray, dn: float, bn1: float, bn2: float | None, vce: str, y=None) -> dict:
        R, w, yy, cl = self.design(dn, p + (2 if bn2 is not None else 1), y)
        kp, kp1 = _n_terms(p), _n_terms(p + 1)
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


def _kth(d: np.ndarray, k: int) -> float:
    k = int(min(max(k, 1), d.size))
    return float(np.partition(d, k - 1)[k - 1])


def rdbwnet(
    y,
    coords: NetworkCoordinates,
    *,
    p: int = 1,
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

    Returns one row per evaluation point with the control/treated bandwidths
    ``h0``/``h1``, the pilots ``dn``, ``bn0``, ``bn1``, the MSE constants and the
    effective sample sizes.
    """
    bwselect = bwselect.lower()
    if bwselect not in BWSELECT:
        raise ValueError(f"bwselect must be one of {BWSELECT}")
    if vce not in ("hc0", "hc1", "hc2", "hc3"):
        vce = "hc1"
    y = np.asarray(y, dtype=float)
    kfun = get_kernel(kernel)
    kname = kernel.lower() if isinstance(kernel, str) else "triangular"
    om = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    D = None if fuzzy is None else np.asarray(fuzzy, dtype=float)
    cl = None if cluster is None else np.asarray(cluster)
    bwcheck = default_min_obs(p + 1) if bwcheck is None else int(bwcheck)
    T, s = coords.treated, coords.s
    n_eff = len(np.unique(coords.xy[np.isfinite(s)], axis=0))
    e0 = np.r_[1.0, np.zeros(_n_terms(p) - 1)]
    expo_b = 1.0 / (2 * p + 6)

    rows = []
    for j in range(coords.n_eval):
        rec = {"eval": j, "arc": coords.eval_arc[j], "status": "ok"}
        rows.append(rec)
        if coords.eval_masked[j]:
            rec["status"] = "masked"
            continue
        d, o = coords.at(j)
        ok = np.isfinite(d) & np.isfinite(o) & np.isfinite(s) & np.isfinite(y) & (om > 0)
        if D is not None:
            ok &= np.isfinite(D)
        masks = [ok & ~T, ok & T]
        if min(m.sum() for m in masks) < bwcheck:
            rec["status"] = "insufficient"
            continue
        sides = [_Side(s[m], o[m], d[m], y[m], om[m], kfun, None if cl is None else cl[m]) for m in masks]
        dn = _rot(s[ok], o[ok], kname, n_eff)
        bounds = [(_kth(sd.d, bwcheck), float(sd.d.max())) for sd in sides]
        dns = [min(max(dn, lo), hi) for lo, hi in bounds]
        thr = [float(np.median(sd.d)) for sd in sides]
        try:
            ybw = [None, None]
            if D is not None:
                itt = fs = 0.0
                for t, sd in enumerate(sides):
                    sgn = 1.0 if t == 1 else -1.0
                    itt += sgn * sd.lm(dns[t], p)[0][0]
                    fs += sgn * sd.lm(dns[t], p, y=D[masks[t]])[0][0]
                if np.isfinite(fs) and abs(fs) > 1e-8:
                    ybw = [y[m] / fs - itt / fs**2 * D[m] for m in masks]
            hc, bn = [], []
            for t, sd in enumerate(sides):
                vq = sd.coeff(dns[t], p)
                bc = sd.consts(p + 1, vq, dns[t], thr[t], None, vce, y=ybw[t])
                if not (bc["V"] > 0 and bc["Reg1"] >= 0 and np.isfinite(bc["B"])):
                    raise np.linalg.LinAlgError("invalid bias-pilot constants")
                b_t = ((2 + 2 * (p + 1)) * bc["V"] / (2 * (bc["B"] ** 2 + scaleregul * bc["Reg1"]))) ** expo_b
                b_t = min(max(b_t, bounds[t][0]), bounds[t][1])
                bn.append(b_t)
                hc.append(sd.consts(p, e0, dns[t], b_t, thr[t], vce, y=ybw[t]))
                if not (hc[-1]["V"] > 0 and hc[-1]["Reg1"] >= 0 and np.isfinite(hc[-1]["B"])):
                    raise np.linalg.LinAlgError("invalid MSE constants")
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            rec["status"] = "singular"
            continue
        rec.update(dn=dn, bn0=bn[0], bn1=bn[1], n0=int(masks[0].sum()), n1=int(masks[1].sum()),
                   hmin0=bounds[0][0], hmin1=bounds[1][0], hmax0=bounds[0][1], hmax1=bounds[1][1])
        for t in (0, 1):
            for k, v in hc[t].items():
                rec[f"{k}{t}"] = v

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
    if bwselect in ("mserd", "imserd"):
        num = 2 * (g["V0"] + g["V1"])
        den = (2 * p + 2) * ((Bc[0] - Bc[1]) ** 2 + scaleregul * (g["Reg10"] + g["Reg11"]))
        h = (num.mean() / den.mean()) ** expo if bwselect == "imserd" else (num / den) ** expo
        h = np.clip(h, np.maximum(g["hmin0"], g["hmin1"]), np.maximum(g["hmax0"], g["hmax1"]))
        out.loc[good, "h0"] = out.loc[good, "h1"] = h
    else:
        for t in (0, 1):
            num = 2 * g[f"V{t}"]
            den = (2 * p + 2) * (Bc[t] ** 2 + scaleregul * g[f"Reg1{t}"])
            h = (num.mean() / den.mean()) ** expo if bwselect == "imsetwo" else (num / den) ** expo
            out.loc[good, f"h{t}"] = np.clip(h, g[f"hmin{t}"], g[f"hmax{t}"])
    out["h"] = np.where(np.isclose(out["h0"], out["h1"]), out["h0"], np.nan)
    return out
