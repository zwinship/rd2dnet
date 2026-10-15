"""Data-driven bandwidth selection for network boundary discontinuity estimation.

A direct plug-in rule in the spirit of rd2d's ``rdbw2d``. At each evaluation
point ``b``:

1. A pilot bandwidth ``h0`` is the distance capturing ``n_side ** ((p+1)/(p+2))``
   network neighbours on each side (at least ``4 * min_obs``).
2. An order ``p+1`` fit at ``h0`` estimates the leading bias term; an order ``p``
   fit estimates the variance.
3. Bias scales as ``h**(p+1)`` and variance as ``h**-2`` (two-dimensional
   coordinates), so the MSE-optimal bandwidth is
   ``h = (Cv / ((p+1) * Cb2)) ** (1 / (2p+4))`` with ``Cb2`` the squared bias
   constant regularised by ``scaleregul`` times its estimated variance.
``"imserd"`` uses one bandwidth for all points from the averaged constants.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._wls import poly_basis, poly_exponents, wls
from .coordinates import NetworkCoordinates
from .kernels import get_kernel


def _kth(d: np.ndarray, k: int) -> float:
    if d.size == 0:
        return np.nan
    k = int(min(max(k, 1), d.size))
    return float(np.partition(d, k - 1)[k - 1])


def default_min_obs(q: int) -> int:
    return 20 + len(poly_exponents(q))


def rdbwnet(
    y,
    coords: NetworkCoordinates,
    *,
    p: int = 1,
    kernel="triangular",
    bwselect: str = "mserd",
    weights=None,
    min_obs: int | None = None,
    scaleregul: float = 1.0,
) -> pd.DataFrame:
    """Select bandwidths (in metric units) at each evaluation point.

    Returns a DataFrame with one row per evaluation point: ``h`` (selected),
    ``h_mse`` (unconstrained plug-in), ``h_min`` (smallest bandwidth with
    ``min_obs`` observations per side), the pilot ``h0``, and the bias/variance
    constants.
    """
    if bwselect not in ("mserd", "imserd"):
        raise ValueError("bwselect must be 'mserd' or 'imserd'")
    y = np.asarray(y, dtype=float)
    kfun = get_kernel(kernel)
    intensity = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    min_obs = default_min_obs(p + 1) if min_obs is None else int(min_obs)
    T = coords.treated
    s = coords.s
    top = [i for i, (a, c) in enumerate(poly_exponents(p + 1)) if a + c == p + 1]

    rows = []
    for j in range(coords.n_eval):
        rec = {"eval": j, "arc": coords.eval_arc[j], "h0": np.nan, "h_min": np.nan, "h_mse": np.nan,
               "bias_const": np.nan, "var_const": np.nan, "status": "ok"}
        rows.append(rec)
        if coords.eval_masked[j]:
            rec["status"] = "masked"
            continue
        d, o = coords.at(j)
        ok = np.isfinite(d) & np.isfinite(o) & np.isfinite(s) & np.isfinite(y) & (intensity > 0)
        sides = [np.flatnonzero(ok & T), np.flatnonzero(ok & ~T)]
        if min(len(i) for i in sides) < min_obs:
            rec["status"] = "insufficient"
            continue
        rec["h_min"] = max(_kth(d[i], min_obs) for i in sides)
        h0 = max(_kth(d[i], max(4 * min_obs, int(np.ceil(len(i) ** ((p + 1) / (p + 2)))))) for i in sides)
        rec["h0"] = h0
        B = VB = V = 0.0
        try:
            for sign, idx in zip((1.0, -1.0), sides):
                w = kfun(d[idx] / h0) * intensity[idx]
                keep = w > 0
                idx, w = idx[keep], w[keep]
                us, uo = s[idx] / h0, o[idx] / h0
                Xq = poly_basis(us, uo, p + 1)
                fq = wls(Xq, y[idx], w)
                cov_q = fq.coef_cov("hc1")[np.ix_(top, top)]
                fp = wls(poly_basis(us, uo, p), y[idx], w)
                a = (fp.Ginv[0] @ (fp.X.T * w)) @ Xq[:, top]
                B += sign * float(a @ fq.beta[top])
                VB += float(a @ cov_q @ a)
                V += float(np.sum(fp.influence(0, "hc1") ** 2))
        except np.linalg.LinAlgError:
            rec["status"] = "singular"
            continue
        rec["bias_const"] = (B**2 + scaleregul * VB) / h0 ** (2 * (p + 1))
        rec["var_const"] = V * h0**2

    out = pd.DataFrame(rows)
    good = out["status"] == "ok"
    expo = 1.0 / (2 * p + 4)
    if bwselect == "mserd":
        out.loc[good, "h_mse"] = (out.loc[good, "var_const"] / ((p + 1) * out.loc[good, "bias_const"])) ** expo
    elif good.any():
        common = (out.loc[good, "var_const"].mean() / ((p + 1) * out.loc[good, "bias_const"].mean())) ** expo
        out.loc[good, "h_mse"] = common
    out["h"] = np.where(good, np.fmax(out["h_mse"], out["h_min"]), np.nan)
    return out
