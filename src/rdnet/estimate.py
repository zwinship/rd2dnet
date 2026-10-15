"""Network-constrained boundary average treatment effect curves (BATEC).

At each boundary evaluation point ``b`` the estimator fits, separately on the
treated and control sides, a weighted local polynomial of order ``p`` in the
network coordinates ``(s_i, o_i(b))`` with kernel weights ``K(d_i(b) / h)``
(times optional intensity weights). ``tau(b)`` is the difference of the two
intercepts. Robust bias-corrected inference uses the order ``q = p + 1`` fit,
as in rd2d.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist

from ._wls import VCE_TYPES, poly_basis, variance, wls
from .bandwidth import default_min_obs, rdbwnet
from .coordinates import NetworkCoordinates
from .kernels import get_kernel

_MAX_CONLEY = 20000


@dataclass
class RDNetResult:
    """Pointwise BATEC estimates along the boundary.

    ``table`` has one row per evaluation point with the order-``p`` estimate
    ``tau`` and its standard error ``se``, the bias-corrected ``tau_bc`` and
    robust ``se_rb``, the robust confidence interval, and (if requested) the
    uniform confidence band ``cb_lo``/``cb_hi``.
    """

    table: pd.DataFrame
    p: int
    q: int
    kernel: str
    vce: str
    level: float
    metric: str
    cb_crit: float | None = None
    bandwidths: pd.DataFrame | None = None
    coords: NetworkCoordinates | None = field(default=None, repr=False)

    def summary(self) -> str:
        t = self.table
        lines = [
            "rdnet: network boundary discontinuity estimates",
            f"  metric={self.metric}  kernel={self.kernel}  p={self.p}  q={self.q}  vce={self.vce}",
            f"  evaluation points: {len(t)}  estimated: {int((t['status'] == 'ok').sum())}",
        ]
        if self.cb_crit is not None:
            lines.append(f"  uniform band critical value: {self.cb_crit:.3f}")
        cols = ["arc", "h", "n0", "n1", "tau", "se", "tau_bc", "ci_lo", "ci_hi", "pvalue", "status"]
        lines.append(t[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()

    def plot(self, ax=None, *, band: bool = True, true_tau=None):
        """Plot the BATEC curve with pointwise CIs and the uniform band."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 4))
        t = self.table
        ax.axhline(0, color="0.6", lw=0.8)
        if band and "cb_lo" in t:
            ax.fill_between(t["arc"], t["cb_lo"], t["cb_hi"], alpha=0.15, label="uniform band")
        ax.errorbar(t["arc"], t["tau_bc"], yerr=[t["tau_bc"] - t["ci_lo"], t["ci_hi"] - t["tau_bc"]],
                    fmt="none", capsize=2, lw=0.8, label=f"{self.level:g}% robust CI")
        ax.plot(t["arc"], t["tau"], "o-", ms=3, label="tau(b)")
        if true_tau is not None:
            ax.plot(t["arc"], true_tau, "k--", lw=1, label="truth")
        ax.set_xlabel("position along boundary")
        ax.set_ylabel("treatment effect")
        ax.legend(frameon=False)
        return ax


def _conley_kernel(coords: NetworkCoordinates, idx: np.ndarray, cutoff: float, distance: str) -> np.ndarray:
    if len(idx) > _MAX_CONLEY:
        raise MemoryError(f"Conley HAC with {len(idx)} local observations is too large; use vce='cluster'")
    if distance == "network":
        D = coords.pairwise(idx, cutoff=cutoff)
    elif distance == "euclidean":
        D = cdist(coords.xy[idx], coords.xy[idx])
    else:
        raise ValueError("conley_distance must be 'network' or 'euclidean'")
    return np.clip(1 - D / cutoff, 0, None)  # Bartlett weights


def rdnet(
    y,
    coords: NetworkCoordinates,
    *,
    h=None,
    p: int = 1,
    q: int | None = None,
    kernel="triangular",
    bwselect: str = "mserd",
    vce: str = "hc1",
    cluster=None,
    conley_cutoff: float | None = None,
    conley_distance: str = "network",
    weights=None,
    pooled: bool = False,
    level: float = 95.0,
    cbands: bool = True,
    repp: int = 2000,
    min_obs: int | None = None,
    seed: int | None = None,
) -> RDNetResult:
    """Estimate the network BATEC curve along the boundary.

    Parameters
    ----------
    y
        Outcome, one value per observation in ``coords``.
    coords
        Output of :func:`rdnet.network_coordinates`.
    h
        Bandwidth in metric units: scalar, one per evaluation point, or ``None``
        to select with :func:`rdbwnet` (``bwselect``).
    p, q
        Polynomial order for estimation and for bias-corrected inference (default ``p + 1``).
    kernel
        ``"triangular"``, ``"epanechnikov"``, ``"uniform"``, ``"gaussian"`` or a callable.
    vce
        ``"hc0"``–``"hc3"``; ``"cluster"`` (with ``cluster``) or ``"conley"`` (with
        ``conley_cutoff``) are implied when those are given.
    conley_cutoff, conley_distance
        Spatial HAC with Bartlett weights over network (default) or Euclidean distance.
    weights
        Non-negative intensity weights (e.g. floor area) multiplying the kernel.
    pooled
        Fit one regression ``y ~ 1 + T + poly(s, o)`` with common slopes instead
        of separate side-specific polynomials.
    """
    y = np.asarray(y, dtype=float)
    if y.shape != (coords.n_obs,):
        raise ValueError("y must have one value per observation")
    q = p + 1 if q is None else int(q)
    if p < 0 or q < p:
        raise ValueError("need 0 <= p <= q")
    kfun = get_kernel(kernel)
    kname = kernel if isinstance(kernel, str) else getattr(kernel, "__name__", "custom")
    intensity = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    if np.any(intensity < 0):
        raise ValueError("weights must be non-negative")
    if cluster is not None:
        vce = "cluster"
        cluster = np.asarray(cluster)
    if conley_cutoff is not None:
        vce = "conley"
    if vce not in VCE_TYPES:
        raise ValueError(f"vce must be one of {VCE_TYPES}")
    if vce == "cluster" and cluster is None:
        raise ValueError("vce='cluster' needs a cluster variable")
    if vce == "conley" and conley_cutoff is None:
        raise ValueError("vce='conley' needs conley_cutoff")
    base_vce = "hc1" if vce == "hc1" else ("hc2" if vce == "hc2" else ("hc3" if vce == "hc3" else "hc0"))
    min_obs = default_min_obs(q) if min_obs is None else int(min_obs)

    J = coords.n_eval
    bw = None
    if h is None:
        bw = rdbwnet(y, coords, p=p, kernel=kernel, bwselect=bwselect, weights=weights, min_obs=min_obs)
        hvec = bw["h"].to_numpy(dtype=float)
    else:
        hvec = np.broadcast_to(np.asarray(h, dtype=float), (J,)).copy()

    z = stats.norm.ppf(0.5 + level / 200)
    T = coords.treated
    s = coords.s
    rows = []
    psi_store: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for j in range(J):
        rec = dict(eval=j, arc=coords.eval_arc[j], node_arc=coords.eval_node_arc[j], h=hvec[j],
                   n0=0, n1=0, tau=np.nan, se=np.nan, tau_bc=np.nan, se_rb=np.nan,
                   ci_lo=np.nan, ci_hi=np.nan, pvalue=np.nan, status="ok")
        rows.append(rec)
        if coords.eval_masked[j]:
            rec["status"] = "masked"
            continue
        hj = hvec[j]
        if not np.isfinite(hj) or hj <= 0:
            rec["status"] = "no bandwidth"
            continue
        d, o = coords.at(j)
        w = kfun(d / hj) * intensity
        use = (w > 0) & np.isfinite(d) & np.isfinite(o) & np.isfinite(s) & np.isfinite(y)
        i1, i0 = np.flatnonzero(use & T), np.flatnonzero(use & ~T)
        rec["n1"], rec["n0"] = len(i1), len(i0)
        if min(len(i1), len(i0)) < min_obs:
            rec["status"] = "insufficient"
            continue
        idx = np.r_[i1, i0]
        us, uo = s / hj, o / hj
        K = _conley_kernel(coords, idx, conley_cutoff, conley_distance) if vce == "conley" else None
        try:
            out = {}
            for order in (p, q):
                if pooled:
                    X = np.column_stack([np.ones(len(idx)), T[idx].astype(float), poly_basis(us[idx], uo[idx], order)[:, 1:]])
                    f = wls(X, y[idx], w[idx])
                    tau, psi = f.beta[1], f.influence(1, base_vce)
                else:
                    f1 = wls(poly_basis(us[i1], uo[i1], order), y[i1], w[i1])
                    f0 = wls(poly_basis(us[i0], uo[i0], order), y[i0], w[i0])
                    tau = f1.beta[0] - f0.beta[0]
                    psi = np.r_[f1.influence(0, base_vce), -f0.influence(0, base_vce)]
                var = variance(psi, vce, cluster=None if cluster is None else cluster[idx], K=K)
                out[order] = (tau, np.sqrt(max(var, 0.0)), psi)
        except np.linalg.LinAlgError:
            rec["status"] = "singular"
            continue
        tau_p, se_p, _ = out[p]
        tau_q, se_q, psi_q = out[q]
        rec.update(tau=tau_p, se=se_p, tau_bc=tau_q, se_rb=se_q,
                   ci_lo=tau_q - z * se_q, ci_hi=tau_q + z * se_q,
                   pvalue=2 * stats.norm.sf(abs(tau_q / se_q)) if se_q > 0 else np.nan)
        psi_store[j] = (idx, psi_q)

    table = pd.DataFrame(rows)
    crit = None
    if cbands and len(psi_store) >= 2:
        crit = _uniform_critical_value(psi_store, coords.n_obs, vce, cluster, level, repp, seed)
        se = table["se_rb"]
        table["cb_lo"] = table["tau_bc"] - crit * se
        table["cb_hi"] = table["tau_bc"] + crit * se
    return RDNetResult(table, p, q, kname, vce, level, repr(coords.network.metric), crit, bw, coords)


def _uniform_critical_value(psi_store, n, vce, cluster, level, repp, seed) -> float:
    """Critical value for a uniform band from the joint influence covariance across points."""
    keys = sorted(psi_store)
    Psi = np.zeros((n, len(keys)))
    for c, j in enumerate(keys):
        idx, psi = psi_store[j]
        Psi[idx, c] = psi
    if vce == "cluster":
        _, codes = np.unique(cluster, return_inverse=True)
        G = np.zeros((codes.max() + 1, len(keys)))
        np.add.at(G, codes, Psi)
        Sigma = G.T @ G
    else:
        if vce == "conley":
            warnings.warn("uniform band uses heteroskedasticity-robust cross-point correlation under vce='conley'",
                          stacklevel=3)
        Sigma = Psi.T @ Psi
    sd = np.sqrt(np.clip(np.diag(Sigma), 1e-300, None))
    corr = Sigma / np.outer(sd, sd)
    evals, evecs = np.linalg.eigh(corr)
    L = evecs * np.sqrt(np.clip(evals, 0, None))
    rng = np.random.default_rng(seed)
    draws = rng.standard_normal((repp, len(keys))) @ L.T
    return float(np.quantile(np.abs(draws).max(axis=1), level / 100))
