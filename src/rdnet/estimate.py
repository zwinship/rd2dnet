"""Network-constrained boundary average treatment effect curves (BATEC).

At each boundary evaluation point ``b`` the estimator fits, separately on the
treated and control sides, a weighted local polynomial of order ``p`` in the
network coordinates ``(s_i, o_i(b))`` with kernel weights ``K(d_i(b) / h)``
(times optional intensity weights). ``tau(b)`` is the difference of the two
intercepts; in fuzzy designs it is the ratio of the outcome and take-up jumps.
Robust bias-corrected inference uses the order ``q = p + 1`` fit, as in rd2d.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from ._wls import VCE_TYPES, poly_basis, poly_exponents, variance, wls
from .bandwidth import default_min_obs, rdbwnet
from .coordinates import NetworkCoordinates
from .kernels import get_kernel


@dataclass
class RDNetResult:
    """Pointwise BATEC estimates along the boundary.

    ``table`` has one row per evaluation point with the order-``p`` estimate
    ``tau`` and its standard error ``se``, the bias-corrected ``tau_bc`` and
    robust ``se_rb``, the robust confidence interval, and (if requested) the
    uniform confidence band ``cb_lo``/``cb_hi``. Fuzzy designs add the
    intention-to-treat jump ``itt`` and first stage ``fs`` (with ``fs_se``).
    """

    table: pd.DataFrame
    p: int
    q: int
    kernel: str
    vce: str
    level: float
    metric: str
    fuzzy: bool = False
    cb_crit: float | None = None
    cb_cov: np.ndarray | None = field(default=None, repr=False)
    bandwidths: pd.DataFrame | None = None
    coords: NetworkCoordinates | None = field(default=None, repr=False)

    def summary(self) -> str:
        t = self.table
        lines = [
            f"rdnet: network boundary discontinuity estimates ({'fuzzy' if self.fuzzy else 'sharp'})",
            f"  metric={self.metric}  kernel={self.kernel}  p={self.p}  q={self.q}  vce={self.vce}",
            f"  evaluation points: {len(t)}  estimated: {int((t['status'] == 'ok').sum())}",
        ]
        if self.cb_crit is not None:
            lines.append(f"  uniform band critical value: {self.cb_crit:.3f}")
        cols = ["arc", "h0", "h1", "n0", "n1", "tau", "se", "tau_bc", "ci_lo", "ci_hi", "pvalue", "status"]
        if self.fuzzy:
            cols.insert(5, "fs")
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


def _bandwidth_matrix(h, J: int) -> np.ndarray:
    """User bandwidths as a (J, 2) array of (control, treated)."""
    h = np.asarray(h, dtype=float)
    if h.ndim == 0:
        return np.full((J, 2), float(h))
    if h.shape == (2,) and J != 2:
        return np.tile(h, (J, 1))
    if h.shape == (J,):
        return np.column_stack([h, h])
    if h.shape == (J, 2):
        return h.copy()
    raise ValueError("h must be a scalar, (h0, h1), one value per evaluation point, or a (J, 2) array")


def _deriv_vector(order: int, deriv, tangvec_j):
    """Target vector over the order-``order`` basis, and the derivative order it picks."""
    exps = poly_exponents(order)
    vec = np.zeros(len(exps))
    if tangvec_j is not None:
        vec[exps.index((1, 0))], vec[exps.index((0, 1))] = tangvec_j
        return vec, 1
    a, c = deriv
    vec[exps.index((a, c))] = math.factorial(a) * math.factorial(c)
    return vec, a + c


def _design(us, uo, Ti, order, Zi, pooled):
    """Joint local design: side-specific polynomials (or pooled slopes) plus common covariates."""
    R = poly_basis(us, uo, order)
    if pooled:
        X = np.column_stack([np.ones(len(us)), Ti, R[:, 1:]])
    else:
        X = np.column_stack([R * Ti[:, None], R * (1.0 - Ti)[:, None]])
    return X if Zi is None else np.column_stack([X, Zi])


def rdnet(
    y,
    coords: NetworkCoordinates,
    *,
    h=None,
    p: int = 1,
    q: int | None = None,
    deriv=(0, 0),
    tangvec=None,
    covs=None,
    kernel="triangular",
    bwselect: str = "mserd",
    vce: str = "hc1",
    cluster=None,
    conley_cutoff: float | None = None,
    conley_distance: str = "network",
    conley_resid_scale: float = 3.0,
    conley_resid_order: int | None = None,
    fuzzy=None,
    weights=None,
    pooled: bool = False,
    level: float = 95.0,
    cbands: bool = True,
    repp: int = 2000,
    min_obs: int | None = None,
    scaleregul: float = 3.0,
    scalebiascrct: float = 1.0,
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
        Bandwidth in metric units: scalar, ``(h_control, h_treated)``, one value
        per evaluation point, a ``(J, 2)`` array, or ``None`` to select with
        :func:`rdbwnet` using ``bwselect`` (``"mserd"``, ``"msetwo"``,
        ``"imserd"``, ``"imsetwo"``, or the coverage-error-optimal ``"cerrd"``,
        ``"certwo"``, ``"icerrd"``, ``"icertwo"``).
    p, q
        Polynomial order for estimation and for bias-corrected inference (default ``p + 1``).
    deriv
        ``(a, c)``: estimate the jump in ``d^(a+c) mu / ds^a do^c`` instead of the
        level. ``(1, 0)`` is the jump in the slope across the boundary (a kink),
        ``(0, 1)`` the jump in the slope along it. Requires ``a + c <= p``.
    tangvec
        Direction ``(t_s, t_o)`` (or one per evaluation point, shape ``(J, 2)``):
        estimate the jump in the directional derivative ``t_s d/ds + t_o d/do``.
    covs
        Covariates (``n x m``) entering the local fit with coefficients common to
        both sides (rd2d's ``covs_eff``). Locally constant covariates are dropped.
    kernel
        ``"triangular"``, ``"epanechnikov"``, ``"uniform"``, ``"gaussian"`` or a callable.
    vce
        ``"hc0"``-``"hc3"``; ``"cluster"`` (with ``cluster``) or ``"conley"`` (with
        ``conley_cutoff``) are implied when those are given.
    conley_cutoff, conley_distance
        Spatial HAC with Bartlett weights over network (default) or Euclidean
        distance. Used for pointwise errors and for the uniform band.
    conley_resid_scale, conley_resid_order
        Conley residuals come from an auxiliary fit of order ``conley_resid_order``
        (default ``q``) at ``conley_resid_scale * h`` (default 3). Residuals from the
        estimation fit itself absorb spatially correlated shocks at scales below
        ``h``, which makes the HAC variance far too small.
    fuzzy
        Treatment take-up ``D``. The estimand becomes the ratio of the outcome
        jump to the take-up jump at each boundary point.
    weights
        Non-negative intensity weights (e.g. floor area) multiplying the kernel.
    pooled
        Fit one regression ``y ~ 1 + T + poly(s, o)`` with common slopes instead
        of separate side-specific polynomials (levels only).
    """
    y = np.asarray(y, dtype=float)
    if y.shape != (coords.n_obs,):
        raise ValueError("y must have one value per observation")
    q = p + 1 if q is None else int(q)
    if p < 0 or q < p:
        raise ValueError("need 0 <= p <= q")
    J = coords.n_eval
    deriv = tuple(int(v) for v in deriv)
    if len(deriv) != 2 or min(deriv) < 0 or sum(deriv) > p:
        raise ValueError("deriv must be (a, c) with a, c >= 0 and a + c <= p")
    tang = None
    if tangvec is not None:
        if sum(deriv) != 0:
            raise ValueError("give either deriv or tangvec, not both")
        if p < 1:
            raise ValueError("tangvec needs p >= 1")
        tang = np.broadcast_to(np.asarray(tangvec, dtype=float), (J, 2))
    if pooled and (sum(deriv) or tang is not None):
        raise ValueError("pooled=True estimates level jumps only")
    kfun = get_kernel(kernel)
    kname = kernel if isinstance(kernel, str) else getattr(kernel, "__name__", "custom")
    intensity = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    if np.any(intensity < 0):
        raise ValueError("weights must be non-negative")
    D = None
    if fuzzy is not None:
        D = np.asarray(fuzzy, dtype=float)
        if D.shape != y.shape:
            raise ValueError("fuzzy must have one value per observation")
    Z = None
    if covs is not None:
        Z = np.asarray(covs, dtype=float)
        if Z.ndim == 1:
            Z = Z[:, None]
        if Z.shape[0] != len(y):
            raise ValueError("covs must have one row per observation")
    if cluster is not None:
        vce = "cluster"
        cluster = np.asarray(cluster)
    if conley_cutoff is not None:
        vce = "conley"
    if vce not in VCE_TYPES:
        raise ValueError(f"vce must be one of {VCE_TYPES}")
    if vce == "cluster" and cluster is None:
        raise ValueError("vce='cluster' needs a cluster variable")
    if vce == "conley" and (conley_cutoff is None or conley_cutoff <= 0):
        raise ValueError("vce='conley' needs a positive conley_cutoff")
    if conley_distance not in ("network", "euclidean"):
        raise ValueError("conley_distance must be 'network' or 'euclidean'")
    base_vce = vce if vce in ("hc1", "hc2", "hc3") else "hc0"
    min_obs = default_min_obs(q) if min_obs is None else int(min_obs)
    aux_order = q if conley_resid_order is None else int(conley_resid_order)

    bw = None
    if h is None:
        bw = rdbwnet(y, coords, p=p, deriv=deriv, tangvec=tangvec, covs=Z, kernel=kernel, bwselect=bwselect,
                     vce=base_vce, cluster=cluster, fuzzy=D, weights=weights, bwcheck=min_obs,
                     scaleregul=scaleregul, scalebiascrct=scalebiascrct)
        H = bw[["h0", "h1"]].to_numpy(dtype=float)
    else:
        H = _bandwidth_matrix(h, J)

    z = stats.norm.ppf(0.5 + level / 200)
    T = coords.treated
    Tf = T.astype(float)
    s = coords.s
    finite = np.isfinite(s) & np.isfinite(y)
    if D is not None:
        finite &= np.isfinite(D)
    if Z is not None:
        finite &= np.all(np.isfinite(Z), axis=1)
    targets = [y] if D is None else [y, D]
    rows = []
    psi_store: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for j in range(J):
        rec = dict(eval=j, arc=coords.eval_arc[j], node_arc=coords.eval_node_arc[j], h0=H[j, 0], h1=H[j, 1],
                   n0=0, n1=0, tau=np.nan, se=np.nan, tau_bc=np.nan, se_rb=np.nan,
                   ci_lo=np.nan, ci_hi=np.nan, pvalue=np.nan, status="ok")
        if D is not None:
            rec.update(itt=np.nan, fs=np.nan, fs_se=np.nan)
        rows.append(rec)
        if coords.eval_masked[j]:
            rec["status"] = "masked"
            continue
        if not (np.all(np.isfinite(H[j])) and np.all(H[j] > 0)):
            rec["status"] = "no bandwidth"
            continue
        d, o = coords.at(j)
        base = np.isfinite(d) & np.isfinite(o) & finite & (intensity > 0)
        hside = np.where(T, H[j, 1], H[j, 0])
        w = kfun(d / hside) * intensity
        use = (w > 0) & base
        i1, i0 = np.flatnonzero(use & T), np.flatnonzero(use & ~T)
        rec["n1"], rec["n0"] = len(i1), len(i0)
        if min(len(i1), len(i0)) < min_obs:
            rec["status"] = "insufficient"
            continue
        idx = np.r_[i1, i0]
        zcols = None
        if Z is not None:
            zcols = np.ptp(Z[idx], axis=0) > 1e-12  # drop locally constant covariates
            if not zcols.any():
                zcols = None
        Zi = None if zcols is None else Z[idx][:, zcols]

        aux = None
        if vce == "conley":
            hs = hside * conley_resid_scale
            wk = kfun(d / hs) * intensity
            m = np.flatnonzero((wk > 0) & base)
            Xa = _design(s[m] / hs[m], o[m] / hs[m], Tf[m], aux_order,
                         None if zcols is None else Z[m][:, zcols], False)
            Xi = _design(s[idx] / hs[idx], o[idx] / hs[idx], Tf[idx], aux_order, Zi, False)
            aux = (Xa, wk[m], m, Xi)

        try:
            out = {}
            for order in (p, q):
                X = _design(s[idx] / hside[idx], o[idx] / hside[idx], Tf[idx], order, Zi, pooled)
                c = np.zeros(X.shape[1])
                if pooled:
                    c[1] = 1.0
                else:
                    vec, v = _deriv_vector(order, deriv, None if tang is None else tang[j])
                    k = len(vec)
                    c[:k], c[k:2 * k] = vec / H[j, 1] ** v, -vec / H[j, 0] ** v
                jumps, psis = [], []
                for target in targets:
                    f = wls(X, target[idx], w[idx])
                    jumps.append(float(c @ f.beta))
                    if aux is None:
                        psis.append(f.influence(c, base_vce))
                    else:
                        Xa, wa, m, Xi = aux
                        r_aux = target[idx] - Xi @ wls(Xa, target[m], wa).beta
                        psis.append(f.loadings(c) * r_aux)
                if D is None:
                    tau, psi, extra = jumps[0], psis[0], {}
                else:
                    itt, fs = jumps
                    if abs(fs) < 1e-8:
                        raise ZeroDivisionError
                    tau = itt / fs
                    psi = (psis[0] - tau * psis[1]) / fs
                    extra = {"itt": itt, "fs": fs, "_psi_fs": psis[1]}
                if vce == "conley":
                    var = float(coords.hac_quadratic(idx, psi, conley_cutoff, conley_distance)[0, 0])
                else:
                    var = variance(psi, vce, cluster=None if cluster is None else cluster[idx])
                out[order] = (tau, np.sqrt(max(var, 0.0)), psi, extra)
        except np.linalg.LinAlgError:
            rec["status"] = "singular"
            continue
        except ZeroDivisionError:
            rec["status"] = "weak first stage"
            continue
        tau_p, se_p, _, extra = out[p]
        tau_q, se_q, psi_q, _ = out[q]
        if "_psi_fs" in extra:
            psi_fs = extra.pop("_psi_fs")
            if vce == "conley":
                v_fs = float(coords.hac_quadratic(idx, psi_fs, conley_cutoff, conley_distance)[0, 0])
            else:
                v_fs = variance(psi_fs, vce, cluster=None if cluster is None else cluster[idx])
            extra["fs_se"] = np.sqrt(max(v_fs, 0.0))
        rec.update(tau=tau_p, se=se_p, tau_bc=tau_q, se_rb=se_q,
                   ci_lo=tau_q - z * se_q, ci_hi=tau_q + z * se_q,
                   pvalue=2 * stats.norm.sf(abs(tau_q / se_q)) if se_q > 0 else np.nan, **extra)
        psi_store[j] = (idx, psi_q)

    table = pd.DataFrame(rows)
    table["h"] = np.where(np.isclose(table["h0"], table["h1"]), table["h0"], np.nan)
    crit = cov = None
    if cbands and len(psi_store) >= 2:
        cov = _joint_covariance(psi_store, coords, vce, cluster, conley_cutoff, conley_distance)
        crit = _uniform_critical_value(cov, level, repp, seed)
        table["cb_lo"] = table["tau_bc"] - crit * table["se_rb"]
        table["cb_hi"] = table["tau_bc"] + crit * table["se_rb"]
    return RDNetResult(table, p, q, kname, vce, level, repr(coords.network.metric), D is not None,
                       crit, cov, bw, coords)


def _joint_covariance(psi_store, coords, vce, cluster, cutoff, distance) -> np.ndarray:
    """Covariance of the bias-corrected estimates across evaluation points."""
    keys = sorted(psi_store)
    union = np.unique(np.concatenate([psi_store[j][0] for j in keys]))
    pos = np.full(coords.n_obs, -1)
    pos[union] = np.arange(len(union))
    Psi = np.zeros((len(union), len(keys)))
    for c, j in enumerate(keys):
        idx, psi = psi_store[j]
        Psi[pos[idx], c] = psi
    if vce == "conley":
        return coords.hac_quadratic(union, Psi, cutoff, distance)
    if vce == "cluster":
        _, codes = np.unique(cluster[union], return_inverse=True)
        G = np.zeros((codes.max() + 1, len(keys)))
        np.add.at(G, codes, Psi)
        g = G.shape[0]
        return G.T @ G * (g / (g - 1) if g > 1 else 1.0)
    return Psi.T @ Psi


def _uniform_critical_value(Sigma, level, repp, seed) -> float:
    """Critical value for a uniform band: quantile of max |Z| with Z ~ N(0, corr)."""
    sd = np.sqrt(np.clip(np.diag(Sigma), 1e-300, None))
    corr = Sigma / np.outer(sd, sd)
    evals, evecs = np.linalg.eigh(corr)
    L = evecs * np.sqrt(np.clip(evals, 0, None))
    rng = np.random.default_rng(seed)
    draws = rng.standard_normal((repp, len(sd))) @ L.T
    return float(np.quantile(np.abs(draws).max(axis=1), level / 100))
