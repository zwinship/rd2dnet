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
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from ._wls import VCE_TYPES, poly_basis, poly_exponents, variance, wls
from .bandwidth import _check_cluster, _check_covs, _check_obs, _check_weights, default_min_obs, rdbwnet
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

    ``cov`` is the ``J x J`` covariance of ``tau_bc`` across evaluation points
    (rd2d's ``params_cov``), with NaN rows and columns for points that were not
    estimated; its diagonal is ``se_rb**2``. ``crit`` is the uniform critical
    value it implies. Both drive the uniform band, :meth:`wbate` and :meth:`lbate`.
    """

    table: pd.DataFrame
    p: int
    q: int
    kernel: str
    vce: str
    level: float
    metric: str
    fuzzy: bool = False
    crit: float | None = None
    cov: np.ndarray | None = field(default=None, repr=False)
    bandwidths: pd.DataFrame | None = None
    coords: NetworkCoordinates | None = field(default=None, repr=False)

    def _ok(self) -> np.ndarray:
        """Boolean mask of evaluation points that were estimated."""
        return (self.table["status"] == "ok").to_numpy()

    def wbate(self, weights=None) -> pd.Series:
        """Weighted boundary average treatment effect ``w' tau`` (rd2d's WBATE).

        ``weights`` has one entry per evaluation point and is normalised to sum
        to one; by default every estimated point gets equal weight. Points that
        were not estimated (masked, insufficient, ...) must get weight 0. The standard error is ``sqrt(w' cov w)``, so it accounts
        for the correlation between nearby evaluation points.
        """
        J = len(self.table)
        w = self._ok().astype(float) if weights is None else np.asarray(weights, dtype=float)
        if w.shape != (J,):
            raise ValueError(f"WBATE weights must have one entry per evaluation point ({J})")
        if not np.all(np.isfinite(w)):
            raise ValueError("WBATE weights must be finite")
        if w.sum() == 0:
            raise ValueError("WBATE weights must have a nonzero sum")
        ok = self._ok()
        bad = np.flatnonzero(~ok & (w != 0))
        if len(bad):
            raise ValueError(f"WBATE puts weight on evaluation points that were not estimated: {bad.tolist()}; "
                             "give them weight 0")
        w = (w / w.sum())[ok]
        t = self.table.loc[ok]
        tau_q = float(w @ t["tau_bc"].to_numpy())
        se = float(np.sqrt(max(w @ self.cov[np.ix_(ok, ok)] @ w, 0.0)))
        z = stats.norm.ppf(0.5 + self.level / 200)
        return pd.Series(dict(tau=float(w @ t["tau"].to_numpy()), tau_bc=tau_q, se_rb=se,
                              ci_lo=tau_q - z * se, ci_hi=tau_q + z * se,
                              pvalue=2 * stats.norm.sf(abs(tau_q / se)) if se > 0 else np.nan), name="WBATE")

    def lbate(self) -> pd.Series:
        """Largest boundary average treatment effect ``max_j tau(b_j)`` (rd2d's LBATE).

        The confidence interval is ``[max_j (tau_bc - crit * se_rb), max_j (tau_bc + crit * se_rb)]``
        over the estimated points, with ``crit`` the uniform critical value, so
        it inherits the uniform band's joint coverage.
        """
        ok = self._ok()
        if not ok.any() or self.crit is None:
            raise ValueError("LBATE needs at least one estimated evaluation point")
        t = self.table.loc[ok]
        lo = t["tau_bc"] - self.crit * t["se_rb"]
        hi = t["tau_bc"] + self.crit * t["se_rb"]
        return pd.Series(dict(tau=t["tau"].max(), tau_bc=t["tau_bc"].max(), se_rb=np.nan,
                              ci_lo=lo.max(), ci_hi=hi.max(), pvalue=np.nan), name="LBATE")

    def summary(self, *, wbate=None, lbate: bool = False) -> str:
        """Text table of the pointwise estimates, optionally with WBATE/LBATE rows."""
        t = self.table
        lines = [
            f"rd2dnet: network boundary discontinuity estimates ({'fuzzy' if self.fuzzy else 'sharp'})",
            f"  metric={self.metric}  kernel={self.kernel}  p={self.p}  q={self.q}  vce={self.vce}",
            f"  evaluation points: {len(t)}  estimated: {int(self._ok().sum())}",
        ]
        if "cb_lo" in t:
            lines.append(f"  uniform band critical value: {self.crit:.3f}")
        cols = ["arc", "h0", "h1", "n0", "n1", "tau", "se", "tau_bc", "ci_lo", "ci_hi", "pvalue", "status"]
        if self.fuzzy:
            cols.insert(5, "fs")
        fmt = lambda v: f"{v:.3f}"  # noqa: E731
        lines.append(t[cols].to_string(index=False, float_format=fmt))
        agg = ([self.wbate(wbate)] if wbate is not None else []) + ([self.lbate()] if lbate else [])
        if agg:
            lines += ["", "aggregate effects (tau_bc and CI are bias-corrected):",
                      pd.DataFrame(agg).to_string(float_format=fmt, na_rep="")]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()

    def plot(self, ax=None, *, band: bool = True, true_tau=None):
        """Plot the BATEC curve with pointwise CIs and the uniform band.

        ``true_tau`` (optional, drawn dashed) is one value per evaluation point, a
        single value, or a function of the arc position.
        """
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
            truth = np.broadcast_to(np.asarray(true_tau(t["arc"].to_numpy()) if callable(true_tau) else true_tau,
                                               dtype=float), (len(t),))
            ax.plot(t["arc"], truth, "k--", lw=1, label="truth")
        ax.set_xlabel("position along boundary")
        ax.set_ylabel("treatment effect")
        ax.legend(frameon=False)
        return ax


def _bandwidth_matrix(h, J: int) -> np.ndarray:
    """User bandwidths as a (J, 2) array of (control, treated)."""
    pair = isinstance(h, tuple) and len(h) == 2
    h = np.asarray(h, dtype=float)
    if np.any(np.isnan(h)) or np.any(h <= 0):
        raise ValueError("bandwidths must be positive")
    if h.ndim == 0:
        return np.full((J, 2), float(h))
    if h.shape == (2,) and (pair or J != 2):
        return np.tile(h, (J, 1))
    if h.shape == (J,):
        return np.column_stack([h, h])
    if h.shape == (J, 2):
        return h.copy()
    raise ValueError("h must be a scalar, (h0, h1), one value per evaluation point, or a (J, 2) array")


def _deriv_vector(order: int, deriv, tangvec_j, dim: int = 2):
    """Target vector over the order-``order`` basis, and the derivative order it picks."""
    exps = poly_exponents(order, dim)
    vec = np.zeros(len(exps))
    if tangvec_j is not None:
        vec[exps.index((1, 0))], vec[exps.index((0, 1))] = tangvec_j
        return vec, 1
    a, c = deriv
    vec[exps.index((a, c))] = math.factorial(a) * math.factorial(c)
    return vec, a + c


def _design(us, uo, Ti, order, Zi, pooled, dim=2):
    """Joint local design: side-specific polynomials (or pooled slopes) plus common covariates."""
    R = poly_basis(us, uo, order, dim)
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
        Output of :func:`rd2dnet.network_coordinates`.
    h
        Bandwidth in metric units: scalar, a tuple ``(h_control, h_treated)``,
        one value per evaluation point (an array), a ``(J, 2)`` array, or
        ``None`` to select with
        :func:`rdbwnet` using ``bwselect`` (``"mserd"``, ``"msetwo"``,
        ``"imserd"``, ``"imsetwo"``, or the coverage-error-optimal ``"cerrd"``,
        ``"certwo"``, ``"icerrd"``, ``"icertwo"``). Selection uses
        ``scaleregul=3``, as rd2d does for estimation, so it differs from a
        call to :func:`rdbwnet` with its defaults; the selected bandwidths are
        stored in ``result.bandwidths``.
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
    cluster
        Cluster labels, one per observation; implies ``vce="cluster"``. Points
        with fewer than two clusters on a side get status ``"few clusters"``.
        With repeated measurements (panels), clusters much smaller than the
        bandwidth, such as the unit itself, can still understate the variance
        of spatially correlated shocks; ``conley_cutoff`` is the safer choice.
    level
        Confidence level in percent for intervals and the uniform band.
    cbands
        Add the uniform confidence band ``cb_lo``/``cb_hi`` to the table. The
        joint covariance and its critical value are stored either way.
    repp
        Gaussian draws used to simulate the uniform critical value.
    min_obs
        Minimum local observations per side at an evaluation point (default
        ``20 +`` the number of order-``q`` terms); points with fewer are
        reported as ``"insufficient"``.
    scaleregul, scalebiascrct
        Bandwidth-selection regularisation and higher-order bias weights,
        passed to :func:`rdbwnet` when ``h`` is selected.
    seed
        Seed for the critical-value simulation (for reproducible bands and LBATE).

    Notes
    -----
    The local fit is a polynomial in ``(s, o)``. Where the network meets the
    boundary at only a few points (rivers, rail lines, highways crossing a
    border), the offset ``o`` hardly varies near an evaluation point and that
    fit is singular; the point is then fitted in ``s`` alone and its ``basis``
    column reads ``"s"`` instead of ``"s,o"``.

    ``status`` is ``"ok"``, ``"masked"``, ``"insufficient"`` (fewer than
    ``min_obs`` observations on a side), ``"few clusters"``, ``"singular"``,
    ``"weak first stage"`` or ``"no bandwidth"``.

    Returns
    -------
    RDNetResult
        Per-evaluation-point table (``tau``, ``se``, ``tau_bc``, ``se_rb``,
        ``ci_lo``/``ci_hi``, ``pvalue``, bandwidths, local sample sizes, the
        local ``basis`` and a ``status``), the joint covariance ``cov``, the uniform critical value
        ``crit``, and the selected bandwidths. Use :meth:`RDNetResult.wbate`
        and :meth:`RDNetResult.lbate` for aggregate effects.

    Examples
    --------
    >>> import rd2dnet
    >>> from rd2dnet.datasets import load_mock_midtown
    >>> city = load_mock_midtown()
    >>> coords = rd2dnet.network_coordinates(city.network(), city.boundary, city.obs[["x", "y"]], n_eval=9)
    >>> res = rd2dnet.rdnet(city.obs["outcome"], coords, seed=0)
    >>> bool(res.wbate(res.table["status"].eq("ok").astype(float))["tau_bc"] > 0)
    True
    >>> res.table[["arc", "tau_bc", "ci_lo", "ci_hi", "status"]]  # doctest: +SKIP
    """
    y = _check_obs(y, coords.n_obs, "y")
    if not 1 < level < 100:
        raise ValueError(f"level is a percentage, e.g. 95 (got {level})")
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
    intensity = np.ones_like(y) if weights is None else _check_weights(weights, coords.n_obs)
    D = None if fuzzy is None else _check_obs(fuzzy, coords.n_obs, "fuzzy")
    Z = None
    if covs is not None:
        Z = _check_covs(covs, coords.n_obs)
    if cluster is not None:
        vce = "cluster"
        cluster = _check_cluster(cluster, coords.n_obs)
    if conley_cutoff is not None:
        if cluster is not None:
            warnings.warn("both cluster and conley_cutoff given: standard errors use Conley; cluster only "
                          "affects bandwidth selection", stacklevel=2)
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
        bases = bw["basis"].tolist() if "basis" in bw else [None] * J
    else:
        bases = [None] * J
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
    # Conley: all variances come from one HAC pass over every local sample at the end.
    hac_cols: list[tuple[int, str, np.ndarray, np.ndarray]] = []
    conley = vce == "conley"

    def local_var(idx, psi, k):
        if conley:
            return np.nan
        v = variance(psi, vce, cluster=None if cluster is None else cluster[idx])
        if vce == "cluster":  # variance() applies g / (g - 1); rd2d's CR1 also has (n - 1) / (n - k)
            v *= _cluster_dof(len(idx), k, cluster[idx])
        return v
    for j in range(J):
        rec = dict(eval=j, arc=coords.eval_arc[j], node_arc=coords.eval_node_arc[j], h0=H[j, 0], h1=H[j, 1],
                   n0=0, n1=0, tau=np.nan, se=np.nan, tau_bc=np.nan, se_rb=np.nan,
                   ci_lo=np.nan, ci_hi=np.nan, pvalue=np.nan, status="ok", basis="")
        if D is not None:
            rec.update(itt=np.nan, fs=np.nan, fs_se=np.nan)
        rows.append(rec)
        if coords.eval_masked[j]:
            rec["status"] = "masked"
            continue
        if not (np.all(np.isfinite(H[j])) and np.all(H[j] > 0)):
            # the selector's reason (insufficient, singular) when it has one
            rec["status"] = bw["status"].iloc[j] if bw is not None and bw["status"].iloc[j] != "ok" else "no bandwidth"
            continue
        d, o = coords.at(j)
        base = np.isfinite(d) & np.isfinite(o) & finite & (intensity > 0)
        hside = np.where(T, H[j, 1], H[j, 0])
        # each side's kernel weights scaled by 1 / h_side^2, as in rd2d: this matters
        # for the coefficients shared across sides (covariates, pooled slopes)
        w = kfun(d / hside) / hside**2 * intensity
        use = (w > 0) & base
        i1, i0 = np.flatnonzero(use & T), np.flatnonzero(use & ~T)
        rec["n1"], rec["n0"] = len(i1), len(i0)
        # count distinct locations: repeated measurements at one site add no design points
        if min(len(np.unique(coords.loc[i1])), len(np.unique(coords.loc[i0]))) < min_obs:
            rec["status"] = "insufficient"
            continue
        idx = np.r_[i1, i0]
        zcols = None
        if Z is not None:
            zcols = np.ptp(Z[idx], axis=0) > 1e-12  # drop locally constant covariates
            if not zcols.any():
                zcols = None
        Zi = None if zcols is None else Z[idx][:, zcols]
        if vce == "cluster":
            g_side = min(len(np.unique(cluster[i1])), len(np.unique(cluster[i0])))
            if g_side < 2:
                rec["status"] = "few clusters"
                continue

        def fit(dim):
            aux = None
            if vce == "conley":
                hs = hside * conley_resid_scale
                wk = kfun(d / hs) / hs**2 * intensity
                m = np.flatnonzero((wk > 0) & base)
                Xa = _design(s[m] / hs[m], o[m] / hs[m], Tf[m], aux_order,
                             None if zcols is None else Z[m][:, zcols], False, dim)
                Xi = _design(s[idx] / hs[idx], o[idx] / hs[idx], Tf[idx], aux_order, Zi, False, dim)
                aux = (Xa, wk[m], m, Xi)
            out, ks = {}, {}
            for order in (p, q):
                X = _design(s[idx] / hside[idx], o[idx] / hside[idx], Tf[idx], order, Zi, pooled, dim)
                c = np.zeros(X.shape[1])
                if pooled:
                    c[1] = 1.0
                else:
                    vec, v = _deriv_vector(order, deriv, None if tang is None else tang[j], dim)
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
                ks[order] = X.shape[1]
                out[order] = (tau, np.sqrt(max(local_var(idx, psi, ks[order]), 0.0)), psi, extra)
            return out, ks

        # Few boundary crossings near this point (a river, rail line or highway
        # crossing the border) leave no variation in the offset o; the fit is
        # then a polynomial in s alone.
        dims = {"s,o": (2,), "s": (1,)}.get(bases[j], (2, 1))
        if deriv[1] > 0 or tang is not None:
            dims = (2,)
        out = None
        try:
            for dim in dims:
                try:
                    out, ks = fit(dim)
                    rec["basis"] = "s,o" if dim == 2 else "s"
                    break
                except np.linalg.LinAlgError:
                    continue
        except ZeroDivisionError:
            rec["status"] = "weak first stage"
            continue
        if out is None:
            rec["status"] = "singular"
            continue
        tau_p, se_p, psi_p, extra = out[p]
        tau_q, se_q, psi_q, _ = out[q]
        if "_psi_fs" in extra:
            psi_fs = extra.pop("_psi_fs")
            extra["fs_se"] = np.sqrt(max(local_var(idx, psi_fs, ks[p]), 0.0))
            if conley:
                hac_cols.append((j, "fs_se", idx, psi_fs))
        if conley:
            hac_cols += [(j, "se", idx, psi_p), (j, "se_rb", idx, psi_q)]
        rec.update(tau=tau_p, se=se_p, tau_bc=tau_q, **extra)
        _set_rb(rec, se_q, z)
        psi_store[j] = (idx, psi_q, ks[q])
    table = pd.DataFrame(rows)
    table["h"] = np.where(np.isclose(table["h0"], table["h1"]), table["h0"], np.nan)
    crit = cov = None
    if conley and hac_cols:
        union = np.unique(np.concatenate([c[2] for c in hac_cols]))
        pos = np.full(coords.n_obs, -1)
        pos[union] = np.arange(len(union))
        Psi = np.zeros((len(union), len(hac_cols)))
        for c, (_, _, idx, psi) in enumerate(hac_cols):
            Psi[pos[idx], c] = psi
        V = coords.hac_quadratic(union, Psi, conley_cutoff, conley_distance)
        se_all = np.sqrt(np.clip(np.diag(V), 0.0, None))
        for c, (j, name, _, _) in enumerate(hac_cols):
            if name == "se_rb":
                _set_rb(rows[j], se_all[c], z)
            else:
                rows[j][name] = se_all[c]
        table = pd.DataFrame(rows)
        table["h"] = np.where(np.isclose(table["h0"], table["h1"]), table["h0"], np.nan)
        qc = [c for c, col in enumerate(hac_cols) if col[1] == "se_rb"]
        qj = [hac_cols[c][0] for c in qc]
        cov = np.full((J, J), np.nan)
        cov[np.ix_(qj, qj)] = V[np.ix_(qc, qc)]
    elif psi_store:
        cov = _joint_covariance(psi_store, coords, vce, cluster, J)
    if cov is not None:
        ok = np.isfinite(np.diag(cov))
        crit = _uniform_critical_value(cov[np.ix_(ok, ok)], level, repp, seed)
        if cbands:
            table["cb_lo"] = table["tau_bc"] - crit * table["se_rb"]
            table["cb_hi"] = table["tau_bc"] + crit * table["se_rb"]
    return RDNetResult(table, p, q, kname, vce, level, repr(coords.network.metric), fuzzy=D is not None,
                       crit=crit, cov=cov, bandwidths=bw, coords=coords)


def _cluster_dof(n: int, k: int, cl) -> float:
    """rd2d's ``(n - 1) / (n - k)`` cluster degrees-of-freedom factor (1 with one cluster or n <= k)."""
    return (n - 1) / (n - k) if n > k and len(np.unique(cl)) > 1 else 1.0


def _set_rb(rec: dict, se_q: float, z: float) -> None:
    """Robust standard error, confidence interval and p-value from the bias-corrected fit."""
    tau_q = rec["tau_bc"]
    rec.update(se_rb=se_q, ci_lo=tau_q - z * se_q, ci_hi=tau_q + z * se_q,
               pvalue=2 * stats.norm.sf(abs(tau_q / se_q)) if se_q > 0 else np.nan)


def _joint_covariance(psi_store, coords, vce, cluster, J) -> np.ndarray:
    """``J x J`` covariance of the bias-corrected estimates, NaN where a point was not estimated.

    Its diagonal reproduces the pointwise variances exactly. For clustering, the
    small-sample factor ``(n - 1) / (n - k) * g / (g - 1)`` uses each point's own
    local sample, as the pointwise standard errors do. (Conley is handled in one
    HAC pass inside :func:`rdnet`.)
    """
    keys = sorted(psi_store)
    union = np.unique(np.concatenate([psi_store[j][0] for j in keys]))
    pos = np.full(coords.n_obs, -1)
    pos[union] = np.arange(len(union))
    Psi = np.zeros((len(union), len(keys)))
    for c, j in enumerate(keys):
        idx, psi, _ = psi_store[j]
        Psi[pos[idx], c] = psi
    if vce == "cluster":
        _, codes = np.unique(cluster[union], return_inverse=True)
        G = np.zeros((codes.max() + 1, len(keys)))
        np.add.at(G, codes, Psi)
        f = np.empty(len(keys))
        for c, j in enumerate(keys):
            idx, _, k = psi_store[j]
            g = len(np.unique(cluster[idx]))
            f[c] = np.sqrt((g / (g - 1) if g > 1 else 1.0) * _cluster_dof(len(idx), k, cluster[idx]))
        V = (G.T @ G) * np.outer(f, f)
    else:
        V = Psi.T @ Psi
    out = np.full((J, J), np.nan)
    out[np.ix_(keys, keys)] = V
    return out


def _uniform_critical_value(Sigma, level, repp, seed) -> float:
    """Critical value for a uniform band: quantile of max |Z| with Z ~ N(0, corr)."""
    if Sigma.shape[0] == 1:  # one point: the normal quantile itself, no simulation
        return float(stats.norm.ppf(0.5 + level / 200))
    sd = np.sqrt(np.clip(np.diag(Sigma), 1e-300, None))
    corr = Sigma / np.outer(sd, sd)
    evals, evecs = np.linalg.eigh(corr)
    L = evecs * np.sqrt(np.clip(evals, 0, None))
    rng = np.random.default_rng(seed)
    draws = rng.standard_normal((repp, len(sd))) @ L.T
    return float(np.quantile(np.abs(draws).max(axis=1), level / 100))
