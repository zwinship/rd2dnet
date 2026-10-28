"""Local weighted least squares in network coordinates and its variance."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

VCE_TYPES = ("hc0", "hc1", "hc2", "hc3", "cluster", "conley")


def poly_exponents(p: int, dim: int = 2) -> list[tuple[int, int]]:
    """Exponents ``(a, c)`` of ``s**a * o**c`` for all terms of total degree <= p.

    With ``dim=1`` the basis is a polynomial in ``s`` alone (no offset terms).
    """
    if dim == 1:
        return [(a, 0) for a in range(p + 1)]
    return [(a, k - a) for k in range(p + 1) for a in range(k, -1, -1)]


def poly_basis(s: np.ndarray, o: np.ndarray, p: int, dim: int = 2) -> np.ndarray:
    """Columns ``s**a * o**c`` in :func:`poly_exponents` order; lower orders are column prefixes."""
    s = np.asarray(s, dtype=float)
    o = np.asarray(o, dtype=float)
    exps = poly_exponents(p, dim)
    out = np.empty((len(exps), len(s))).T  # column-major: each term is contiguous
    sp, op = [np.ones_like(s)], [np.ones_like(o)]
    for _ in range(p):
        sp.append(sp[-1] * s)
        op.append(op[-1] * o)
    for k, (a, c) in enumerate(exps):
        np.multiply(sp[a], op[c], out=out[:, k])
    return out


def _pinv_sym(A: np.ndarray, rcond: float = 1e-10) -> tuple[np.ndarray, int]:
    """Pseudo-inverse and rank of a symmetric PSD matrix (eigenvalues above ``rcond`` times the largest)."""
    if A.size == 0:
        return A.copy(), 0
    lam, V = np.linalg.eigh(A)
    top = max(float(lam.max()), 0.0)
    keep = lam > rcond * top if top > 0 else np.zeros(len(lam), dtype=bool)
    return (V[:, keep] / lam[keep]) @ V[:, keep].T, int(keep.sum())


@dataclass
class Fit:
    beta: np.ndarray
    Ginv: np.ndarray
    resid: np.ndarray
    X: np.ndarray
    w: np.ndarray

    @cached_property
    def leverage(self) -> np.ndarray:
        """Diagonal of the weighted hat matrix (computed only for HC2/HC3)."""
        return self.w * np.einsum("ij,jk,ik->i", self.X, self.Ginv, self.X)

    @property
    def n(self) -> int:
        return len(self.resid)

    @property
    def k(self) -> int:
        return len(self.beta)

    def _row(self, coef) -> np.ndarray:
        return self.Ginv[coef] if np.ndim(coef) == 0 else self.Ginv @ np.asarray(coef, dtype=float)

    def loadings(self, coef=0) -> np.ndarray:
        """Weights ``a_i`` with ``c' beta = sum_i a_i y_i`` (``coef`` an index or a contrast ``c``)."""
        return (self.X @ self._row(coef)) * self.w

    def influence(self, coef=0, vce: str = "hc1") -> np.ndarray:
        """Per-observation influence ``psi_i`` with ``c' beta - truth ~ sum_i psi_i``."""
        e = self.resid
        if vce == "hc1":
            e = e * np.sqrt(self.n / max(self.n - self.k, 1))
        elif vce == "hc2":
            e = e / np.sqrt(np.clip(1 - self.leverage, 1e-8, None))
        elif vce == "hc3":
            e = e / np.clip(1 - self.leverage, 1e-8, None)
        return (self.X @ self._row(coef)) * self.w * e


    def cr2(self, coef, cluster) -> tuple[np.ndarray, float]:
        """Bias-reduced cluster adjustment (CR2) and Bell-McCaffrey degrees of freedom.

        Returns ``R`` with ``psi_i = R_i * resid_i`` (the CR2 influence; the
        variance is the sum over clusters of the squared cluster sums, with no
        small-sample factor) and the Satterthwaite degrees of freedom of that
        variance. ``R_g = D_g^(-1/2) a_g`` with ``a`` the loadings and
        ``D_g = ((I - H)(I - H)')_gg`` the cluster's block, ``H = X B`` and
        ``B = (X'WX)^(-1) X'W`` (Bell and McCaffrey 2002, Pustejovsky and
        Tipton 2018, working model the identity). ``D_g - I = U M U'`` with
        ``U = [X_g, W_g X_g]``, so ``D_g^(-1/2)`` comes from a small
        eigenproblem; directions with eigenvalue 0 (leverage one) get weight 0.
        """
        X, w, Gi = self.X, self.w, self.Ginv
        k = X.shape[1]
        a = self.loadings(coef)
        _, codes = np.unique(cluster, return_inverse=True)
        G = codes.max() + 1
        Q = Gi @ ((X.T * w**2) @ X) @ Gi  # B B'
        M = np.block([[Q, -Gi], [-Gi, np.zeros((k, k))]])
        order = np.argsort(codes, kind="stable")
        bounds = np.r_[0, np.cumsum(np.bincount(codes, minlength=G))]
        R = np.empty_like(a)
        for g in range(G):
            m = order[bounds[g]:bounds[g + 1]]
            U = np.hstack([X[m], X[m] * w[m, None]])
            Uo, Rm = np.linalg.qr(U)
            lam, V = np.linalg.eigh(Rm @ M @ Rm.T)
            ev = 1.0 + lam  # eigenvalues of D_g on span(U); 1 elsewhere
            f = np.where(ev > 1e-10 * max(ev.max(), 1.0), 1.0 / np.sqrt(np.clip(ev, 1e-300, None)), 0.0) - 1.0
            am = a[m]
            R[m] = am + Uo @ (V @ (f * (V.T @ (Uo.T @ am))))
        # Bell-McCaffrey: P = Rb' (I - H)(I - H)' Rb with Rb the n x G block matrix of R
        Xr = np.zeros((G, k))
        np.add.at(Xr, codes, X * R[:, None])
        Br = np.zeros((G, k))
        np.add.at(Br, codes, X * (R * w)[:, None])
        Br = Br @ Gi
        P = np.diag(np.bincount(codes, weights=R * R, minlength=G)) - Xr @ Br.T - Br @ Xr.T + Xr @ Q @ Xr.T
        den = float(np.sum(P * P))
        df = float(np.trace(P) ** 2 / den) if den > 0 else np.inf
        return R, df


def wls(X: np.ndarray, y: np.ndarray, w: np.ndarray) -> Fit:
    XtW = X.T * w
    G = XtW @ X
    # conditioning of the column-equilibrated Gram matrix: covariates in their own
    # units (income in dollars) must not make a well-posed fit look singular
    dg = np.sqrt(np.clip(np.diag(G), 1e-300, None))
    if not np.all(np.isfinite(G)) or np.linalg.cond(G / np.outer(dg, dg)) > 1e12:
        raise np.linalg.LinAlgError("local design matrix is singular")
    D = np.outer(dg, dg)
    Ginv = np.linalg.inv(G / D) / D
    Ginv = (Ginv + Ginv.T) / 2
    beta = Ginv @ (XtW @ y)
    resid = y - X @ beta
    return Fit(beta, Ginv, resid, X, w)


def variance(psi: np.ndarray, vce: str, *, cluster: np.ndarray | None = None, K: np.ndarray | None = None) -> float:
    """Variance of ``sum(psi)`` under independence, clustering or a spatial HAC kernel."""
    if vce == "cluster":
        _, codes = np.unique(cluster, return_inverse=True)
        sums = np.bincount(codes, weights=psi)
        g = len(sums)
        return float(np.sum(sums**2) * (g / (g - 1) if g > 1 else 1.0))
    if vce == "conley":
        return float(psi @ K @ psi)
    return float(np.sum(psi**2))
