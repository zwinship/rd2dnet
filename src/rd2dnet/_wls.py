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
