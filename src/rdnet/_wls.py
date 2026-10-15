"""Local weighted least squares in network coordinates and its variance."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

VCE_TYPES = ("hc0", "hc1", "hc2", "hc3", "cluster", "conley")


def poly_exponents(p: int) -> list[tuple[int, int]]:
    """Exponents ``(a, c)`` of ``s**a * o**c`` for all terms of total degree <= p."""
    return [(a, k - a) for k in range(p + 1) for a in range(k, -1, -1)]


def poly_basis(s: np.ndarray, o: np.ndarray, p: int) -> np.ndarray:
    return np.column_stack([s**a * o**c for a, c in poly_exponents(p)])


@dataclass
class Fit:
    beta: np.ndarray
    Ginv: np.ndarray
    resid: np.ndarray
    leverage: np.ndarray
    X: np.ndarray
    w: np.ndarray

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
            e = e / np.sqrt(np.clip(1 - self.leverage, 1e-12, None))
        elif vce == "hc3":
            e = e / np.clip(1 - self.leverage, 1e-12, None)
        return (self.X @ self._row(coef)) * self.w * e

    def coef_cov(self, vce: str = "hc1") -> np.ndarray:
        """Sandwich covariance of all coefficients (heteroskedasticity-robust)."""
        e = self.resid
        if vce == "hc1":
            e = e * np.sqrt(self.n / max(self.n - self.k, 1))
        S = self.X * (self.w * e)[:, None]
        return self.Ginv @ (S.T @ S) @ self.Ginv


def wls(X: np.ndarray, y: np.ndarray, w: np.ndarray) -> Fit:
    XtW = X.T * w
    G = XtW @ X
    if np.linalg.cond(G) > 1e12:
        raise np.linalg.LinAlgError("local design matrix is singular")
    Ginv = np.linalg.inv(G)
    Ginv = (Ginv + Ginv.T) / 2
    beta = Ginv @ (XtW @ y)
    resid = y - X @ beta
    lev = w * np.einsum("ij,jk,ik->i", X, Ginv, X)
    return Fit(beta, Ginv, resid, lev, X, w)


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
