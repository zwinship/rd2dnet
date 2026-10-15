"""Distance-decay kernels evaluated at ``u = d / h`` with ``d`` a network distance."""

from __future__ import annotations

from typing import Callable

import numpy as np


def triangular(u):
    u = np.abs(np.asarray(u, dtype=float))
    return np.where(u <= 1, 1 - u, 0.0)


def epanechnikov(u):
    u = np.abs(np.asarray(u, dtype=float))
    return np.where(u <= 1, 0.75 * (1 - u**2), 0.0)


def uniform(u):
    u = np.abs(np.asarray(u, dtype=float))
    return np.where(u <= 1, 0.5, 0.0)


def gaussian(u):
    u = np.asarray(u, dtype=float)
    with np.errstate(invalid="ignore", over="ignore"):
        out = np.exp(-0.5 * u**2) / np.sqrt(2 * np.pi)
    return np.where(np.isfinite(u), out, 0.0)


KERNELS: dict[str, Callable] = {
    "triangular": triangular,
    "tri": triangular,
    "epanechnikov": epanechnikov,
    "epa": epanechnikov,
    "uniform": uniform,
    "uni": uniform,
    "gaussian": gaussian,
    "gau": gaussian,
}


def get_kernel(kernel) -> Callable:
    if callable(kernel):
        return kernel
    try:
        return KERNELS[kernel.lower()]
    except (KeyError, AttributeError):
        raise ValueError(f"unknown kernel {kernel!r}; choose from {sorted(set(KERNELS))}") from None
