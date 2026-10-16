"""Boundary geometry: the treatment threshold, side assignment and masking."""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import shapely
from shapely.geometry import LineString
from shapely.geometry.base import BaseGeometry

from ._parallel import chunked


class Boundary:
    """A boundary line separating treated from control locations.

    Parameters
    ----------
    line
        The boundary as a shapely ``LineString`` or a sequence of ``(x, y)``
        vertices, in the same projected coordinates as the network.
    treated_region
        Polygon containing the treated side. Points inside it are treated.
    treated_side
        Alternatively ``"left"`` or ``"right"`` of the line's direction.
    mask
        Geometries (e.g. a park polygon) or ``(start, end)`` arc-length intervals
        along the line. Evaluation points touching them are skipped: use this
        for segments with no mass on one side (Central Park on 60th Street).
    tol
        Distance at which a network node counts as lying on the boundary.
    """

    def __init__(
        self,
        line,
        *,
        treated_region: BaseGeometry | None = None,
        treated_side: str | None = None,
        mask: Iterable | None = None,
        tol: float = 1e-6,
    ):
        self.line = line if isinstance(line, LineString) else LineString(line)
        if (treated_region is None) == (treated_side is None):
            raise ValueError("give exactly one of treated_region or treated_side")
        if treated_side not in (None, "left", "right"):
            raise ValueError("treated_side must be 'left' or 'right'")
        self.treated_region = treated_region
        self.treated_side = treated_side
        self.tol = float(tol)
        self.mask_geoms: list[BaseGeometry] = []
        self.mask_arcs: list[tuple[float, float]] = []
        for m in mask or ():
            if isinstance(m, BaseGeometry):
                self.mask_geoms.append(m)
            else:
                a, b = m
                self.mask_arcs.append((min(a, b), max(a, b)))
        shapely.prepare(self.line)

    def __repr__(self) -> str:
        side = f"region" if self.treated_region is not None else self.treated_side
        return f"Boundary(length={self.length:.1f}, treated={side}, masks={len(self.mask_geoms) + len(self.mask_arcs)})"

    @property
    def length(self) -> float:
        return float(self.line.length)

    @staticmethod
    def _xy(xy) -> np.ndarray:
        xy = np.asarray(xy, dtype=float)
        return xy.reshape(-1, 2)

    def distance(self, xy) -> np.ndarray:
        """Euclidean distance from each point to the boundary line."""
        xy = self._xy(xy)
        return chunked(lambda a: shapely.distance(self.line, shapely.points(a)), len(xy), xy)

    def project(self, xy) -> np.ndarray:
        """Arc-length position of each point's projection onto the line."""
        xy = self._xy(xy)
        return chunked(lambda a: shapely.line_locate_point(self.line, shapely.points(a)), len(xy), xy)

    def interpolate(self, arc) -> np.ndarray:
        """(x, y) coordinates at the given arc-length positions."""
        pts = shapely.line_interpolate_point(self.line, np.atleast_1d(arc))
        return shapely.get_coordinates(pts)

    def is_treated(self, xy) -> np.ndarray:
        xy = self._xy(xy)
        if self.treated_region is not None:
            region, tol = self.treated_region, self.tol
            return chunked(lambda a: shapely.contains_xy(region, a[:, 0], a[:, 1])
                           | (shapely.distance(region, shapely.points(a)) <= tol), len(xy), xy)
        arc = self.project(xy)
        eps = min(1e-3, self.length * 1e-6)
        a = self.interpolate(np.clip(arc - eps, 0, self.length))
        b = self.interpolate(np.clip(arc + eps, 0, self.length))
        foot = self.interpolate(arc)
        t = b - a
        rel = xy - foot
        left = t[:, 0] * rel[:, 1] - t[:, 1] * rel[:, 0] > 0
        return left if self.treated_side == "left" else ~left

    def is_masked(self, arc) -> np.ndarray:
        arc = np.atleast_1d(np.asarray(arc, dtype=float))
        out = np.zeros(arc.shape, dtype=bool)
        for a, b in self.mask_arcs:
            out |= (arc >= a) & (arc <= b)
        if self.mask_geoms:
            pts = shapely.line_interpolate_point(self.line, arc)
            for g in self.mask_geoms:
                out |= shapely.distance(g, pts) <= self.tol
        return out

    def grid(self, n: int | None = None, spacing: float | None = None, start: float = 0.0, end: float | None = None) -> np.ndarray:
        """Equally spaced arc-length positions along the boundary."""
        end = self.length if end is None else end
        if spacing is not None:
            return np.arange(start, end + 1e-9, spacing)
        return np.linspace(start, end, 25 if n is None else n)
