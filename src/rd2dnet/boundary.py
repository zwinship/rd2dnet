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
        Distance (in the network's coordinate units) at which a network node
        counts as lying on the boundary. Increase it when the boundary follows
        a street but its vertices are not exactly on the street's nodes.

    Points exactly on the line count as treated with ``treated_region`` (and
    are assigned by ``treated_side`` otherwise); pass ``treated=`` to
    :func:`network_coordinates` to assign them yourself. Everything outside
    ``treated_region`` is control, so restrict the sample to the two
    neighbouring areas first.
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
        if isinstance(line, BaseGeometry) and line.geom_type in ("Polygon", "MultiPolygon"):
            raise TypeError("line is a polygon: use Boundary.from_region(polygon) to make its outline the "
                            "boundary with the inside treated")
        if isinstance(line, BaseGeometry) and line.geom_type == "MultiLineString":
            from shapely.ops import linemerge

            line = linemerge(line)
            if line.geom_type != "LineString":
                raise ValueError("the boundary has several disconnected pieces; pass one connected line")
        if isinstance(line, BaseGeometry) and line.geom_type == "LinearRing":
            line = LineString(line.coords)
        if not isinstance(line, LineString):
            try:
                pts = np.asarray(line, dtype=float)
            except (TypeError, ValueError):  # e.g. a list of shapely Points
                pts = None
            if pts is not None and (pts.ndim != 2 or pts.shape[0] < 2 or pts.shape[1] != 2):
                raise ValueError("the boundary line needs at least two (x, y) vertices")
            line = LineString(line if pts is None else pts)
        self.line = line
        if not np.isfinite(self.line.length) or self.line.length <= 0:
            raise ValueError("the boundary line has zero length (its vertices coincide) or non-finite vertices")
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
                try:
                    a, b = (float(v) for v in m)
                except (TypeError, ValueError):
                    raise ValueError(f"each mask must be a polygon or an (arc_start, arc_end) pair, got {m!r}") from None
                self.mask_arcs.append((min(a, b), max(a, b)))
        shapely.prepare(self.line)

    @classmethod
    def from_points(cls, b, X=None, t=None, *, treated_side: str | None = None, **kwargs) -> "Boundary":
        """Boundary through ordered points ``b``, the evaluation points of rd2d.

        ``b`` is a ``(J, 2)`` array of points along the boundary, in order. The
        treated side is ``treated_side`` if given; otherwise it is inferred
        from observation locations ``X`` and treatment indicators ``t``: the
        side of the line on which most treated observations lie.
        """
        b = np.asarray(b, dtype=float).reshape(-1, 2)
        keep = np.r_[True, np.any(np.diff(b, axis=0) != 0, axis=1)]
        b = b[keep]
        if len(b) < 2:
            raise ValueError("b needs at least two distinct points to define a boundary line")
        if treated_side is None:
            if X is None or t is None:
                raise ValueError("give treated_side, or X and t to infer it")
            treated_side = _infer_side(LineString(b), X, t)
        return cls(b, treated_side=treated_side, **kwargs)

    @classmethod
    def from_region(cls, region, *, within: BaseGeometry | None = None, **kwargs) -> "Boundary":
        """Boundary formed by the outline of the treated region (a school zone,
        a city, a low-emission zone): points inside ``region`` are treated.

        ``within`` keeps only the part of the outline inside that geometry,
        for example the stretch where the zone borders its neighbour rather
        than the coast or the edge of the study area. The kept part must be a
        single connected line.
        """
        from shapely.ops import linemerge

        if hasattr(region, "geometry") and not isinstance(region, BaseGeometry):  # GeoSeries / GeoDataFrame
            geoms = region.geometry
            if len(geoms) == 0:
                raise ValueError("region is empty")
            region = geoms.iloc[0] if len(geoms) == 1 else shapely.union_all(np.asarray(geoms))
        if not isinstance(region, BaseGeometry) or region.geom_type not in ("Polygon", "MultiPolygon"):
            raise TypeError("region must be a Polygon or MultiPolygon (or a GeoSeries of them, in the network's CRS)")
        if region.geom_type == "MultiPolygon":
            region_line = shapely.boundary(region)
        else:
            region_line = region.exterior
        line = region_line if within is None else region_line.intersection(within)
        if line.geom_type == "MultiLineString":
            line = linemerge(line)
        if line.geom_type == "LinearRing":
            line = LineString(line.coords)
        if line.geom_type != "LineString" or line.is_empty:
            raise ValueError("the boundary must be a single connected line; narrow it with within=")
        return cls(line, treated_region=region, **kwargs)

    def __repr__(self) -> str:
        side = f"region" if self.treated_region is not None else self.treated_side
        return f"Boundary(length={self.length:.1f}, treated={side}, masks={len(self.mask_geoms) + len(self.mask_arcs)})"

    @property
    def closed(self) -> bool:
        """Whether the line is a closed ring (a zone outline); offsets then wrap around it."""
        return bool(self.line.is_closed)

    @property
    def length(self) -> float:
        """Length of the boundary line."""
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
        """Whether each point lies on the treated side (inside ``treated_region`` or on ``treated_side``)."""
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
        """Whether each arc-length position falls in a masked stretch of the boundary."""
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
        ring = self.closed and end is None and start == 0.0
        end = self.length if end is None else end
        if spacing is not None:
            out = np.arange(start, end + 1e-9, spacing)
            # on a closed ring the end is the start again
            return out[out < end - 1e-9] if ring else out
        k = 25 if n is None else n
        return np.linspace(start, end, k, endpoint=not ring)


def _infer_side(line: LineString, X, t) -> str:
    """The side ('left' or 'right') of ``line`` on which treated observations lie."""
    X = np.asarray(X, dtype=float).reshape(-1, 2)
    t = np.asarray(t).astype(bool).ravel()
    if X.shape[0] != t.shape[0]:
        raise ValueError("X and t must have the same number of observations")
    if t.all() or not t.any():
        raise ValueError("t must contain both treated and control observations")
    bl = Boundary(line, treated_side="left")
    left = bl.is_treated(X)
    # observations beyond the ends of the line have no clear side: use those whose
    # nearest point on the line is inside it, when there are enough of both groups
    arc = bl.project(X)
    inner = (arc > 1e-9 * bl.length) & (arc < bl.length * (1 - 1e-9))
    if t[inner].sum() >= 5 and (~t[inner]).sum() >= 5:
        left, t = left[inner], t[inner]
    agree = float(np.mean(left == t))
    side, share = ("left", agree) if agree >= 0.5 else ("right", 1 - agree)
    if share < 0.9:
        import warnings

        warnings.warn(f"only {share:.0%} of observations lie on the side of the boundary that matches t; "
                      "check that b is ordered along the boundary or pass treated_side/boundary", stacklevel=3)
    return side
