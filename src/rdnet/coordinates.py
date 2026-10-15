"""Network-native coordinates relative to boundary evaluation points.

For observation ``i`` and evaluation point ``b`` on the boundary:

``s_i``
    Signed network distance from ``i`` to the boundary: the metric distance to
    the nearest boundary crossing node ``c_i`` (plus off-network snapping
    depth), negative on the control side. It does not depend on ``b``.
``o_i(b)``
    Signed lateral offset along the boundary: the metric distance from ``b`` to
    ``c_i``, signed by whether ``c_i`` lies before or after ``b`` on the line.
``d_i(b)``
    Network distance from ``b`` to ``i``, used for kernel weights.

``(s_i, o_i(b))`` play the role of ``(X1 - b1, X2 - b2)`` in rd2d's location-based
local polynomial.
"""

from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd

from .boundary import Boundary
from .network import NetworkGraph, Snap, project_lonlat


def _points_array(points) -> np.ndarray:
    if hasattr(points, "geometry") and hasattr(points, "crs"):  # GeoDataFrame / GeoSeries
        geom = points.geometry if hasattr(points, "columns") else points
        return np.column_stack([geom.x.to_numpy(), geom.y.to_numpy()])
    if isinstance(points, pd.DataFrame):
        return points[["x", "y"]].to_numpy(dtype=float)
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("points must be an (n, 2) array, a DataFrame with x/y columns, or a GeoSeries")
    return arr


class NetworkCoordinates:
    """Observations expressed in network coordinates around a boundary.

    Build with :func:`network_coordinates`. Distances from each evaluation point
    are computed on demand and cached.
    """

    def __init__(
        self,
        network: NetworkGraph,
        boundary: Boundary,
        xy: np.ndarray,
        treated: np.ndarray,
        snap: Snap,
        depth_cost: np.ndarray,
        s: np.ndarray,
        crossing_node: np.ndarray,
        crossing_arc: np.ndarray,
        crossing: Snap,
        node_treated: np.ndarray,
        node_on_boundary: np.ndarray,
        node_dist: np.ndarray,
        node_source: np.ndarray,
        eval_arc: np.ndarray,
        eval_snap: Snap,
        eval_masked: np.ndarray,
    ):
        self.network = network
        self.boundary = boundary
        self.xy = xy
        self.treated = treated
        self.snap = snap
        self.depth_cost = depth_cost
        self.s = s
        self.crossing_node = crossing_node
        #: Arc position where each observation reaches the boundary.
        self.crossing_arc = crossing_arc
        #: Network location of that crossing (a boundary node, or the point itself
        #: when it sits on a street segment running along the boundary).
        self.crossing = crossing
        self.node_treated = node_treated
        self.node_on_boundary = node_on_boundary
        self.node_dist = node_dist
        self.node_source = node_source
        self.eval_arc = eval_arc
        self.eval_snap = eval_snap
        # Nearest network node to each evaluation point (for display only).
        self.eval_node = np.where(eval_snap.cost_u <= eval_snap.cost_v, eval_snap.u, eval_snap.v)
        cu, cv = eval_snap.cost_u, eval_snap.cost_v
        tot = cu + cv
        frac = np.divide(cu, tot, out=np.zeros_like(cu), where=tot > 0)
        xy = network.coords[eval_snap.u] * (1 - frac)[:, None] + network.coords[eval_snap.v] * frac[:, None]
        #: Arc position of each evaluation point after snapping onto the network.
        self.eval_node_arc = boundary.project(xy)
        self.eval_xy = xy
        self.eval_masked = eval_masked
        self._rows: dict[int, np.ndarray] = {}

    def __repr__(self) -> str:
        return (
            f"NetworkCoordinates(n={self.n_obs}, treated={int(self.treated.sum())}, "
            f"eval_points={self.n_eval}, masked={int(self.eval_masked.sum())})"
        )

    @property
    def n_obs(self) -> int:
        return len(self.s)

    @property
    def n_eval(self) -> int:
        return len(self.eval_arc)

    def eval_table(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "arc": self.eval_arc,
                "node": [self.network.nodes[k] for k in self.eval_node],
                "node_arc": self.eval_node_arc,
                "x": self.eval_xy[:, 0],
                "y": self.eval_xy[:, 1],
                "masked": self.eval_masked,
            }
        )

    def node_row(self, j: int) -> np.ndarray:
        """Metric distance from evaluation point ``j`` to every network node."""
        if j not in self._rows:
            e = self.eval_snap
            D = self.network.distances_from([e.u[j], e.v[j]])
            self._rows[j] = np.minimum(e.cost_u[j] + D[0], e.cost_v[j] + D[1])
        return self._rows[j]

    def at(self, j: int) -> tuple[np.ndarray, np.ndarray]:
        """Kernel distance ``d_i(b_j)`` and lateral offset ``o_i(b_j)`` for every observation."""
        D = self.node_row(j)
        d = self.snap.to_nodes(D[None, :], self.depth_cost)[0]
        sign = np.sign(self.crossing_arc - self.eval_node_arc[j])
        c = self.crossing
        reach = np.minimum(c.cost_u + D[c.u], c.cost_v + D[c.v])
        o = np.where(self.crossing_node >= 0, sign * reach, np.nan)
        return d, o

    def frame(self, j: int | None = None) -> pd.DataFrame:
        """Per-observation coordinates (and, given ``j``, ``d`` and ``o`` at that point)."""
        out = pd.DataFrame(
            {
                "x": self.xy[:, 0],
                "y": self.xy[:, 1],
                "treated": self.treated,
                "s": self.s,
                "crossing_arc": self.crossing_arc,
                "depth": self.snap.depth,
            }
        )
        if j is not None:
            out["d"], out["o"] = self.at(j)
        return out

    def pairwise(self, idx: np.ndarray, cols: np.ndarray | None = None, cutoff: float = np.inf) -> np.ndarray:
        """Network distance from observations ``idx`` (rows) to ``cols`` (default: ``idx``).

        Distances beyond ``cutoff`` may be reported as ``inf``.
        """
        idx = np.asarray(idx)
        cols = idx if cols is None else np.asarray(cols)
        sn = self.snap
        row_nodes, row_inv = np.unique(np.r_[sn.u[idx], sn.v[idx]], return_inverse=True)
        Dn = self.network.distances_from(row_nodes, limit=cutoff)
        iu, iv = row_inv[: len(idx)], row_inv[len(idx):]
        cu_r, cv_r = sn.cost_u[idx], sn.cost_v[idx]
        cu_c, cv_c = sn.cost_u[cols], sn.cost_v[cols]
        best = np.full((len(idx), len(cols)), np.inf)
        for ri, cr in ((iu, cu_r), (iv, cv_r)):
            Drow = Dn[ri]
            for cn, cc in ((sn.u[cols], cu_c), (sn.v[cols], cv_c)):
                best = np.minimum(best, cr[:, None] + Drow[:, cn] + cc[None, :])
        dr, dc = self.depth_cost[idx], self.depth_cost[cols]
        out = best + dr[:, None] + dc[None, :]
        same_edge = (sn.u[idx][:, None] == sn.u[cols][None, :]) & (sn.v[idx][:, None] == sn.v[cols][None, :])
        along = np.abs(cu_r[:, None] - cu_c[None, :]) + dr[:, None] + dc[None, :]
        out = np.where(same_edge, np.minimum(out, along), out)
        out[idx[:, None] == cols[None, :]] = 0.0
        return out

    def hac_quadratic(self, idx: np.ndarray, Psi: np.ndarray, cutoff: float, distance: str = "network",
                      max_cells: int = 20_000_000) -> np.ndarray:
        """Conley (Bartlett) HAC quadratic form ``Psi' K Psi`` for observations ``idx``.

        ``Psi`` has one row per entry of ``idx`` and one column per estimand; the
        result is the joint covariance matrix of the column sums. The kernel
        matrix is built block by block so memory stays bounded.
        """
        idx = np.asarray(idx)
        Psi = np.asarray(Psi, dtype=float).reshape(len(idx), -1)
        if distance not in ("network", "euclidean"):
            raise ValueError("distance must be 'network' or 'euclidean'")
        out = np.zeros((Psi.shape[1], Psi.shape[1]))
        block = max(1, max_cells // max(len(idx), 1))
        for start in range(0, len(idx), block):
            rows = idx[start:start + block]
            if distance == "network":
                Dm = self.pairwise(rows, idx, cutoff=cutoff)
            else:
                diff = self.xy[rows][:, None, :] - self.xy[idx][None, :, :]
                Dm = np.sqrt((diff**2).sum(-1))
            K = np.clip(1.0 - Dm / cutoff, 0.0, None)
            out += Psi[start:start + block].T @ (K @ Psi)
        return (out + out.T) / 2


def network_coordinates(
    network: NetworkGraph,
    boundary: Boundary,
    points,
    *,
    treated: Sequence[bool] | None = None,
    eval_points=None,
    n_eval: int = 25,
    snap: str = "edge",
    lonlat: bool = False,
    split_edges: bool = True,
) -> NetworkCoordinates:
    """Translate point locations into network coordinates around a boundary.

    Parameters
    ----------
    network
        The network, carrying its metric.
    boundary
        The boundary line and treated side.
    points
        Observation locations: ``(n, 2)`` array, DataFrame with ``x``/``y``, or a
        GeoSeries/GeoDataFrame of points, in the network's projected CRS (or
        lon/lat with ``lonlat=True``).
    treated
        Treatment indicator. Defaults to the side of the boundary each point is on.
    eval_points
        Arc-length positions along the boundary, or ``(k, 2)`` coordinates to
        project. Defaults to ``n_eval`` equally spaced points spanning the
        boundary nodes. Each evaluation point is placed on the nearest street
        segment running along the boundary, or at the nearest boundary node
        when no such segment exists.
    snap
        ``"edge"`` (attach to the nearest street segment) or ``"node"``.
    split_edges
        Insert nodes where edges cross the boundary (recommended).
    """
    xy = _points_array(points)
    if lonlat:
        if network.crs is None:
            raise ValueError("network has no CRS to project lon/lat into")
        xy = project_lonlat(xy[:, 0], xy[:, 1], network.crs)

    if split_edges:
        network = network.with_boundary_crossings(boundary)

    node_on = boundary.distance(network.coords) <= boundary.tol
    node_T = boundary.is_treated(network.coords) & ~node_on
    boundary_nodes = np.flatnonzero(node_on)
    if boundary_nodes.size == 0:
        raise ValueError("no network nodes lie on the boundary; check the CRS or increase Boundary(tol=...)")

    T = boundary.is_treated(xy) if treated is None else np.asarray(treated, dtype=bool)
    if T.shape != (len(xy),):
        raise ValueError("treated must have one entry per point")

    # Snap each observation to the network on its own side of the boundary.
    n = len(xy)
    u = np.empty(n, np.int64)
    v = np.empty(n, np.int64)
    cu = np.empty(n)
    cv = np.empty(n)
    depth = np.empty(n)
    for side_mask, allowed in ((T, node_T | node_on), (~T, ~node_T | node_on)):
        if side_mask.any():
            sn = network.snap(xy[side_mask], allowed=allowed, mode=snap)
            u[side_mask], v[side_mask] = sn.u, sn.v
            cu[side_mask], cv[side_mask] = sn.cost_u, sn.cost_v
            depth[side_mask] = sn.depth
    snapped = Snap(u, v, cu, cv, depth)
    depth_cost = depth * network.snap_scale

    node_dist, node_src = network.nearest_source(boundary_nodes)
    via_u = cu + node_dist[u]
    via_v = cv + node_dist[v]
    dist = depth_cost + np.minimum(via_u, via_v)
    crossing = np.where(via_u <= via_v, node_src[u], node_src[v])
    # A point snapped onto a segment that runs along the boundary is already there.
    on_edge = node_on[u] & node_on[v] & (u != v)
    dist = np.where(on_edge, depth_cost, dist)
    crossing = np.where(on_edge, np.where(cu <= cv, u, v), crossing)
    unreachable = ~np.isfinite(dist) | (crossing < 0)
    if unreachable.any():
        warnings.warn(f"{int(unreachable.sum())} observations cannot reach the boundary and are dropped", stacklevel=2)
    s = np.where(unreachable, np.nan, np.where(T, dist, -dist))
    crossing = np.where(unreachable, -1, crossing)
    safe = np.maximum(crossing, 0)
    cross_loc = Snap(
        np.where(on_edge, u, safe),
        np.where(on_edge, v, safe),
        np.where(on_edge, cu, 0.0),
        np.where(on_edge, cv, 0.0),
        np.zeros(n),
    )
    cross_arc = np.where(on_edge, boundary.project(xy), boundary.project(network.coords[safe]))
    cross_arc = np.where(unreachable, np.nan, cross_arc)

    bnode_arc = boundary.project(network.coords[boundary_nodes])
    if eval_points is None:
        eval_arc = np.linspace(bnode_arc.min(), bnode_arc.max(), n_eval)
    else:
        ep = np.asarray(eval_points, dtype=float)
        eval_arc = boundary.project(ep) if ep.ndim == 2 else ep.ravel()
    eval_snap = _snap_eval_points(network, boundary, eval_arc, node_on, boundary_nodes, bnode_arc)
    masked = boundary.is_masked(eval_arc)

    return NetworkCoordinates(
        network,
        boundary,
        xy,
        T,
        snapped,
        depth_cost,
        s,
        crossing,
        cross_arc,
        cross_loc,
        node_T,
        node_on,
        node_dist,
        node_src,
        eval_arc,
        eval_snap,
        masked,
    )


def _snap_eval_points(network, boundary, eval_arc, node_on, boundary_nodes, bnode_arc) -> Snap:
    """Place evaluation points on boundary street segments (or boundary nodes)."""
    k = len(eval_arc)
    nearest = boundary_nodes[np.abs(bnode_arc[None, :] - eval_arc[:, None]).argmin(axis=1)]
    out = Snap(nearest, nearest.copy(), np.zeros(k), np.zeros(k), np.zeros(k))
    if not network.metric.supports_edge_snapping:
        return out
    E = network.edge_frame()
    u, v = E["u"].to_numpy(), E["v"].to_numpy()
    along = node_on[u] & node_on[v] & np.isfinite(E["cost"].to_numpy())
    if not along.any():
        return out
    E = E[along]
    au = boundary.project(network.coords[E["u"].to_numpy()])
    av = boundary.project(network.coords[E["v"].to_numpy()])
    lo, hi = np.minimum(au, av), np.maximum(au, av)
    for j, a in enumerate(eval_arc):
        inside = np.flatnonzero((lo <= a) & (a <= hi) & (hi > lo))
        if inside.size == 0:
            continue
        r = inside[0]
        frac = (a - au[r]) / (av[r] - au[r])
        cost = float(E["cost"].iloc[r])
        out.u[j], out.v[j] = E["u"].iloc[r], E["v"].iloc[r]
        out.cost_u[j], out.cost_v[j] = frac * cost, (1 - frac) * cost
    return out
