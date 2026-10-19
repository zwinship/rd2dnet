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
from concurrent.futures import ThreadPoolExecutor
from typing import Sequence

import numpy as np
import pandas as pd

from .boundary import Boundary
from ._parallel import n_threads
from .network import NetworkGraph, Snap, project_lonlat


def _points_array(points, crs=None) -> np.ndarray:
    if hasattr(points, "geometry") and hasattr(points, "crs"):  # GeoDataFrame / GeoSeries
        geom = points.geometry if hasattr(points, "columns") else points
        if crs is not None and geom.crs is not None and geom.crs != crs:
            geom = geom.to_crs(crs)  # match the network's CRS
        return np.column_stack([geom.x.to_numpy(), geom.y.to_numpy()])
    if isinstance(points, pd.DataFrame):
        if {"x", "y"} <= set(points.columns):
            return points[["x", "y"]].to_numpy(dtype=float)
        if points.shape[1] == 2:
            return points.to_numpy(dtype=float)
        raise ValueError("a DataFrame of points needs x and y columns (or exactly two columns)")
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("points must be an (n, 2) array, a DataFrame with x/y columns, or a GeoSeries")
    return arr


class NetworkCoordinates:
    """Observations expressed in network coordinates around a boundary.

    Build with :func:`network_coordinates`. Distances from each evaluation point
    are computed on demand and cached.

    Attributes
    ----------
    s
        Signed network distance to the boundary (positive on the treated side,
        NaN if the boundary cannot be reached).
    treated
        Treatment indicator per observation.
    xy
        Observation coordinates.
    crossing_arc
        Arc-length position where each observation reaches the boundary.
    eval_arc, eval_node_arc, eval_xy, eval_masked
        Requested and snapped evaluation points, and whether each is masked.
    loc
        Location id of each observation (repeated coordinates share one).
    network, boundary
        The (boundary-split) network and the boundary.
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
        loc: np.ndarray | None = None,
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
        # Observations at the same location (and side) share all coordinates:
        # ``loc`` maps each observation to its location and ``_first`` holds one
        # observation per location, so per-location work is done once.
        n = len(s)
        self.loc = np.arange(n) if loc is None else np.asarray(loc)
        first = np.full(int(self.loc.max()) + 1 if n else 0, -1)
        first[self.loc[::-1]] = np.arange(n)[::-1]
        self._first = f = first
        self._lsnap = Snap(snap.u[f], snap.v[f], snap.cost_u[f], snap.cost_v[f], snap.depth[f])
        self._ldepth = depth_cost[f]
        self._lcross = Snap(crossing.u[f], crossing.v[f], crossing.cost_u[f], crossing.cost_v[f], crossing.depth[f])
        self._lcross_arc = crossing_arc[f]
        self._lcross_node = crossing_node[f]

    def __repr__(self) -> str:
        return (
            f"NetworkCoordinates(n={self.n_obs}, treated={int(self.treated.sum())}, "
            f"eval_points={self.n_eval}, masked={int(self.eval_masked.sum())})"
        )

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        return len(self.s)

    @property
    def n_eval(self) -> int:
        """Number of evaluation points."""
        return len(self.eval_arc)

    @property
    def n_locations(self) -> int:
        """Number of distinct observation locations (repeated measurements share one)."""
        return len(self._first)

    def eval_table(self) -> pd.DataFrame:
        """Evaluation points: requested and snapped arc position, nearest node, coordinates, mask."""
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

    def node_rows(self, js) -> None:
        """Compute and cache :meth:`node_row` for several evaluation points in one search."""
        js = [j for j in js if j not in self._rows]
        if not js:
            return
        e = self.eval_snap
        D = self.network.distances_from(np.r_[e.u[js], e.v[js]])
        k = len(js)
        for r, j in enumerate(js):
            self._rows[j] = np.minimum(e.cost_u[j] + D[r], e.cost_v[j] + D[k + r])

    def _along_eval_edge(self, loc: Snap, j: int) -> np.ndarray:
        """Direct distance for locations on the same segment as evaluation point ``j`` (else inf)."""
        e = self.eval_snap
        u0, v0 = e.u[j], e.v[j]
        out = np.full(len(loc.u), np.inf)
        if u0 == v0:
            return out
        same = (loc.u == u0) & (loc.v == v0)
        rev = (loc.u == v0) & (loc.v == u0)
        out[same] = np.abs(loc.cost_u[same] - e.cost_u[j])
        out[rev] = np.abs(loc.cost_u[rev] - e.cost_v[j])
        return out

    def at(self, j: int) -> tuple[np.ndarray, np.ndarray]:
        """Kernel distance ``d_i(b_j)`` and lateral offset ``o_i(b_j)`` for every observation."""
        D = self.node_row(j)
        sn = self._lsnap
        d = sn.to_nodes(D[None, :], self._ldepth)[0]
        d = np.minimum(d, self._ldepth + self._along_eval_edge(sn, j))
        diff = self._lcross_arc - self.eval_node_arc[j]
        if self.boundary.closed:  # the shorter way round the ring
            L = self.boundary.length
            diff = (diff + L / 2) % L - L / 2
        sign = np.sign(diff)
        c = self._lcross
        reach = np.minimum(c.cost_u + D[c.u], c.cost_v + D[c.v])
        reach = np.minimum(reach, self._along_eval_edge(c, j))
        o = np.where(self._lcross_node >= 0, sign * reach, np.nan)
        return d[self.loc], o[self.loc]

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
        result is the joint covariance matrix of the column sums, with
        ``K_ik = max(0, 1 - d_ik / cutoff)``. Only pairs closer than ``cutoff``
        are visited: rows of ``Psi`` are first summed by location (``K``
        depends on location only), node pairs within ``cutoff`` come from
        Dijkstra searches bounded at ``cutoff``, and the kernel is built in
        sparse row blocks. Time and memory grow with the number of close
        pairs, not with ``len(idx)**2``.
        """
        idx = np.asarray(idx)
        Psi = np.asarray(Psi, dtype=float).reshape(len(idx), -1)
        if distance not in ("network", "euclidean"):
            raise ValueError("distance must be 'network' or 'euclidean'")
        locs, inv = np.unique(self.loc[idx], return_inverse=True)
        P = np.zeros((len(locs), Psi.shape[1]))
        np.add.at(P, inv, Psi)
        threads = min(n_threads(len(locs)), 8)

        def part(block):
            r0, r1, make_K = block
            return P[r0:r1].T @ (make_K() @ P)

        # Row blocks are independent (numpy and scipy.sparse do the work in C,
        # releasing the GIL), so they run in threads; the memory budget per
        # block is divided by the thread count.
        blocks = self._kernel_blocks(locs, cutoff, distance, max(1, int(max_cells) // threads))
        if threads > 1:
            with ThreadPoolExecutor(threads) as ex:
                parts = list(ex.map(part, blocks))
        else:
            parts = [part(b) for b in blocks]
        out = np.sum(parts, axis=0) if parts else np.zeros((Psi.shape[1], Psi.shape[1]))
        return (out + out.T) / 2

    def _kernel_blocks(self, locs, cutoff, distance, max_cells):
        """Yield ``(r0, r1, make_K)``; ``make_K()`` builds sparse Bartlett kernel rows ``r0:r1``."""
        from scipy import sparse

        m = len(locs)
        rows_per = max(1, int(max_cells) // max(m, 1))
        if distance == "euclidean":
            from scipy.spatial import cKDTree

            pts = self.xy[self._first[locs]]
            tree = cKDTree(pts)
            def make_euclid(r0, r1):
                sdm = cKDTree(pts[r0:r1]).sparse_distance_matrix(tree, cutoff, output_type="ndarray")
                k = 1.0 - sdm["v"] / cutoff
                keep = k > 0
                return sparse.csr_matrix((k[keep], (sdm["i"][keep], sdm["j"][keep])), shape=(r1 - r0, m))

            for r0 in range(0, m, rows_per):
                r1 = min(m, r0 + rows_per)
                yield r0, r1, (lambda r0=r0, r1=r1: make_euclid(r0, r1))
            return

        sn = self._lsnap
        eu, ev = sn.u[locs], sn.v[locs]
        cu, cv, dep = sn.cost_u[locs], sn.cost_v[locs], self._ldepth[locs]
        nodes, ninv = np.unique(np.r_[eu, ev], return_inverse=True)
        pu, pv = ninv[:m], ninv[m:]
        nn = len(nodes)
        # Node pairs (a, b) with network distance a -> b below the cutoff, stored
        # as D + cutoff so that zero distances survive in sparse form. Sources are
        # taken in Z-order blocks so each block's bounded search stays local.
        npos = np.full(self.network.n_nodes, -1, dtype=np.int64)
        npos[nodes] = np.arange(nn)
        ordered = np.argsort(_zorder(self.network.coords[nodes]), kind="stable")
        engine = self.network.engine
        ra, rb, dab = [], [], []
        for s0 in range(0, nn, 512):
            blk = ordered[s0:s0 + 512]
            i, node, dist = engine.pairs_within(nodes[blk], cutoff)
            t = npos[node]
            keep = t >= 0
            ra.append(blk[i[keep]])
            rb.append(t[keep])
            dab.append(dist[keep] + cutoff)
        ND = sparse.csr_matrix((np.concatenate(dab), (np.concatenate(ra), np.concatenate(rb))), shape=(nn, nn))
        rng = np.arange(m)
        Au = sparse.csr_matrix((np.ones(m), (rng, pu)), shape=(m, nn))
        Av = sparse.csr_matrix((np.ones(m), (rng, pv)), shape=(m, nn))
        AuT, AvT = Au.T.tocsr(), Av.T.tocsr()
        # Each location has one u and one v endpoint, so (A_a ND A_b')[i, j] is the
        # single term D(a_i, b_j) + cutoff. Rows are blocked by their candidate count.
        pat = ND.copy()
        pat.data[:] = 1.0
        both = Au + Av
        per_row = both @ (pat @ (both.T @ np.ones(m)))
        budget = max(1, int(max_cells) // 4)
        cuts = np.searchsorted(np.cumsum(per_row), np.arange(budget, per_row.sum() + budget, budget), side="right")
        bounds = np.unique(np.r_[0, np.clip(cuts, 1, m), m])
        def make_net(r0, r1):
            K = None
            for Aa, ca, a_is_u in ((Au, cu, True), (Av, cv, False)):
                left = Aa[r0:r1] @ ND
                for ATb, cb, b_is_u in ((AuT, cu, True), (AvT, cv, False)):
                    M = (left @ ATb).tocoo()
                    i, j = M.row + r0, M.col
                    d = (M.data - cutoff) + ca[i] + cb[j]
                    if a_is_u and b_is_u:  # every same-edge pair and i == j shows up here (D(u, u) = 0)
                        same = (eu[i] == eu[j]) & (ev[i] == ev[j])
                        d = np.where(same, np.minimum(d, np.abs(cu[i] - cu[j])), d)
                    d = d + dep[i] + dep[j]
                    if a_is_u and b_is_u:
                        d[i == j] = 0.0
                    k = 1.0 - d / cutoff
                    keep = k > 0
                    Kc = sparse.csr_matrix((k[keep], (i[keep] - r0, j[keep])), shape=(r1 - r0, m))
                    K = Kc if K is None else K.maximum(Kc)  # kernel at the shortest distance
            return K

        for r0, r1 in zip(bounds[:-1], bounds[1:]):
            yield r0, r1, (lambda r0=r0, r1=r1: make_net(r0, r1))


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
    extensive: Sequence[str] = ("length", "travel_time"),
) -> NetworkCoordinates:
    """Translate point locations into network coordinates around a boundary.

    Parameters
    ----------
    network
        The network, carrying its metric.
    boundary
        The boundary line and treated side.
    points
        Observation locations: ``(n, 2)`` array, DataFrame with ``x``/``y`` (or
        exactly two columns), or a GeoSeries/GeoDataFrame of points, in the
        network's projected CRS (or lon/lat with ``lonlat=True``). A
        GeoSeries/GeoDataFrame in another CRS is reprojected to the network's.
    treated
        Treatment indicator (boolean or 0/1). Defaults to the side of the
        boundary each point is on (points on the line count as treated with a
        ``treated_region``).
    eval_points
        Arc-length positions along the boundary, or ``(k, 2)`` coordinates to
        project. Defaults to ``n_eval`` equally spaced points spanning the
        boundary nodes (around the whole ring for a closed boundary). Each evaluation point is placed on the nearest street
        segment running along the boundary, or at the nearest boundary node
        when no such segment exists.
    snap
        ``"edge"`` (attach to the nearest street segment) or ``"node"``.
    lonlat
        Points are longitude/latitude; project them into the network's CRS
        (requires pyproj).
    split_edges
        Insert nodes where edges cross the boundary (recommended).
    extensive
        Edge attributes divided in proportion when an edge is split at the
        boundary (the metric's own cost attribute always is); others are copied.

    Returns
    -------
    NetworkCoordinates
        Signed boundary distances ``s``, treatment indicators, where each
        observation reaches the boundary, and the evaluation points; kernel
        distances and lateral offsets at evaluation point ``j`` come from
        :meth:`NetworkCoordinates.at`. Repeated locations are processed once.

    Notes
    -----
    On a directed network (``direction="out"``), distances run from the
    boundary to each location, following edge directions: on a river drawn
    downstream, locations upstream of the boundary cannot be reached. Use
    ``direction="in"`` for distances from each location to the boundary, or
    ``"both"`` to ignore direction.
    """
    xy = _points_array(points, network.crs)
    if lonlat:
        if network.crs is None:
            raise ValueError("network has no CRS to project lon/lat into")
        xy = project_lonlat(xy[:, 0], xy[:, 1], network.crs)

    if not np.isfinite(xy).all():
        raise ValueError("point coordinates must be finite")
    if split_edges:
        network = network.with_boundary_crossings(boundary, extensive)

    bdist = boundary.distance(network.coords)
    node_on = bdist <= boundary.tol
    near = (bdist > boundary.tol) & (bdist <= 1e-4 * boundary.length)
    if near.sum() >= 3:
        warnings.warn(f"{int(near.sum())} network nodes lie within {1e-4 * boundary.length:.3g} of the boundary but "
                      f"outside Boundary(tol={boundary.tol:g}); if the boundary follows a street, increase tol "
                      "(in the network's coordinate units) so they count as on it", stacklevel=2)
    node_T = boundary.is_treated(network.coords) & ~node_on
    boundary_nodes = np.flatnonzero(node_on)
    if boundary_nodes.size == 0:
        raise ValueError("no network nodes lie on the boundary; check the CRS or increase Boundary(tol=...)")

    T = None if treated is None else _as_treatment(treated, len(xy))
    # Repeated locations (panels, mass points) are processed once.
    loc, first = _unique_locations(xy, T)
    xy_obs, xy = xy, xy[first]
    T = boundary.is_treated(xy) if T is None else T[first]

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

    bnode_arc = boundary.project(network.coords[boundary_nodes])
    if eval_points is None:
        if boundary.closed:  # around the ring, without repeating the start
            eval_arc = (bnode_arc.min() + np.arange(n_eval) * boundary.length / n_eval) % boundary.length
        else:
            eval_arc = np.linspace(bnode_arc.min(), bnode_arc.max(), n_eval)
    else:
        ep = np.asarray(eval_points, dtype=float)
        eval_arc = boundary.project(ep) if ep.ndim == 2 else ep.ravel()
        if eval_arc.size == 0:
            raise ValueError("eval_points is empty")
        outside = (eval_arc < bnode_arc.min() - 1e-9) | (eval_arc > bnode_arc.max() + 1e-9)
        if outside.any() and not boundary.closed:
            warnings.warn(f"{int(outside.sum())} evaluation points lie beyond the first or last place the network "
                          "meets the boundary; they are moved to the nearest boundary node", stacklevel=2)
    eval_snap = _snap_eval_points(network, boundary, eval_arc, node_on, boundary_nodes, bnode_arc)

    if network.needs_point_nodes:
        # One-way edges and class-change penalties: put every location and
        # evaluation point on a node, so that paths respect both.
        k_obs = len(u)
        pu = np.r_[u, eval_snap.u]
        pv = np.r_[v, eval_snap.v]
        pc = np.r_[cu, eval_snap.cost_u]
        tot = pc + np.r_[cv, eval_snap.cost_v]
        frac = np.divide(pc, tot, out=np.zeros_like(pc), where=tot > 0)
        on_edge_ends = node_on[pu] & node_on[pv]
        network, pnode = network._with_point_nodes(pu, pv, frac, extensive)
        on_new = boundary.distance(network.coords) <= boundary.tol
        on_new[pnode[on_edge_ends]] = True
        node_on = on_new
        node_T = boundary.is_treated(network.coords) & ~node_on
        boundary_nodes = np.flatnonzero(node_on)
        zero = np.zeros(len(pnode))
        u, v, cu, cv = pnode[:k_obs], pnode[:k_obs].copy(), zero[:k_obs], zero[:k_obs].copy()
        snapped = Snap(u, v, cu, cv, depth)
        e = pnode[k_obs:]
        eval_snap = Snap(e, e.copy(), zero[k_obs:], zero[k_obs:].copy(), np.zeros(len(e)))

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
        n_rows = int(np.isin(loc, np.flatnonzero(unreachable)).sum())
        warnings.warn(f"{n_rows} observations ({int(unreachable.sum())} locations) have no network path to the "
                      "boundary and are dropped (disconnected parts of the network, or edge directions)", stacklevel=2)
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

    masked = boundary.is_masked(eval_arc)
    spot = np.column_stack([np.minimum(eval_snap.u, eval_snap.v), np.maximum(eval_snap.u, eval_snap.v),
                            np.round(np.where(eval_snap.u <= eval_snap.v, eval_snap.cost_u, eval_snap.cost_v), 6)])
    n_spots = len(np.unique(spot[~masked], axis=0))
    if n_spots < 0.75 * (~masked).sum():  # a few repeats near the ends are harmless
        warnings.warn(f"the {int((~masked).sum())} evaluation points fall on only {n_spots} distinct places (the network "
                      "meets the boundary at few points); their estimates repeat, so pass one evaluation point "
                      "per crossing, e.g. eval_points=<crossing coordinates>", stacklevel=2)

    def per_obs(snp_or_arr):
        if isinstance(snp_or_arr, Snap):
            a = snp_or_arr
            return Snap(a.u[loc], a.v[loc], a.cost_u[loc], a.cost_v[loc], a.depth[loc])
        return snp_or_arr[loc]

    return NetworkCoordinates(
        network,
        boundary,
        xy_obs,
        per_obs(T),
        per_obs(snapped),
        per_obs(depth_cost),
        per_obs(s),
        per_obs(crossing),
        per_obs(cross_arc),
        per_obs(cross_loc),
        node_T,
        node_on,
        node_dist,
        node_src,
        eval_arc,
        eval_snap,
        masked,
        loc=loc,
    )


def _as_treatment(treated, n: int) -> np.ndarray:
    """Boolean treatment indicator from bool or 0/1 input, with checks."""
    t = np.asarray(treated).ravel()
    if t.dtype != bool:
        try:
            tf = t.astype(float)
        except (TypeError, ValueError):
            raise ValueError("treated must be boolean or 0/1") from None
        if np.isnan(tf).any() or not np.isin(tf, (0, 1)).all():
            raise ValueError("treated must be boolean or 0/1 without missing values")
        t = tf == 1
    if t.shape != (n,):
        raise ValueError(f"treated must have one entry per point ({n}), got {t.size}")
    return t


def _zorder(xy: np.ndarray) -> np.ndarray:
    """Morton (Z-order) key of 2-D points: nearby points get nearby keys."""
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    q = ((xy - lo) / np.maximum(hi - lo, 1e-12) * 65535).astype(np.uint64)
    key = np.zeros(len(xy), dtype=np.uint64)
    for b in range(16):
        bit = np.uint64(1 << b)
        key |= ((q[:, 0] & bit) << np.uint64(b)) | ((q[:, 1] & bit) << np.uint64(b + 1))
    return key


def _unique_locations(xy: np.ndarray, T: np.ndarray | None):
    """``(loc, first)``: location id of each point and one point per location (by x, y and side).

    Locations are numbered in order of first appearance.
    """
    n = len(xy)
    keys = (xy[:, 1], xy[:, 0]) if T is None else (T, xy[:, 1], xy[:, 0])
    order = np.lexsort(keys)
    a = xy[order]
    new = np.ones(n, dtype=bool)
    new[1:] = (a[1:, 0] != a[:-1, 0]) | (a[1:, 1] != a[:-1, 1])
    if T is not None:
        t = T[order]
        new[1:] |= t[1:] != t[:-1]
    loc = np.empty(n, dtype=np.int64)
    loc[order] = np.cumsum(new) - 1
    first = order[new]
    srt = np.argsort(first)
    rank = np.empty(len(first), dtype=np.int64)
    rank[srt] = np.arange(len(first))
    return rank[loc], first[srt]


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
