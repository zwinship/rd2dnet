"""Network metrics: how distance (dissimilarity) is measured on a network.

A metric turns a graph into a shortest-path engine. Three families are provided:

``EdgeCost``
    Distance is the sum of per-edge costs along the cheapest path. The cost can
    be an edge attribute (``"length"``, ``"travel_time"``) or any callable
    ``f(u, v, data) -> float``.

``TransitionCost``
    Edge costs plus a penalty for *changing* edge class (e.g. highway -> local
    street). Edges whose class is listed in ``connectors`` (e.g. ramps) are
    transparent: the penalty is charged between the last non-connector class
    and the next one, so highway -> ramp -> highway stays cheap while
    highway -> ramp -> residential is expensive.

``Precomputed``
    A user-supplied node-by-node dissimilarity matrix.

All costs must be non-negative; ``inf`` marks an untraversable edge.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass
from typing import Any, Callable, Hashable, Mapping, Sequence

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra

_EPS = 1e-12

EdgeCostFn = Callable[[Hashable, Hashable, Mapping[str, Any]], float]


def _scalar(value: Any) -> Any:
    """OSM attributes are sometimes lists (merged ways); use the first entry."""
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _edge_list(graph, node_index):
    """``(ui, vi, u_labels, v_labels, datas)`` for every stored edge (self-loops dropped)."""
    E = [(u, v, d) for u, v, d in graph.edges(data=True) if u != v]
    ui = np.fromiter((node_index[e[0]] for e in E), dtype=np.int64, count=len(E))
    vi = np.fromiter((node_index[e[1]] for e in E), dtype=np.int64, count=len(E))
    return ui, vi, [e[0] for e in E], [e[1] for e in E], [e[2] for e in E]


def _orient(graph, direction, ui, vi, *cols):
    """Expand stored edges to traversable directed arcs: ``(tails, heads, *cols)``."""
    both = (not graph.is_directed()) or direction == "both"
    if both:
        return (np.r_[ui, vi], np.r_[vi, ui], *(np.r_[c, c] for c in cols))
    if direction == "out":
        return (ui, vi, *cols)
    return (vi, ui, *cols)


def _check_costs(c: np.ndarray, u_labels, v_labels) -> np.ndarray:
    bad = np.flatnonzero(np.isnan(c) | (c < 0))
    if bad.size:
        k = bad[0]
        raise ValueError(f"metric produced invalid cost {c[k]!r} for edge ({u_labels[k]!r}, {v_labels[k]!r}); "
                         "costs must be >= 0")
    return c


def _min_csr(rows, cols, vals, n) -> sparse.csr_matrix:
    """Sparse matrix keeping the minimum value among duplicate entries."""
    rows = np.asarray(rows, dtype=np.int64)
    cols = np.asarray(cols, dtype=np.int64)
    vals = np.asarray(vals, dtype=float)
    if rows.size == 0:
        return sparse.csr_matrix((n, n))
    key = rows * n + cols
    order = np.lexsort((vals, key))
    key = key[order]
    keep = np.r_[True, key[1:] != key[:-1]]
    sel = order[keep]
    # csgraph ignores stored zeros on some code paths; keep zero-cost edges alive.
    return sparse.csr_matrix(
        (np.maximum(vals[sel], _EPS), (rows[sel], cols[sel])), shape=(n, n)
    )


def _check_cost(c: float, where: str) -> float:
    c = float(c)
    if np.isnan(c) or c < 0:
        raise ValueError(f"metric produced invalid cost {c!r} for {where}; costs must be >= 0")
    return c


@dataclass
class Engine:
    """Shortest-path engine over a (possibly state-expanded) graph.

    ``start_state[v]`` is the state for standing at node ``v``; ``state_node[s]``
    is the node a state sits at. For plain edge-cost metrics states are nodes.
    """

    csr: sparse.csr_matrix
    n_nodes: int
    start_state: np.ndarray
    state_node: np.ndarray
    identity: bool

    def distances(self, sources, limit: float = np.inf) -> np.ndarray:
        """Distances from each source node to every node, shape (k, n_nodes)."""
        sources = np.atleast_1d(np.asarray(sources, dtype=np.int64))
        D = np.atleast_2d(
            dijkstra(self.csr, directed=True, indices=self.start_state[sources], limit=limit)
        )
        if self.identity:
            return D
        out = np.full((len(sources), self.n_nodes), np.inf)
        np.minimum.at(out.T, self.state_node, D.T)
        return out

    def pairs_within(self, sources, limit: float):
        """``(i, node, dist)`` for every node closer than ``limit`` to ``sources[i]``.

        Paths shorter than ``limit`` only visit nodes within ``limit`` of their
        source, so one bounded multi-source search finds the region the
        sources can reach and the per-source searches run on that induced
        subgraph. Cost scales with the reached region, not the whole network;
        pass spatially clustered sources.
        """
        sources = np.atleast_1d(np.asarray(sources, dtype=np.int64))
        st = self.start_state[sources]
        reach = dijkstra(self.csr, directed=True, indices=st, limit=limit, min_only=True)
        sub = np.flatnonzero(reach < limit)
        pos = np.full(self.csr.shape[0], -1, dtype=np.int64)
        pos[sub] = np.arange(len(sub))
        S = self.csr[sub][:, sub]
        D = np.atleast_2d(dijkstra(S, directed=True, indices=pos[st], limit=limit))
        i, c = np.nonzero(D < limit)
        d = D[i, c]
        node = self.state_node[sub[c]]
        if not self.identity:  # several states per node: keep the shortest
            order = np.lexsort((d, node, i))
            i, node, d = i[order], node[order], d[order]
            first = np.r_[True, (i[1:] != i[:-1]) | (node[1:] != node[:-1])]
            i, node, d = i[first], node[first], d[first]
        return i, node, d

    def nearest(self, sources):
        """Distance from the nearest source to every node, and which source it was."""
        sources = np.atleast_1d(np.asarray(sources, dtype=np.int64))
        dist, _, src = dijkstra(
            self.csr,
            directed=True,
            indices=self.start_state[sources],
            min_only=True,
            return_predecessors=True,
        )
        src_node = np.where(src >= 0, self.state_node[np.maximum(src, 0)], -1)
        if self.identity:
            return dist, src_node
        order = np.lexsort((dist, self.state_node))
        sn = self.state_node[order]
        first = np.r_[True, sn[1:] != sn[:-1]]
        out_d = np.full(self.n_nodes, np.inf)
        out_s = np.full(self.n_nodes, -1, dtype=np.int64)
        out_d[sn[first]] = dist[order[first]]
        out_s[sn[first]] = src_node[order[first]]
        return out_d, out_s


@dataclass
class DenseEngine:
    """Engine backed by a full node-by-node distance matrix."""

    D: np.ndarray

    @property
    def n_nodes(self) -> int:
        return self.D.shape[0]

    def distances(self, sources, limit: float = np.inf) -> np.ndarray:
        out = self.D[np.atleast_1d(sources)].astype(float, copy=True)
        out[out > limit] = np.inf
        return out

    def pairs_within(self, sources, limit: float):
        D = self.D[np.atleast_1d(sources)]
        i, node = np.nonzero(D < limit)
        return i, node, D[i, node]

    def nearest(self, sources):
        sources = np.atleast_1d(np.asarray(sources, dtype=np.int64))
        sub = self.D[sources]
        arg = np.argmin(sub, axis=0)
        dist = sub[arg, np.arange(sub.shape[1])]
        return dist, np.where(np.isfinite(dist), sources[arg], -1)


class NetworkMetric:
    """Base class for network metrics."""

    #: Whether boundary-crossing edges may be split (costs scale with the split).
    supports_edge_splitting: bool = True
    #: Whether observations may be snapped onto edges (needs per-edge costs).
    supports_edge_snapping: bool = True

    def compile(self, graph, node_index: Mapping[Hashable, int], direction: str):
        raise NotImplementedError

    def edge_cost(self, u, v, data) -> float:
        """Cost of traversing a single edge, used when snapping onto edges."""
        raise NotImplementedError

    def edge_costs(self, u_labels, v_labels, datas) -> np.ndarray:
        """Vectorised :meth:`edge_cost` over many edges."""
        return np.array([self.edge_cost(u, v, d) for u, v, d in zip(u_labels, v_labels, datas)], dtype=float)

    @property
    def default_snap_scale(self) -> float:
        """Metric cost per unit of straight-line snapping distance."""
        return 0.0


class EdgeCost(NetworkMetric):
    """Path cost = sum of edge costs.

    Parameters
    ----------
    cost
        Edge attribute name, or a callable ``f(u, v, data)``. A callable that
        also accepts a ``graph`` keyword receives the networkx graph, which
        gives access to node attributes (see :func:`node_contrast`).
    snap_scale
        Cost per metre of off-network snapping distance. Defaults to 1 when the
        cost is ``"length"`` and 0 otherwise.
    """

    def __init__(self, cost: str | EdgeCostFn = "length", *, snap_scale: float | None = None, name: str | None = None):
        self.cost = cost
        self._snap_scale = snap_scale
        self.name = name or (cost if isinstance(cost, str) else getattr(cost, "__name__", "custom"))
        self._graph = None
        self._wants_graph = callable(cost) and "graph" in inspect.signature(cost).parameters

    def __repr__(self) -> str:
        return f"EdgeCost({self.name!r})"

    @property
    def default_snap_scale(self) -> float:
        if self._snap_scale is not None:
            return float(self._snap_scale)
        return 1.0 if self.cost == "length" else 0.0

    def edge_cost(self, u, v, data) -> float:
        if isinstance(self.cost, str):
            if self.cost not in data:
                raise KeyError(f"edge ({u!r}, {v!r}) has no attribute {self.cost!r}")
            return float(_scalar(data[self.cost]))
        if self._wants_graph:
            return float(self.cost(u, v, data, graph=self._graph))
        return float(self.cost(u, v, data))

    def edge_costs(self, u_labels, v_labels, datas) -> np.ndarray:
        if isinstance(self.cost, str):
            attr = self.cost
            try:
                raw = [d[attr] for d in datas]
            except KeyError:
                k = next(i for i, d in enumerate(datas) if attr not in d)
                raise KeyError(f"edge ({u_labels[k]!r}, {v_labels[k]!r}) has no attribute {attr!r}") from None
            try:
                return np.asarray(raw, dtype=float)
            except (TypeError, ValueError):  # OSM list-valued attributes
                return np.array([float(_scalar(r)) for r in raw])
        return super().edge_costs(u_labels, v_labels, datas)

    def compile(self, graph, node_index, direction):
        self._graph = graph
        ui, vi, ul, vl, datas = _edge_list(graph, node_index)
        c = _check_costs(self.edge_costs(ul, vl, datas), ul, vl)
        fin = np.isfinite(c)
        rows, cols, vals = _orient(graph, direction, ui[fin], vi[fin], c[fin])
        n = len(node_index)
        idx = np.arange(n)
        return Engine(_min_csr(rows, cols, vals, n), n, idx, idx, identity=True)


class TransitionCost(NetworkMetric):
    """Edge costs plus penalties for switching edge class.

    Parameters
    ----------
    edge_cost
        Base per-edge cost (attribute name or callable), as in :class:`EdgeCost`.
    class_attr
        Edge attribute holding the edge class (OSM's ``"highway"`` by default).
    penalty
        ``{(from_class, to_class): cost}`` or a callable ``f(from_class, to_class)``.
        Missing pairs get ``default_penalty`` when the classes differ and 0 otherwise.
    connectors
        Classes that do not reset the remembered class (ramps, links, crosswalks).
    """

    def __init__(
        self,
        edge_cost: str | EdgeCostFn = "length",
        *,
        class_attr: str = "highway",
        penalty: Mapping[tuple[Any, Any], float] | Callable[[Any, Any], float] | None = None,
        default_penalty: float = 0.0,
        connectors: Sequence[Any] = (),
        default_class: Any = "unclassified",
        snap_scale: float | None = None,
    ):
        self._base = EdgeCost(edge_cost, snap_scale=snap_scale)
        self.class_attr = class_attr
        self.penalty = penalty or {}
        self.default_penalty = float(default_penalty)
        self.connectors = set(connectors)
        self.default_class = default_class

    def __repr__(self) -> str:
        return f"TransitionCost({self._base.name!r}, class_attr={self.class_attr!r})"

    @property
    def default_snap_scale(self) -> float:
        return self._base.default_snap_scale

    def edge_cost(self, u, v, data) -> float:
        return self._base.edge_cost(u, v, data)

    def transition_penalty(self, a, b) -> float:
        if callable(self.penalty):
            c = self.penalty(a, b)
        elif (a, b) in self.penalty:
            c = self.penalty[(a, b)]
        else:
            c = 0.0 if a == b else self.default_penalty
        return _check_cost(c, f"transition {a!r} -> {b!r}")

    def edge_costs(self, u_labels, v_labels, datas) -> np.ndarray:
        return self._base.edge_costs(u_labels, v_labels, datas)

    def compile(self, graph, node_index, direction):
        self._base._graph = graph
        n = len(node_index)
        ui, vi, ul, vl, datas = _edge_list(graph, node_index)
        c = _check_costs(self.edge_costs(ul, vl, datas), ul, vl)
        cls_raw = np.empty(len(datas), dtype=object)
        cls_raw[:] = [_scalar(d.get(self.class_attr, self.default_class)) for d in datas]
        fin = np.isfinite(c)
        tail, head, cost, ecls = _orient(graph, direction, ui[fin], vi[fin], c[fin], cls_raw[fin])
        m = len(tail)
        is_conn = np.fromiter((x in self.connectors for x in ecls), dtype=bool, count=m)

        classes = sorted({x for x, cn in zip(ecls, is_conn) if not cn}, key=str)
        class_id = {x: i + 1 for i, x in enumerate(classes)}  # memory 0 = nothing yet
        n_mem = len(classes) + 1
        pen = np.zeros((n_mem, n_mem))
        for a in classes:
            for b in classes:
                pen[class_id[a], class_id[b]] = self.transition_penalty(a, b)
        cid = np.fromiter((class_id.get(x, 0) for x in ecls), dtype=np.int64, count=m)

        # States: one "standing at node" state per node, one per non-connector
        # arc, and one per (connector arc, remembered class).
        nc = np.flatnonzero(~is_conn)
        cn = np.flatnonzero(is_conn)
        nc_state = np.full(m, -1, dtype=np.int64)
        nc_state[nc] = n + np.arange(len(nc))
        conn_base = np.full(m, -1, dtype=np.int64)
        conn_base[cn] = n + len(nc) + np.arange(len(cn)) * n_mem
        n_states = n + len(nc) + len(cn) * n_mem
        state_node = np.empty(n_states, dtype=np.int64)
        state_node[:n] = np.arange(n)
        state_node[n:n + len(nc)] = head[nc]
        state_node[n + len(nc):] = np.repeat(head[cn], n_mem)

        # Every state (at node, with remembered class) leaves along each out-arc of its node.
        from_state = np.r_[np.arange(n), nc_state[nc], (conn_base[cn][:, None] + np.arange(n_mem)).ravel()]
        from_node = state_node[from_state]
        from_mem = np.r_[np.zeros(n, dtype=np.int64), cid[nc], np.tile(np.arange(n_mem), len(cn))]
        order = np.argsort(tail, kind="stable")
        indptr = np.r_[0, np.cumsum(np.bincount(tail, minlength=n))]
        deg = indptr[from_node + 1] - indptr[from_node]
        rows = np.repeat(from_state, deg)
        mem = np.repeat(from_mem, deg)
        start = np.repeat(indptr[from_node], deg)
        k = order[start + (np.arange(deg.sum()) - np.repeat(np.cumsum(deg) - deg, deg))]
        conn_k = is_conn[k]
        cols = np.where(conn_k, conn_base[k] + mem, nc_state[k])
        vals = cost[k] + np.where(conn_k, 0.0, pen[mem, cid[k]])
        return Engine(_min_csr(rows, cols, vals, n_states), n, np.arange(n), state_node, identity=False)


class Precomputed(NetworkMetric):
    """User-supplied dissimilarity matrix between network nodes.

    Parameters
    ----------
    matrix
        Square array; ``matrix[i, j]`` is the distance from node ``nodes[i]`` to ``nodes[j]``.
    nodes
        Node labels for the rows/columns. Defaults to the graph's node order.
    """

    supports_edge_splitting = False
    supports_edge_snapping = False

    def __init__(self, matrix, nodes: Sequence[Hashable] | None = None):
        self.matrix = np.asarray(matrix, dtype=float)
        if self.matrix.ndim != 2 or self.matrix.shape[0] != self.matrix.shape[1]:
            raise ValueError("precomputed metric must be a square matrix")
        if np.any(self.matrix < 0) or np.any(np.isnan(self.matrix)):
            raise ValueError("precomputed distances must be non-negative and not NaN")
        self.nodes = list(nodes) if nodes is not None else None

    def __repr__(self) -> str:
        return f"Precomputed(shape={self.matrix.shape})"

    def compile(self, graph, node_index, direction):
        n = len(node_index)
        if self.matrix.shape[0] != n:
            raise ValueError(
                f"precomputed matrix is {self.matrix.shape[0]}x{self.matrix.shape[0]} "
                f"but the graph has {n} nodes"
            )
        if self.nodes is None:
            return DenseEngine(self.matrix)
        perm = np.array([self.nodes.index(node) for node in node_index])
        return DenseEngine(self.matrix[np.ix_(perm, perm)])


def as_metric(spec) -> NetworkMetric:
    """Coerce a metric specification (name, callable, matrix or metric) to a metric."""
    if isinstance(spec, NetworkMetric):
        return spec
    if isinstance(spec, str) or callable(spec):
        return EdgeCost(spec)
    if isinstance(spec, np.ndarray):
        return Precomputed(spec)
    raise TypeError(f"cannot interpret {spec!r} as a network metric")


# --- convenience constructors -------------------------------------------------

def travel_time(speed_attr: str = "speed_kph", length_attr: str = "length", default_kph: float | None = None) -> EdgeCost:
    """Seconds to traverse each edge, from a length (m) and speed (km/h) attribute."""

    def _cost(u, v, data):
        speed = _scalar(data.get(speed_attr, default_kph))
        if speed is None:
            raise KeyError(f"edge ({u!r}, {v!r}) has no {speed_attr!r} and no default_kph")
        speed = float(speed)
        return np.inf if speed <= 0 else float(_scalar(data[length_attr])) / (speed / 3.6)

    _cost.__name__ = "travel_time"
    return EdgeCost(_cost, name="travel_time")


def node_contrast(attr: str, weight: float, base: str = "length", snap_scale: float | None = None) -> EdgeCost:
    """Edge cost = ``base`` + ``weight * |attr(u) - attr(v)|``.

    Use for land-use or other node-level dissimilarity: moving between nodes of
    very different character costs more than the physical length alone.
    """

    def _cost(u, v, data, graph):
        a = float(graph.nodes[u].get(attr, 0.0))
        b = float(graph.nodes[v].get(attr, 0.0))
        return float(_scalar(data[base])) + weight * abs(a - b)

    _cost.__name__ = f"{base}+{weight}*|d{attr}|"
    return EdgeCost(_cost, snap_scale=snap_scale if snap_scale is not None else (1.0 if base == "length" else 0.0))
