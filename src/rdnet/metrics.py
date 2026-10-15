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


def _directed_edges(graph, node_index, direction):
    """Yield ``(ui, vi, u, v, data)`` for every traversable edge direction."""
    both = (not graph.is_directed()) or direction == "both"
    for u, v, data in graph.edges(data=True):
        ui, vi = node_index[u], node_index[v]
        if ui == vi:
            continue
        if both or direction == "out":
            yield ui, vi, u, v, data
        if both or direction == "in":
            yield vi, ui, v, u, data


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

    def compile(self, graph, node_index, direction):
        self._graph = graph
        rows, cols, vals = [], [], []
        for ui, vi, u, v, data in _directed_edges(graph, node_index, direction):
            c = _check_cost(self.edge_cost(u, v, data), f"edge ({u!r}, {v!r})")
            if np.isinf(c):
                continue
            rows.append(ui)
            cols.append(vi)
            vals.append(c)
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

    def compile(self, graph, node_index, direction):
        self._base._graph = graph
        n = len(node_index)
        edges = []  # (tail, head, cost, class, is_connector)
        for ui, vi, u, v, data in _directed_edges(graph, node_index, direction):
            c = _check_cost(self.edge_cost(u, v, data), f"edge ({u!r}, {v!r})")
            if np.isinf(c):
                continue
            cls = _scalar(data.get(self.class_attr, self.default_class))
            edges.append((ui, vi, c, cls, cls in self.connectors))

        classes = sorted({e[3] for e in edges if not e[4]}, key=str)
        class_id = {c: i + 1 for i, c in enumerate(classes)}  # memory 0 = nothing yet
        n_mem = len(classes) + 1
        pen = np.zeros((n_mem, n_mem))
        for a in classes:
            for b in classes:
                pen[class_id[a], class_id[b]] = self.transition_penalty(a, b)

        # States: one "standing at node" state per node, one per non-connector
        # edge, and one per (connector edge, remembered class).
        state_node = list(range(n))
        nc_state: dict[int, int] = {}
        conn_state: dict[tuple[int, int], int] = {}
        for k, e in enumerate(edges):
            if not e[4]:
                nc_state[k] = len(state_node)
                state_node.append(e[1])
        for k, e in enumerate(edges):
            if e[4]:
                for m in range(n_mem):
                    conn_state[(k, m)] = len(state_node)
                    state_node.append(e[1])

        out_edges: list[list[int]] = [[] for _ in range(n)]
        for k, e in enumerate(edges):
            out_edges[e[0]].append(k)

        rows, cols, vals = [], [], []

        def emit(from_state: int, at_node: int, mem: int) -> None:
            for k in out_edges[at_node]:
                _, _, c, cls, is_conn = edges[k]
                rows.append(from_state)
                if is_conn:
                    cols.append(conn_state[(k, mem)])
                    vals.append(c)
                else:
                    cols.append(nc_state[k])
                    vals.append(c + pen[mem, class_id[cls]])

        for v in range(n):
            emit(v, v, 0)
        for k, s in nc_state.items():
            emit(s, edges[k][1], class_id[edges[k][3]])
        for (k, m), s in conn_state.items():
            emit(s, edges[k][1], m)

        n_states = len(state_node)
        return Engine(
            _min_csr(rows, cols, vals, n_states),
            n,
            np.arange(n),
            np.asarray(state_node, dtype=np.int64),
            identity=False,
        )


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
