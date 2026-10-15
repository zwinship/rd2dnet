"""Network container: graph + metric, snapping, boundary splitting and loaders."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Hashable, Sequence

import networkx as nx
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import LineString, Point
from shapely.ops import substring

from .metrics import NetworkMetric, as_metric

_DIRECTIONS = ("both", "out", "in")


@dataclass
class Snap:
    """Where each point attaches to the network.

    A point reaches node ``u`` at cost ``cost_u`` and node ``v`` at ``cost_v``
    (along its edge), after an off-network ``depth`` (Euclidean units). Points
    snapped to a node have ``u == v`` and zero edge costs.
    """

    u: np.ndarray
    v: np.ndarray
    cost_u: np.ndarray
    cost_v: np.ndarray
    depth: np.ndarray

    def to_nodes(self, D: np.ndarray, depth_cost: np.ndarray) -> np.ndarray:
        """Distance from each source row of ``D`` (k, n_nodes) to each point, shape (k, n_points)."""
        D = np.atleast_2d(D)
        return depth_cost + np.minimum(self.cost_u + D[:, self.u], self.cost_v + D[:, self.v])


class NetworkGraph:
    """A spatial network together with the metric that defines distance on it.

    Parameters
    ----------
    graph
        A networkx graph whose nodes carry projected coordinates.
    metric
        How distance is measured: an edge attribute name (``"length"``), a callable
        ``f(u, v, data)``, a :class:`~rd2dnet.metrics.NetworkMetric`, or a square
        node-by-node matrix.
    direction
        ``"both"`` treats every edge as two-way (pedestrian access, the default);
        ``"out"`` follows edge direction; ``"in"`` reverses it.
    snap_scale
        Metric cost per unit of off-network snapping distance. Defaults to the
        metric's own default (1 for length, 0 otherwise).
    """

    def __init__(
        self,
        graph: nx.Graph,
        metric="length",
        *,
        direction: str = "both",
        x_attr: str = "x",
        y_attr: str = "y",
        crs=None,
        snap_scale: float | None = None,
    ):
        if direction not in _DIRECTIONS:
            raise ValueError(f"direction must be one of {_DIRECTIONS}")
        self.graph = graph
        self.metric: NetworkMetric = as_metric(metric)
        self.direction = direction
        self.x_attr, self.y_attr = x_attr, y_attr
        self.crs = crs if crs is not None else graph.graph.get("crs")
        self.nodes: list[Hashable] = list(graph.nodes)
        self.node_index = {n: i for i, n in enumerate(self.nodes)}
        try:
            self.coords = np.array(
                [[graph.nodes[n][x_attr], graph.nodes[n][y_attr]] for n in self.nodes], dtype=float
            ).reshape(-1, 2)
        except KeyError as e:
            raise KeyError(f"every node needs {x_attr!r} and {y_attr!r} coordinates") from e
        self.snap_scale = self.metric.default_snap_scale if snap_scale is None else float(snap_scale)
        self._engine = None
        self._edges: pd.DataFrame | None = None

    def __repr__(self) -> str:
        return (
            f"NetworkGraph(nodes={self.n_nodes}, edges={self.graph.number_of_edges()}, "
            f"metric={self.metric!r}, direction={self.direction!r})"
        )

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def engine(self):
        if self._engine is None:
            self._engine = self.metric.compile(self.graph, self.node_index, self.direction)
        return self._engine

    def _replace_graph(self, graph: nx.Graph) -> "NetworkGraph":
        return NetworkGraph(
            graph,
            self.metric,
            direction=self.direction,
            x_attr=self.x_attr,
            y_attr=self.y_attr,
            crs=self.crs,
            snap_scale=self.snap_scale,
        )

    # --- distances ---------------------------------------------------------

    def distances_from(self, sources, limit: float = np.inf) -> np.ndarray:
        """Metric distance from each source node index to every node, shape (k, n_nodes)."""
        return self.engine.distances(sources, limit=limit)

    def nearest_source(self, sources):
        """For every node: distance to the nearest of ``sources`` and that source's index."""
        return self.engine.nearest(sources)

    def node_distance(self, a: Hashable, b: Hashable) -> float:
        """Metric distance between two node labels."""
        return float(self.distances_from([self.node_index[a]])[0, self.node_index[b]])

    # --- edges and snapping ------------------------------------------------

    def edge_frame(self) -> pd.DataFrame:
        """One row per undirected street segment: endpoints, length, geometry, metric cost."""
        if self._edges is not None:
            return self._edges
        self.engine  # metrics that need the graph bind it during compile
        seen: set[tuple[int, int]] = set()
        rows = []
        for u, v, data in self.graph.edges(data=True):
            ui, vi = self.node_index[u], self.node_index[v]
            if ui == vi or (min(ui, vi), max(ui, vi)) in seen:
                continue
            seen.add((min(ui, vi), max(ui, vi)))
            geom = data.get("geometry")
            if geom is None:
                geom = LineString([self.coords[ui], self.coords[vi]])
            elif Point(geom.coords[0]).distance(Point(self.coords[ui])) > Point(geom.coords[-1]).distance(
                Point(self.coords[ui])
            ):
                geom = LineString(geom.coords[::-1])
            length = float(data.get("length", geom.length))
            if self.metric.supports_edge_snapping:
                cost = self.metric.edge_cost(u, v, data)
            else:
                cost = np.nan
            rows.append((ui, vi, length, cost, geom))
        self._edges = pd.DataFrame(rows, columns=["u", "v", "length", "cost", "geometry"])
        return self._edges

    def snap(self, xy, allowed: np.ndarray | None = None, mode: str = "edge") -> Snap:
        """Attach points to the network.

        ``allowed`` is a boolean mask over nodes; with ``mode="edge"`` a point may
        only snap to edges whose two endpoints are allowed.
        """
        xy = np.asarray(xy, dtype=float).reshape(-1, 2)
        allowed = np.ones(self.n_nodes, bool) if allowed is None else np.asarray(allowed, bool)
        if mode == "edge" and not self.metric.supports_edge_snapping:
            mode = "node"
        if mode == "node":
            from scipy.spatial import cKDTree

            cand = np.flatnonzero(allowed)
            if cand.size == 0:
                raise ValueError("no network nodes available for snapping")
            dist, k = cKDTree(self.coords[cand]).query(xy)
            node = cand[k]
            z = np.zeros(len(xy))
            return Snap(node, node.copy(), z, z.copy(), dist)
        if mode != "edge":
            raise ValueError("mode must be 'edge' or 'node'")

        E = self.edge_frame()
        ok = allowed[E["u"].to_numpy()] & allowed[E["v"].to_numpy()] & np.isfinite(E["cost"].to_numpy())
        E = E[ok]
        if len(E) == 0:
            raise ValueError("no network edges available for snapping")
        geoms = E["geometry"].to_numpy()
        tree = shapely.STRtree(geoms)
        pts = shapely.points(xy)
        (pi, ei), dist = tree.query_nearest(pts, return_distance=True, all_matches=False)
        order = np.argsort(pi)
        ei, dist = ei[order], dist[order]
        frac = shapely.line_locate_point(geoms[ei], pts, normalized=True)
        cost = E["cost"].to_numpy()[ei]
        return Snap(
            E["u"].to_numpy()[ei],
            E["v"].to_numpy()[ei],
            frac * cost,
            (1 - frac) * cost,
            dist,
        )

    # --- boundary preparation ----------------------------------------------

    def with_boundary_crossings(self, boundary, extensive: Sequence[str] = ("length", "travel_time")) -> "NetworkGraph":
        """Split every edge that crosses the boundary line at the crossing point.

        After splitting, every path between the two sides passes through a node
        on the boundary. Attributes listed in ``extensive`` are divided in
        proportion to the split; all others are copied.
        """
        if not self.metric.supports_edge_splitting:
            warnings.warn(
                f"{self.metric!r} cannot split edges; make sure boundary nodes exist in the graph",
                stacklevel=2,
            )
            return self
        T = boundary.is_treated(self.coords)
        on = boundary.distance(self.coords) <= boundary.tol
        G = self.graph.copy()
        multi = G.is_multigraph()
        edges = list(G.edges(keys=True, data=True)) if multi else [(u, v, None, d) for u, v, d in G.edges(data=True)]
        crossing_ids: dict[tuple[float, float], Hashable] = {}
        n_split = 0
        for u, v, key, data in edges:
            iu, iv = self.node_index[u], self.node_index[v]
            if on[iu] or on[iv] or T[iu] == T[iv]:
                continue
            geom = data.get("geometry")
            if geom is None:
                geom = LineString([self.coords[iu], self.coords[iv]])
            hit = geom.intersection(boundary.line)
            if hit.is_empty:
                continue  # crosses the side partition away from this boundary line
            parts = [p for p in shapely.get_parts(hit) if p.geom_type == "Point"] or [
                Point(shapely.get_coordinates(hit)[0])
            ]
            point = min(parts, key=lambda p: geom.project(p))
            t = geom.project(point) / geom.length if geom.length > 0 else 0.5
            key_xy = (round(point.x, 6), round(point.y, 6))
            if key_xy not in crossing_ids:
                cid = f"rdnet_crossing_{len(crossing_ids)}"
                G.add_node(cid, **{self.x_attr: point.x, self.y_attr: point.y, "rdnet_crossing": True})
                crossing_ids[key_xy] = cid
            cid = crossing_ids[key_xy]
            first, second = dict(data), dict(data)
            for a in extensive:
                if a in data:
                    first[a] = float(data[a]) * t
                    second[a] = float(data[a]) * (1 - t)
            if "geometry" in data:
                first["geometry"] = substring(geom, 0, t, normalized=True)
                second["geometry"] = substring(geom, t, 1, normalized=True)
            if multi:
                G.remove_edge(u, v, key)
            else:
                G.remove_edge(u, v)
            G.add_edge(u, cid, **first)
            G.add_edge(cid, v, **second)
            n_split += 1
        out = self._replace_graph(G)
        out.n_crossings_inserted = n_split
        return out

    # --- constructors --------------------------------------------------------

    @classmethod
    def from_osmnx(cls, graph, metric="length", *, project: bool = True, **kwargs) -> "NetworkGraph":
        """Wrap an OSMnx graph, projecting it to a metric CRS if needed."""
        if project and not _is_projected(graph):
            import osmnx as ox

            graph = ox.project_graph(graph)
        return cls(graph, metric, crs=graph.graph.get("crs"), **kwargs)

    @classmethod
    def from_place(
        cls,
        query: str | None = None,
        *,
        bbox: tuple[float, float, float, float] | None = None,
        network_type: str = "walk",
        metric="length",
        offline_fallback: bool = True,
        **kwargs,
    ) -> "NetworkGraph":
        """Download an OpenStreetMap network with OSMnx.

        If OSMnx is missing or the download fails and ``offline_fallback`` is
        true, the bundled mock Midtown network is returned with a warning, so
        scripts and tests keep working without internet access.
        """
        try:
            import osmnx as ox

            if bbox is not None:
                G = ox.graph_from_bbox(bbox, network_type=network_type)
            elif query is not None:
                G = ox.graph_from_place(query, network_type=network_type)
            else:
                raise ValueError("give a place query or a bbox")
        except ValueError:
            raise
        except Exception as e:  # network errors, missing osmnx, Overpass failures
            if not offline_fallback:
                raise
            from .datasets import load_mock_midtown

            warnings.warn(f"OSM download failed ({e!r}); using the bundled mock Midtown network", stacklevel=2)
            return load_mock_midtown().network(metric, **kwargs)
        return cls.from_osmnx(G, metric, **kwargs)

    @classmethod
    def from_adjacency(cls, adjacency, coords, *, nodes: Sequence[Hashable] | None = None, directed: bool = False, metric="weight", **kwargs) -> "NetworkGraph":
        """Build from a (sparse or dense) weighted adjacency matrix and node coordinates.

        Each edge gets ``weight`` (the matrix entry) and ``length`` (straight-line).
        """
        from scipy import sparse

        A = sparse.coo_matrix(adjacency)
        coords = np.asarray(coords, dtype=float)
        nodes = list(nodes) if nodes is not None else list(range(A.shape[0]))
        G = nx.DiGraph() if directed else nx.Graph()
        for n, (x, y) in zip(nodes, coords):
            G.add_node(n, x=x, y=y)
        for i, j, w in zip(A.row, A.col, A.data):
            if i == j:
                continue
            G.add_edge(nodes[i], nodes[j], weight=float(w), length=float(np.hypot(*(coords[i] - coords[j]))))
        return cls(G, metric, **kwargs)

    @classmethod
    def from_igraph(cls, g, metric="weight", *, x_attr: str = "x", y_attr: str = "y", **kwargs) -> "NetworkGraph":
        """Build from a python-igraph graph whose vertices carry coordinates."""
        G = g.to_networkx()
        return cls(G, metric, x_attr=x_attr, y_attr=y_attr, **kwargs)


def _is_projected(graph) -> bool:
    crs = graph.graph.get("crs")
    if crs is None:
        return True
    return "4326" not in str(crs) and "longlat" not in str(crs).lower()


def project_lonlat(lon, lat, to_crs) -> np.ndarray:
    """Project longitude/latitude to ``to_crs`` (requires pyproj)."""
    try:
        from pyproj import Transformer
    except ImportError as e:  # pragma: no cover - depends on environment
        raise ImportError("projecting lon/lat requires pyproj (pip install pyproj)") from e
    tr = Transformer.from_crs("EPSG:4326", to_crs, always_xy=True)
    x, y = tr.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.column_stack([x, y])
