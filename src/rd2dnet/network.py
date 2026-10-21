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

from ._parallel import chunked
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
        metric's own default: 1 for length, and otherwise the snapped edge's
        cost per unit length (so with travel time the walk to the street is
        charged at that street's seconds per metre); 0 for a precomputed matrix.
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
        if not np.isfinite(self.coords).all():
            bad = self.nodes[int(np.flatnonzero(~np.isfinite(self.coords).all(axis=1))[0])]
            raise ValueError(f"node coordinates must be finite (node {bad!r} has {self.coords[self.node_index[bad]]})")
        if _looks_lonlat(self.coords, self.crs):
            warnings.warn("the network is in longitude/latitude, so distances are in degrees; project it to a "
                          "metric CRS first (NetworkGraph.from_osmnx and from_lines do this)", stacklevel=2)
        cost = getattr(self.metric, "cost", None)
        if isinstance(cost, str) and graph.number_of_edges():
            missing = sum(1 for *_, d in graph.edges(data=True) if cost not in d)
            if missing:
                have = sorted({k for *_, d in graph.edges(data=True) for k in d if k != "geometry"})
                raise KeyError(f"{missing} edges have no {cost!r} attribute for the metric; edge attributes "
                               f"available: {have}")
        self.snap_scale = self.metric.default_snap_scale if snap_scale is None else float(snap_scale)
        if self.snap_scale is not None and not (np.isfinite(self.snap_scale) and self.snap_scale >= 0):
            raise ValueError(f"snap_scale must be finite and non-negative, got {snap_scale!r}")
        self._engine = None
        self._edges: pd.DataFrame | None = None
        self._split: tuple | None = None  # (boundary, split network) cache

    def __repr__(self) -> str:
        return (
            f"NetworkGraph(nodes={self.n_nodes}, edges={self.graph.number_of_edges()}, "
            f"metric={self.metric!r}, direction={self.direction!r})"
        )

    @property
    def n_nodes(self) -> int:
        """Number of network nodes."""
        return len(self.nodes)

    @property
    def engine(self):
        """Shortest-path engine compiled from the metric (built on first use)."""
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
        from .metrics import _edge_list

        ui, vi, ul, vl, datas = _edge_list(self.graph, self.node_index)
        # Keep one edge per undirected node pair, in graph order: the cheapest
        # of parallel edges, which is the one shortest paths use.
        key = np.minimum(ui, vi) * self.n_nodes + np.maximum(ui, vi)
        if self.metric.supports_edge_snapping and len(np.unique(key)) < len(key):
            cost_all = np.asarray(self.metric.edge_costs(ul, vl, datas), dtype=float)
            order = np.lexsort((cost_all, key))
            first = np.sort(order[np.r_[True, key[order][1:] != key[order][:-1]]])
        else:
            _, first = np.unique(key, return_index=True)
            first = np.sort(first)
        ui, vi = ui[first], vi[first]
        datas = [datas[k] for k in first]
        ul, vl = [ul[k] for k in first], [vl[k] for k in first]
        geom = np.array([d.get("geometry") for d in datas], dtype=object)
        missing = np.array([g is None for g in geom], dtype=bool)
        if missing.any():
            seg = np.stack([self.coords[ui[missing]], self.coords[vi[missing]]], axis=1)
            geom[missing] = shapely.linestrings(seg)
        given = np.flatnonzero(~missing)
        if given.size:  # orient supplied geometries from u to v
            g = geom[given]
            pu = shapely.points(self.coords[ui[given]])
            flip = shapely.distance(shapely.get_point(g, 0), pu) > shapely.distance(shapely.get_point(g, -1), pu)
            geom[given[flip]] = shapely.reverse(g[flip])
        length = np.array([d.get("length", np.nan) for d in datas], dtype=object)
        length = np.array([np.nan if x is None else float(x if not isinstance(x, (list, tuple)) else x[0])
                           for x in length], dtype=float)
        nolen = np.isnan(length)
        length[nolen] = shapely.length(geom[nolen])
        if self.metric.supports_edge_snapping:
            cost = self.metric.edge_costs(ul, vl, datas)
        else:
            cost = np.full(len(ui), np.nan)
        self._edges = pd.DataFrame({"u": ui, "v": vi, "length": length, "cost": cost, "geometry": geom})
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

        def nearest(pts_xy):
            pts = shapely.points(pts_xy)
            (pi, ei), dist = tree.query_nearest(pts, return_distance=True, all_matches=False)
            order = np.argsort(pi)
            ei, dist = ei[order], dist[order]
            return ei, dist, shapely.line_locate_point(geoms[ei], pts, normalized=True)

        ei, dist, frac = chunked(nearest, len(xy), xy)
        cost = E["cost"].to_numpy()[ei]
        return Snap(
            E["u"].to_numpy()[ei],
            E["v"].to_numpy()[ei],
            frac * cost,
            (1 - frac) * cost,
            dist,
        )

    def snap_rates(self, snap: "Snap") -> np.ndarray:
        """Cost per unit of off-network distance for each snapped point.

        ``snap_scale`` when it is a number; otherwise the snapped edge's cost
        divided by its length (for points snapped to a node, the mean over the
        node's edges).
        """
        n = len(snap.u)
        if self.snap_scale is not None:
            return np.full(n, float(self.snap_scale))
        E = self.edge_frame()
        rate = np.asarray(E["cost"], dtype=float) / np.asarray(E["length"], dtype=float)
        good = np.isfinite(rate) & (np.asarray(E["length"], dtype=float) > 0)
        fallback = float(np.median(rate[good])) if good.any() else 0.0
        rate = np.where(good, rate, fallback)
        eu, ev = E["u"].to_numpy(), E["v"].to_numpy()
        key = np.minimum(eu, ev) * self.n_nodes + np.maximum(eu, ev)
        lookup = dict(zip(key.tolist(), rate.tolist()))
        node_sum = np.bincount(np.r_[eu, ev], weights=np.r_[rate, rate], minlength=self.n_nodes)
        node_cnt = np.bincount(np.r_[eu, ev], minlength=self.n_nodes)
        node_rate = np.where(node_cnt > 0, node_sum / np.maximum(node_cnt, 1), fallback)
        u, v = np.asarray(snap.u), np.asarray(snap.v)
        k = np.minimum(u, v) * self.n_nodes + np.maximum(u, v)
        out = np.array([lookup.get(kk, np.nan) for kk in k.tolist()], dtype=float)
        at_node = (u == v) | np.isnan(out)
        out[at_node] = node_rate[u[at_node]]
        return out

    # --- boundary preparation ----------------------------------------------

    def with_boundary_crossings(self, boundary, extensive: Sequence[str] = ("length", "travel_time")) -> "NetworkGraph":
        """Split every edge that crosses the boundary line at the crossing point.

        After splitting, every path between the two sides passes through a node
        on the boundary. Attributes listed in ``extensive`` are divided in
        proportion to the split, as is the attribute the metric reads (for
        example ``"weight"``); all others are copied.
        """
        if not self.metric.supports_edge_splitting:
            warnings.warn(
                f"{self.metric!r} cannot split edges; make sure boundary nodes exist in the graph",
                stacklevel=2,
            )
            return self
        if self._split is not None and self._split[0] is boundary:
            return self._split[1]
        # the metric's own cost attribute is always split in proportion
        extensive = set(extensive) | _cost_attributes(self.metric)
        T = boundary.is_treated(self.coords)
        on = boundary.distance(self.coords) <= boundary.tol
        multi = self.graph.is_multigraph()
        edges = (list(self.graph.edges(keys=True, data=True)) if multi
                 else [(u, v, None, d) for u, v, d in self.graph.edges(data=True)])
        iu = np.fromiter((self.node_index[e[0]] for e in edges), dtype=np.int64, count=len(edges))
        iv = np.fromiter((self.node_index[e[1]] for e in edges), dtype=np.int64, count=len(edges))
        cand = np.flatnonzero(~on[iu] & ~on[iv] & (T[iu] != T[iv]))
        if cand.size == 0:  # nothing to split: reuse this network (and its compiled engine)
            self.n_crossings_inserted = 0
            self._split = (boundary, self)
            return self
        G = self.graph.copy()
        crossing_ids: dict[tuple[float, float], Hashable] = {}
        next_cross = 0
        n_split = 0
        for k in cand:
            u, v, key, data = edges[k]
            iu_k, iv_k = iu[k], iv[k]
            geom = data.get("geometry")
            if geom is None:
                geom = LineString([self.coords[iu_k], self.coords[iv_k]])
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
                cid = _fresh_label(G, "rdnet_crossing", next_cross)
                next_cross = int(cid.rsplit("_", 1)[1]) + 1
                G.add_node(cid, **{**_between(self.graph.nodes[u], self.graph.nodes[v], t),
                                   self.x_attr: point.x, self.y_attr: point.y, "rdnet_crossing": True})
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
        self._split = (boundary, out)
        return out

    @property
    def needs_point_nodes(self) -> bool:
        """Whether locations on edges must become nodes for exact distances.

        A location partway along an edge is normally reached through the edge's
        end nodes plus part of the edge's cost. That ignores one-way edges and
        the penalty for switching onto the location's edge class, so directed
        networks and transition costs place each location on a node instead.
        """
        if not self.metric.supports_edge_splitting:
            return False
        directed = self.graph.is_directed() and self.direction != "both"
        return directed or hasattr(self.metric, "transition_penalty")

    def _with_point_nodes(self, u, v, frac, extensive: Sequence[str] = ("length", "travel_time")):
        """Insert a node at each location ``frac`` of the way from node ``u`` to node ``v``.

        Every stored edge between ``u`` and ``v`` (both directions, parallel
        edges) is split there, with its extensive attributes and the metric's
        cost attribute divided in proportion. Returns the new network and the
        node index of each location (an end node when ``frac`` is 0 or 1).
        """
        u = np.asarray(u, dtype=np.int64)
        v = np.asarray(v, dtype=np.int64)
        frac = np.clip(np.asarray(frac, dtype=float), 0.0, 1.0)
        node_label = np.empty(len(u), dtype=object)
        at_u = (u == v) | (frac <= 1e-9)
        at_v = ~at_u & (frac >= 1 - 1e-9)
        node_label[at_u] = [self.nodes[k] for k in u[at_u]]
        node_label[at_v] = [self.nodes[k] for k in v[at_v]]
        mid = np.flatnonzero(~at_u & ~at_v)
        if mid.size == 0:
            return self, np.array([self.node_index[x] for x in node_label], dtype=np.int64)
        extensive = set(extensive) | _cost_attributes(self.metric)
        # positions measured from the smaller node index, so that a location
        # snapped as (u, v, f) and one snapped as (v, u, 1 - f) coincide
        lo, hi = np.minimum(u[mid], v[mid]), np.maximum(u[mid], v[mid])
        t = np.where(u[mid] == lo, frac[mid], 1 - frac[mid])
        cuts: dict[tuple[int, int], dict[float, Hashable]] = {}
        G = self.graph.copy()
        n_new = 0
        for k, a, b, tk in zip(mid, lo, hi, t):
            spots = cuts.setdefault((int(a), int(b)), {})
            key = round(float(tk), 12)
            if key not in spots:
                xa, xb = self.coords[a], self.coords[b]
                xy = xa + (xb - xa) * key
                pid = _fresh_label(G, "rdnet_point", n_new)
                n_new = int(pid.rsplit("_", 1)[1]) + 1
                G.add_node(pid, **{**_between(self.graph.nodes[self.nodes[a]], self.graph.nodes[self.nodes[b]], key),
                                   self.x_attr: float(xy[0]), self.y_attr: float(xy[1]), "rdnet_point": True})
                spots[key] = pid
            node_label[k] = spots[key]
        multi = self.graph.is_multigraph()
        for (a, b), spots in cuts.items():
            la, lb = self.nodes[a], self.nodes[b]
            ts = sorted(spots)
            stored = []
            for p, q in ((la, lb), (lb, la)) if G.is_directed() else ((la, lb),):
                if not G.has_edge(p, q):
                    continue
                if multi:
                    stored += [(p, q, key, dict(d)) for key, d in G[p][q].items()]
                else:
                    stored.append((p, q, None, dict(G[p][q])))
            for p, q, key, data in stored:
                forward = p == la  # edge runs from the smaller-index node
                pos = ts if forward else [1 - x for x in reversed(ts)]
                ids = [spots[x] for x in (ts if forward else list(reversed(ts)))]
                geom = data.get("geometry")
                if geom is not None:
                    gp = shapely.get_point(geom, 0)
                    if Point(self.graph.nodes[p][self.x_attr], self.graph.nodes[p][self.y_attr]).distance(
                            gp) > Point(self.graph.nodes[q][self.x_attr], self.graph.nodes[q][self.y_attr]).distance(gp):
                        geom = shapely.reverse(geom)
                if multi:
                    G.remove_edge(p, q, key)
                else:
                    G.remove_edge(p, q)
                chain = [p, *ids, q]
                cut = [0.0, *pos, 1.0]
                for i in range(len(chain) - 1):
                    share = cut[i + 1] - cut[i]
                    piece = dict(data)
                    for attr in extensive:
                        if attr in data and data[attr] is not None:
                            piece[attr] = float(data[attr]) * share
                    if geom is not None:
                        piece["geometry"] = substring(geom, cut[i], cut[i + 1], normalized=True)
                    G.add_edge(chain[i], chain[i + 1], **piece)
        out = self._replace_graph(G)
        return out, np.array([out.node_index[x] for x in node_label], dtype=np.int64)

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
        dist: float = 1000.0,
        min_area: float = 250_000.0,
        **kwargs,
    ) -> "NetworkGraph":
        """Download an OpenStreetMap network with OSMnx.

        ``query`` is geocoded with Nominatim and its top match is used. If
        that match is an area of at least ``min_area`` square metres (a city,
        borough or district), the streets inside its outline are downloaded.
        Otherwise (many neighbourhoods, such as Midtown Manhattan, are mapped
        as points, and a query can also match a single building) the streets
        within ``dist`` metres of it are downloaded, with a warning.

        If OSMnx is missing or the download fails and ``offline_fallback`` is
        true, the bundled mock Midtown network is returned with a warning, so
        scripts and tests keep working without internet access.
        """
        try:
            import osmnx as ox

            if bbox is not None:
                G = ox.graph_from_bbox(bbox, network_type=network_type)
            elif query is not None:
                # Top match whatever its geometry: osmnx's default skips points and
                # can land on a much smaller feature (a building) further down the list.
                place = ox.geocode_to_gdf(query, which_result=1)
                geom = place.geometry.iloc[0]
                is_area = geom.geom_type in ("Polygon", "MultiPolygon")
                area = float(ox.projection.project_gdf(place).area.iloc[0]) if is_area else 0.0
                if area >= min_area:
                    G = ox.graph_from_polygon(geom, network_type=network_type)
                else:
                    name = place["display_name"].iloc[0] if "display_name" in place else query
                    what = f"an area of {area:,.0f} m2" if is_area else f"a {geom.geom_type}"
                    warnings.warn(f"{query!r} geocodes to {what} ({name}); using the streets within "
                                  f"{dist:g} m of it", stacklevel=2)
                    point = (float(place["lat"].iloc[0]), float(place["lon"].iloc[0]))
                    G = ox.graph_from_point(point, dist=dist, network_type=network_type)
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
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"adjacency must be square, got shape {A.shape}")
        if coords.shape != (A.shape[0], 2):
            raise ValueError(f"coords must be ({A.shape[0]}, 2), one row per node, got shape {coords.shape}")
        nodes = list(nodes) if nodes is not None else list(range(A.shape[0]))
        if len(nodes) != A.shape[0]:
            raise ValueError(f"nodes must have {A.shape[0]} labels, got {len(nodes)}")
        G = nx.DiGraph() if directed else nx.Graph()
        for n, (x, y) in zip(nodes, coords):
            G.add_node(n, x=x, y=y)
        for i, j, w in zip(A.row, A.col, A.data):
            if i == j:
                continue
            G.add_edge(nodes[i], nodes[j], weight=float(w), length=float(np.hypot(*(coords[i] - coords[j]))))
        return cls(G, metric, **kwargs)

    @classmethod
    def from_lines(
        cls,
        lines,
        metric="length",
        *,
        tol: float = 1e-6,
        split_at_intersections: bool = False,
        directed: bool = False,
        coords: Sequence[str] = ("x1", "y1", "x2", "y2"),
        to_crs="utm",
        crs=None,
        **kwargs,
    ) -> "NetworkGraph":
        """Build from line geometries, such as a GIS layer of roads, rivers, rail or pipes.

        Parameters
        ----------
        lines
            A GeoDataFrame or GeoSeries of (Multi)LineStrings, a list of shapely
            lines or ``(m, 2)`` vertex arrays, or a DataFrame of straight segments
            with endpoint columns ``coords`` (``x1, y1, x2, y2``). Other columns
            of a (Geo)DataFrame (speed, lanes, flow, ...) become edge attributes,
            so they can be used in a metric.
        tol
            Line ends closer than this are joined into one node. Increase it for
            layers whose lines do not quite meet.
        split_at_intersections
            Also split lines where they cross each other, for layers that are not
            split at junctions. Leave it off when crossings without a junction
            (bridges, tunnels) matter.
        directed
            Lines are drawn in the direction of travel or flow (one-way streets,
            rivers). Combine with ``direction="out"`` (downstream) or ``"in"``.
        to_crs
            For lines in longitude/latitude (a GeoDataFrame in a geographic CRS,
            or plain coordinates with a geographic ``crs``): ``"utm"`` (default)
            projects them to their local UTM zone so distances are in metres;
            any other CRS projects to that; ``None`` keeps degrees.
        crs
            The CRS of plain (non-geopandas) coordinates, if known.
        **kwargs
            Passed to :class:`NetworkGraph` (``direction``, ``snap_scale``).

        Each edge gets ``geometry`` and ``length`` (the line's length, unless a
        ``length`` column is supplied). Closed rings are split in two so they
        can be routed through.
        """
        geoms, attrs, crs = _lines_input(lines, coords, to_crs, crs)
        if split_at_intersections:
            geoms, attrs = _node_lines(geoms, attrs)
        geoms = np.asarray(shapely.force_2d(geoms), dtype=object)
        ok = ~shapely.is_empty(geoms) & (shapely.length(geoms) > 0)
        if not ok.all():
            warnings.warn(f"{int((~ok).sum())} empty or zero-length lines were dropped", stacklevel=2)
            geoms, attrs = geoms[ok], [a for a, k in zip(attrs, ok) if k]
        # closed rings: split at the midpoint so both halves are routable
        rings = shapely.equals(shapely.get_point(geoms, 0), shapely.get_point(geoms, -1))
        if rings.any():
            out_g, out_a = [], []
            for g, a, r in zip(geoms, attrs, rings):
                parts = [substring(g, 0, 0.5, normalized=True), substring(g, 0.5, 1, normalized=True)] if r else [g]
                out_g += parts
                out_a += [a] * len(parts)
            geoms, attrs = np.asarray(out_g, dtype=object), out_a
        ends = np.vstack([shapely.get_coordinates(shapely.get_point(geoms, 0)),
                          shapely.get_coordinates(shapely.get_point(geoms, -1))])
        label = _merge_close(ends, tol)
        m = len(geoms)
        G = nx.MultiDiGraph() if directed else nx.MultiGraph()
        G.graph["crs"] = crs
        first = {}
        for k, lab in enumerate(label):
            first.setdefault(lab, k)
        for lab, k in sorted(first.items()):
            G.add_node(int(lab), x=float(ends[k, 0]), y=float(ends[k, 1]))
        lengths = shapely.length(geoms)
        for k in range(m):
            data = dict(attrs[k])
            data.setdefault("length", float(lengths[k]))
            data["geometry"] = geoms[k]
            G.add_edge(int(label[k]), int(label[m + k]), **data)
        return cls(G, metric, crs=crs, **kwargs)

    @classmethod
    def from_frames(
        cls,
        nodes,
        edges,
        metric="length",
        *,
        node_id: str = "id",
        x: str = "x",
        y: str = "y",
        source: str = "from",
        target: str = "to",
        directed: bool = False,
        crs=None,
        **kwargs,
    ) -> "NetworkGraph":
        """Build from a node table and an edge table (the format most GIS and
        network tools export).

        ``nodes`` has an id column ``node_id`` and coordinates ``x``, ``y``
        (if ``node_id`` is not a column, the index is used). ``edges`` has the
        endpoint ids ``source`` and ``target``; every other column becomes an
        edge attribute. Edges without a ``length`` get the straight-line
        length between their endpoints. With ``directed=True`` edges run from
        ``source`` to ``target``; combine with ``direction="out"``.
        """
        nodes = pd.DataFrame(nodes)
        edges = pd.DataFrame(edges)
        for col in (x, y):
            if col not in nodes.columns:
                raise KeyError(f"nodes has no {col!r} column (set x= and y= to the coordinate columns)")
        for col in (source, target):
            if col not in edges.columns:
                raise KeyError(f"edges has no {col!r} column (set source= and target= to the endpoint columns)")
        ids = nodes[node_id].to_numpy() if node_id in nodes.columns else nodes.index.to_numpy()
        if pd.Index(ids).has_duplicates:
            raise ValueError("node ids must be unique")
        xy = nodes[[x, y]].to_numpy(dtype=float)
        if not np.isfinite(xy).all():
            raise ValueError("node coordinates must be finite")
        pos = {k: i for i, k in enumerate(ids.tolist())}
        missing = set(edges[source].tolist()) | set(edges[target].tolist())
        missing -= pos.keys()
        if missing:
            raise KeyError(f"edges refer to node ids missing from nodes: {sorted(map(str, missing))[:5]}")
        G = nx.MultiDiGraph() if directed else nx.MultiGraph()
        G.graph["crs"] = crs
        for k, (px, py) in zip(ids.tolist(), xy):
            G.add_node(k, x=float(px), y=float(py))
        extra = [c for c in edges.columns if c not in (source, target)]
        rows = edges[extra].to_dict("records") if extra else [{} for _ in range(len(edges))]
        for a, b, data in zip(edges[source].tolist(), edges[target].tolist(), rows):
            data = {k: v for k, v in data.items() if not (np.isscalar(v) and pd.isna(v))}
            if "length" not in data:
                data["length"] = float(np.hypot(*(xy[pos[a]] - xy[pos[b]])))
            G.add_edge(a, b, **data)
        return cls(G, metric, crs=crs, **kwargs)

    @classmethod
    def from_igraph(cls, g, metric="weight", *, x_attr: str = "x", y_attr: str = "y", **kwargs) -> "NetworkGraph":
        """Build from a python-igraph graph whose vertices carry coordinates.

        With the default ``metric="weight"``, a graph without a ``weight`` edge
        attribute uses ``length`` instead.
        """
        G = g.to_networkx()
        if metric == "weight" and "weight" not in g.es.attributes() and "length" in g.es.attributes():
            metric = "length"
        return cls(G, metric, x_attr=x_attr, y_attr=y_attr, **kwargs)


def _fresh_label(G, prefix: str, k: int) -> str:
    """``prefix_k``, or the next free ``prefix_k'`` when the graph already has that node.

    Networks that were split before (a coordinates object's network, reused for a
    placebo boundary) already carry such labels.
    """
    while f"{prefix}_{k}" in G:
        k += 1
    return f"{prefix}_{k}"


def _between(a: dict, b: dict, t: float) -> dict:
    """Node attributes for a node a fraction ``t`` of the way from a node with ``a`` to one with ``b``.

    Numeric attributes both ends share are interpolated, so that costs built
    from node attributes (``node_contrast``) add up along a split edge.
    """
    out = {}
    for k, x in a.items():
        y = b.get(k)
        if isinstance(x, (int, float, np.number)) and isinstance(y, (int, float, np.number))                 and not isinstance(x, bool) and not isinstance(y, bool):
            out[k] = float(x) + (float(y) - float(x)) * t
    return out


def _cost_attributes(metric) -> set:
    """Edge attributes a metric reads as additive costs (split in proportion when an edge is cut)."""
    out = set()
    for m in (metric, getattr(metric, "_base", None)):
        cost = getattr(m, "cost", None)
        if isinstance(cost, str):
            out.add(cost)
    return out


def _is_projected(graph) -> bool:
    crs = graph.graph.get("crs")
    if crs is None:
        return True
    return "4326" not in str(crs) and "longlat" not in str(crs).lower()


def _looks_lonlat(xy: np.ndarray, crs) -> bool:
    """A geographic CRS, or (with no CRS) coordinates that look like lon/lat degrees."""
    if crs is not None:
        return not _crs_projected(crs)
    if len(xy) < 2 or not np.isfinite(xy).all():
        return False
    span = np.ptp(xy, axis=0)
    return bool(np.all(np.abs(xy[:, 0]) <= 180) and np.all(np.abs(xy[:, 1]) <= 90)
                and span.max() < 2 and np.abs(xy).max() > 2 and np.any(xy != np.round(xy)))


def _crs_projected(crs) -> bool:
    text = str(crs)
    return "4326" not in text and "longlat" not in text.lower() and "4269" not in text


def project_lonlat(lon, lat, to_crs) -> np.ndarray:
    """Project longitude/latitude to ``to_crs`` (requires pyproj)."""
    try:
        from pyproj import Transformer
    except ImportError as e:  # pragma: no cover - depends on environment
        raise ImportError("projecting lon/lat requires pyproj (pip install pyproj)") from e
    tr = Transformer.from_crs("EPSG:4326", to_crs, always_xy=True)
    x, y = tr.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.column_stack([x, y])


def _project_plain_lines(geoms, crs, to_crs):
    """Project lines given as plain coordinates in a geographic ``crs`` (as GeoDataFrames are)."""
    try:
        from pyproj import CRS, Transformer
    except ImportError as e:  # pragma: no cover - depends on environment
        raise ImportError("projecting lon/lat lines requires pyproj (pip install pyproj)") from e
    src = CRS.from_user_input(crs)
    if not src.is_geographic:
        return geoms, src.to_string()
    if to_crs == "utm":
        lon, lat = shapely.get_coordinates(geoms).mean(axis=0)
        zone = int((lon + 180) // 6) + 1
        target = CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)
    else:
        target = CRS.from_user_input(to_crs)
    tr = Transformer.from_crs(src, target, always_xy=True)
    out = shapely.transform(geoms, lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1])))
    return out, target.to_string()


def _lines_input(lines, coords, to_crs, crs):
    """Geometries, per-line attribute dicts and CRS from any supported line input."""
    if hasattr(lines, "geometry") and hasattr(lines, "crs"):  # GeoDataFrame / GeoSeries
        gdf = lines if hasattr(lines, "columns") else lines.to_frame("geometry")
        crs = gdf.crs
        if crs is not None and crs.is_geographic and to_crs is not None:
            target = gdf.estimate_utm_crs() if to_crs == "utm" else to_crs
            gdf = gdf.to_crs(target)
            crs = gdf.crs
        gdf = gdf.explode(index_parts=False) if (gdf.geom_type == "MultiLineString").any() else gdf
        geoms = np.asarray(gdf.geometry.to_numpy(), dtype=object)
        rest = pd.DataFrame(gdf.drop(columns=gdf.geometry.name))
        attrs = rest.to_dict("records") if rest.shape[1] else [{} for _ in range(len(gdf))]
        crs = None if crs is None else crs.to_string()
    elif isinstance(lines, pd.DataFrame):
        missing = [c for c in coords if c not in lines.columns]
        if missing:
            raise KeyError(f"segment table has no {missing} columns (set coords= to the endpoint columns)")
        xy = lines[list(coords)].to_numpy(dtype=float)
        geoms = shapely.linestrings(xy.reshape(-1, 2, 2))
        rest = lines.drop(columns=list(coords))
        attrs = rest.to_dict("records") if rest.shape[1] else [{} for _ in range(len(lines))]
    else:
        geoms = np.asarray([g if isinstance(g, shapely.Geometry) else LineString(np.asarray(g, float)) for g in lines],
                           dtype=object)
        parts = shapely.get_parts(geoms, return_index=True)
        geoms, idx = parts
        attrs = [{} for _ in idx]
    attrs = [{k: v for k, v in a.items() if not (np.isscalar(v) and pd.isna(v))} for a in attrs]
    if not (hasattr(lines, "geometry") and hasattr(lines, "crs")) and crs is not None and to_crs is not None:
        geoms, crs = _project_plain_lines(np.asarray(geoms, dtype=object), crs, to_crs)
    bad = ~np.isin(shapely.get_type_id(geoms), (1, 2))  # LineString, LinearRing
    if bad.any():
        raise TypeError("every geometry must be a LineString or MultiLineString")
    return np.asarray(geoms, dtype=object), attrs, crs


def _node_lines(geoms, attrs):
    """Split lines at every crossing, keeping each piece's source attributes."""
    pieces = shapely.get_parts(shapely.unary_union(geoms))
    tree = shapely.STRtree(geoms)
    mid = shapely.line_interpolate_point(pieces, 0.5, normalized=True)
    src = tree.nearest(mid)
    out_a = []
    for p, k in zip(pieces, src):
        a = dict(attrs[k])
        if "length" in a:  # a supplied length is for the whole line; rescale to the piece
            a["length"] = float(a["length"]) * float(p.length) / float(geoms[k].length)
        out_a.append(a)
    return np.asarray(pieces, dtype=object), out_a


def _merge_close(xy: np.ndarray, tol: float) -> np.ndarray:
    """Label points so that points within ``tol`` of each other share a label."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    pairs = cKDTree(xy).query_pairs(max(tol, 0.0), output_type="ndarray")
    n = len(xy)
    A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)) if len(pairs) else coo_matrix((n, n))
    _, lab = connected_components(A, directed=False)
    # relabel in order of first appearance
    _, first, inv = np.unique(lab, return_index=True, return_inverse=True)
    order = np.argsort(np.argsort(first))
    return order[inv]
