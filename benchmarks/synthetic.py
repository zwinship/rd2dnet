"""Large synthetic networks and observation designs for scale testing.

Not part of the package API. Three network families:

``grid_city``
    Manhattan-style grid (default 250 x 200 = 50,000 nodes): jittered nodes,
    a share of segments removed, arterials every 10th street/avenue, a
    motorway with ramps, a park on the boundary (masked). The boundary runs
    along a street, so its nodes lie on it.
``suburb``
    Irregular planar network: Delaunay triangulation of random points, long
    edges and a share of the rest removed (the minimum spanning tree is always
    kept, so it stays connected), giving cul-de-sacs. The boundary is a sine
    curve that cuts edges at arbitrary points.
``river``
    Dendritic river network: minimum spanning tree of random points, edges
    directed downstream towards an outlet. The boundary is a jurisdiction line
    across the basin. Distance is along the channel.

Observation designs: ``parcels`` (uniform over the area), ``gps`` (on the
network with metre-scale jitter), ``panel`` (``n_sites`` fixed locations, each
observed ``n_obs // n_sites`` times: high-frequency sensor data).
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import Delaunay
from shapely.geometry import LineString, Polygon, box

from rd2dnet import Boundary

SPEED_KPH = {"residential": 25.0, "primary": 45.0, "motorway": 90.0, "motorway_link": 40.0, "river": 3.6}


@dataclass
class Synthetic:
    name: str
    graph: nx.Graph
    boundary: Boundary
    bbox: tuple[float, float, float, float]
    holes: list  # polygons with no observations (parks)
    edges: np.ndarray  # (E, 4): x0, y0, x1, y1 for on-network sampling

    def tau(self, arc):
        arc = np.asarray(arc, dtype=float)
        return 1.0 + 0.5 * np.sin(2 * np.pi * arc / self.boundary.length)


def _add_edge(G, a, b, xy, cls):
    length = float(np.hypot(*(xy[a] - xy[b])))
    speed = SPEED_KPH[cls]
    G.add_edge(int(a), int(b), length=max(length, 1e-3), highway=cls, speed_kph=speed,
               travel_time=length / (speed / 3.6))


def _finish(name, G, xy, boundary, bbox, holes):
    for i in G.nodes:
        G.nodes[i]["x"], G.nodes[i]["y"] = float(xy[i, 0]), float(xy[i, 1])
    e = np.array([(xy[u, 0], xy[u, 1], xy[v, 0], xy[v, 1]) for u, v in G.edges()])
    return Synthetic(name, G, boundary, bbox, holes, e)


def grid_city(n_streets: int = 250, n_avenues: int = 200, street: float = 80.0, avenue: float = 120.0,
              jitter: float = 8.0, drop: float = 0.05, seed: int = 0) -> Synthetic:
    rng = np.random.default_rng(seed)
    nid = lambda i, j: i * n_streets + j  # noqa: E731  (avenue i, street j)
    xy = np.zeros((n_avenues * n_streets + n_streets, 2))
    ii, jj = np.meshgrid(np.arange(n_avenues), np.arange(n_streets), indexing="ij")
    xy[: n_avenues * n_streets, 0] = ii.ravel() * avenue
    xy[: n_avenues * n_streets, 1] = jj.ravel() * street
    jb = n_streets // 2  # boundary street
    by = jb * street
    off = np.abs(jj.ravel() - jb) > 0  # keep the boundary street straight
    xy[: n_avenues * n_streets][off] += rng.normal(0, jitter, (off.sum(), 2))
    park_i = (n_avenues // 3, n_avenues // 3 + 8)
    park_j = (jb, jb + 25)
    in_park = lambda i, j: park_i[0] < i < park_i[1] and park_j[0] < j < park_j[1]  # noqa: E731
    G = nx.Graph(name="grid_city", crs="local-metres")
    G.add_nodes_from(nid(i, j) for i in range(n_avenues) for j in range(n_streets) if not in_park(i, j))
    keep = rng.uniform(size=(n_avenues, n_streets, 2)) > drop
    for i in range(n_avenues):
        for j in range(n_streets):
            if in_park(i, j):
                continue
            if i + 1 < n_avenues and not in_park(i + 1, j) and (keep[i, j, 0] or j == jb):
                _add_edge(G, nid(i, j), nid(i + 1, j), xy, "primary" if j % 10 == 0 else "residential")
            if j + 1 < n_streets and not in_park(i, j + 1) and keep[i, j, 1]:
                _add_edge(G, nid(i, j), nid(i, j + 1), xy, "primary" if i % 10 == 0 else "residential")
    hw0 = n_avenues * n_streets
    for j in range(n_streets):
        xy[hw0 + j] = (-300.0, j * street)
        G.add_node(hw0 + j)
        if j:
            _add_edge(G, hw0 + j - 1, hw0 + j, xy, "motorway")
        if j % 12 == 5:
            _add_edge(G, hw0 + j, nid(0, j), xy, "motorway_link")
    G.remove_nodes_from([n for n in list(G.nodes) if G.degree(n) == 0])
    W = (n_avenues - 1) * avenue
    park = box(park_i[0] * avenue, by, park_i[1] * avenue, park_j[1] * street)
    bnd = Boundary(LineString([(-400.0, by), (W + 100.0, by)]),
                   treated_region=box(-1e6, -1e6, 1e6, by), mask=[park])
    return _finish("grid_city", G, xy, bnd, (0.0, 0.0, W, (n_streets - 1) * street), [park])


def _delaunay_edges(xy):
    tri = Delaunay(xy)
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    return e, np.hypot(*(xy[e[:, 0]] - xy[e[:, 1]]).T)


def suburb(n_nodes: int = 50_000, density_per_km2: float = 400.0, drop: float = 0.35, seed: int = 0) -> Synthetic:
    rng = np.random.default_rng(seed)
    side = np.sqrt(n_nodes / density_per_km2) * 1000.0
    xy = rng.uniform(0, side, (n_nodes, 2))
    e, L = _delaunay_edges(xy)
    A = sparse.coo_matrix((L, (e[:, 0], e[:, 1])), shape=(n_nodes, n_nodes))
    mst = minimum_spanning_tree(A).tocoo()
    in_mst = set(zip(np.minimum(mst.row, mst.col), np.maximum(mst.row, mst.col)))
    long = L > np.quantile(L, 0.9)
    q80 = np.quantile(L, 0.8)
    G = nx.Graph(name="suburb", crs="local-metres")
    G.add_nodes_from(range(n_nodes))
    u = rng.uniform(size=len(e))
    for k, (a, b) in enumerate(e):
        if (a, b) in in_mst or (not long[k] and u[k] > drop):
            _add_edge(G, a, b, xy, "primary" if L[k] > q80 else "residential")
    xs = np.linspace(-50, side + 50, 200)
    ys = side / 2 + 0.08 * side * np.sin(2 * np.pi * xs / side * 1.5)
    line = LineString(np.column_stack([xs, ys]))
    region = Polygon(list(zip(xs, ys)) + [(side + 50, -1e6), (-50, -1e6)])
    return _finish("suburb", G, xy, Boundary(line, treated_region=region), (0.0, 0.0, side, side), [])


def river(n_nodes: int = 50_000, width: float = 40_000.0, height: float = 60_000.0, seed: int = 0) -> Synthetic:
    rng = np.random.default_rng(seed)
    xy = rng.uniform((0, 0), (width, height), (n_nodes, 2))
    xy[0] = (width / 2, 0.0)  # outlet
    e, L = _delaunay_edges(xy)
    A = sparse.coo_matrix((L, (e[:, 0], e[:, 1])), shape=(n_nodes, n_nodes))
    T = minimum_spanning_tree(A)
    T = T + T.T
    # Orient every reach downstream (towards the outlet) by BFS from the outlet.
    from scipy.sparse.csgraph import breadth_first_order

    _, pred = breadth_first_order(T, 0, directed=False, return_predecessors=True)
    G = nx.DiGraph(name="river", crs="local-metres")
    G.add_nodes_from(range(n_nodes))
    for child in range(1, n_nodes):
        _add_edge(G, child, int(pred[child]), xy, "river")
    by = height * 0.45
    bnd = Boundary(LineString([(-100.0, by), (width + 100.0, by)]), treated_region=box(-1e6, -1e6, 1e6, by))
    return _finish("river", G, xy, bnd, (0.0, 0.0, width, height), [])


def observations(net: Synthetic, n: int, design: str = "parcels", n_sites: int = 2000, seed: int = 1):
    """Return ``(xy, y)``: locations and an outcome with effect ``net.tau(arc)`` on the treated side."""
    rng = np.random.default_rng(seed)
    x0, y0, x1, y1 = net.bbox
    m = n if design != "panel" else n_sites
    if design == "parcels":
        pts = []
        while sum(len(p) for p in pts) < m:
            c = rng.uniform((x0, y0), (x1, y1), (m, 2))
            for h in net.holes:
                c = c[~((c[:, 0] > h.bounds[0]) & (c[:, 0] < h.bounds[2]) & (c[:, 1] > h.bounds[1]) & (c[:, 1] < h.bounds[3]))]
            pts.append(c)
        xy = np.vstack(pts)[:m]
    else:  # gps / panel sites: on the network with metre-scale jitter
        E = net.edges
        L = np.hypot(E[:, 2] - E[:, 0], E[:, 3] - E[:, 1])
        k = rng.choice(len(E), size=m, p=L / L.sum())
        t = rng.uniform(size=m)[:, None]
        xy = E[k, :2] * (1 - t) + E[k, 2:] * t + rng.normal(0, 5.0, (m, 2))
    if design == "panel":
        xy = np.repeat(xy, int(np.ceil(n / n_sites)), axis=0)[:n]
    T = net.boundary.is_treated(xy)
    arc = net.boundary.project(xy)
    scale = max(x1 - x0, y1 - y0)
    f = 0.5 * np.cos(2 * np.pi * xy[:, 0] / scale) + 0.3 * (xy[:, 1] - net.boundary.line.centroid.y) / scale
    y = f + net.tau(arc) * T + rng.normal(0, 1.0, len(xy))
    return xy, y
