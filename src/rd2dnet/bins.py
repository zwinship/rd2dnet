"""Nearest-crossing bins: boundary discontinuity estimates per boundary crossing.

An alternative to the default location-based design (:func:`rdnet` on the
output of :func:`network_coordinates`, which pools every observation near an
evaluation point and fits a polynomial in ``(s, o)``). Here, following
Winship, *Network Boundary Discontinuity Designs via Nearest-Intersection
Partitioning*:

1. **Bins.** Every observation belongs to the bin of the boundary crossing
   nearest to it along the network (``coords.crossing_node``); its signed
   distance to that crossing is ``R = s``. Where a street runs along the
   boundary, each of its nodes is a crossing, and locations on that street
   belong to the nearer one.
2. **Edge similarity.** Each street segment gets a feature vector: edge
   attributes you name (capacity, speed, ...) and its local topology (adjacent
   edges, edges two steps away, the larger end degree, mean clustering of its
   ends), computed on the network before it was split at the boundary. ``S``
   is the Mahalanobis distance between the features of an observation's
   segment and those of the crossing's segment(s), with the covariance of the
   features over all segments.
3. **Bin fits.** :func:`rdnet` on :func:`bin_coordinates` fits each bin
   separately, weighting by the product kernel ``K1(R / h_R) * K2(S / h_S)``.
   Branches (the edge where the path from the crossing first passes a
   junction) get their own slopes in ``R`` on each side, and ``R * S`` enters
   with a slope per side; the intercepts, and so the jump, are common to the
   bin. Bias correction, bandwidth selection (``h_R``), cluster and Conley
   errors, fuzzy designs and covariates work as in the default design.
4. **Aggregation.** ``res.wbate("density")`` weights bins by their local
   sample size, ``res.wbate("gls")`` by the inverse of the joint covariance;
   both use ``se = sqrt(w' cov w)``.

:func:`rdnet_fe` pools all bins in one regression with bin fixed effects and
inverse-propensity weights, for boundaries with many sparsely populated
crossings. :func:`similarity_sensitivity` reports the aggregate effect over a
range of ``h_S``.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.sparse.csgraph import dijkstra

from ._wls import poly_basis, variance, wls
from .bandwidth import _check_cluster, _check_covs, _check_obs, _check_weights, rdbwnet
from .coordinates import NetworkCoordinates
from .kernels import get_kernel
from .metrics import _min_csr, _scalar
from .network import Snap

TOPOLOGY_FEATURES = ("n_adjacent", "n_second", "degree_max", "clustering")


def similarity_bandwidth(k: int, level: float = 0.95) -> float:
    """Rule-of-thumb ``h_S``: ``sqrt(2 * chi2_k(level))``.

    If two segments' features were independent draws with covariance
    ``Omega``, their difference would have covariance ``2 Omega`` and ``S**2 / 2``
    would follow a chi-squared distribution with ``k`` degrees of freedom
    (``k`` the rank of ``Omega``). The default keeps all but the 5% most
    dissimilar pairs under that reference, with weights falling linearly in
    ``S`` (triangular ``K2``). It is a scale-free starting point, not an
    optimal choice; report :func:`similarity_sensitivity` alongside it.
    """
    if k <= 0:
        return np.inf
    return float(np.sqrt(2.0 * stats.chi2.ppf(level, k)))


class BinCoordinates(NetworkCoordinates):
    """Network coordinates with one evaluation point per nearest-crossing bin.

    Build with :func:`bin_coordinates` and pass to :func:`rdnet`,
    :func:`rdbwnet`, :func:`rdnet_fe` or :func:`similarity_sensitivity`.
    Evaluation point ``j`` is bin ``j``: its local sample is the bin's
    members, at kernel distance ``|s|`` (their distance to the crossing).

    Attributes
    ----------
    bin
        Bin of each observation (-1 if it cannot reach the boundary).
    branch
        Branch of each observation: the network edge where the path from the
        crossing first passes a junction (-1 before any junction: the trunk).
    similarity
        Edge dissimilarity ``S`` to the bin's crossing segment (NaN when
        ``similarity=False``).
    h_S, kernel_S
        Similarity bandwidth and kernel; ``kernel_weights`` is
        ``K2(S / h_S) / K2(0)``.
    features
        Feature names behind ``S``; ``reference`` holds each bin's crossing
        features (one row per bin).
    """

    def __repr__(self) -> str:
        return (f"BinCoordinates(n={self.n_obs}, bins={self.n_eval}, features={len(self.features)}, "
                f"h_S={self.h_S:.4g})")

    def at(self, j: int) -> tuple[np.ndarray, np.ndarray]:
        """Kernel distance ``|s|`` and offset 0 for the members of bin ``j`` (inf / NaN for everyone else)."""
        member = self.bin == j
        return np.where(member, np.abs(self.s), np.inf), np.where(member, 0.0, np.nan)

    def node_rows(self, js) -> None:
        """Nothing to precompute: kernel distances are ``|s|``."""

    def local_slopes(self, j: int, rows: np.ndarray) -> np.ndarray | None:
        """Extra slope regressors for bin ``j`` at ``rows``, multiplied by ``s / h`` in the fit.

        One indicator per branch of the bin (the trunk is the base; when every
        member is on some branch, the last indicator is aliased and dropped),
        then ``S``.
        """
        cols = [(self.branch[rows] == b).astype(float) for b in self._bin_branches[j]]
        if self.has_similarity:
            cols.append(self.similarity[rows])
        return np.column_stack(cols) if cols else None

    @property
    def has_similarity(self) -> bool:
        return self._use_similarity and len(self.features) > 0

    def with_similarity_bandwidth(self, h_S: float) -> "BinCoordinates":
        """The same bins with another ``h_S`` (``np.inf``: no similarity kernel)."""
        out = object.__new__(BinCoordinates)
        out.__dict__.update(self.__dict__)
        out._set_h_S(h_S)
        return out

    def _set_h_S(self, h_S) -> None:
        if not (np.isinf(h_S) or h_S > 0) or np.isnan(h_S):
            raise ValueError(f"h_S must be positive (np.inf for no similarity kernel), got {h_S!r}")
        self.h_S = float(h_S)
        if self.has_similarity and np.isfinite(self.h_S):
            k = get_kernel(self.kernel_S)
            k0 = float(k(0.0))
            S = np.where(np.isfinite(self.similarity), self.similarity, 0.0)
            self.kernel_weights = np.asarray(k(S / self.h_S), dtype=float) / k0
        else:
            self.kernel_weights = None

    def bins_table(self) -> pd.DataFrame:
        """One row per bin: position, crossing nodes, members per side, branches and crossing features."""
        T = self.treated
        out = pd.DataFrame({
            "bin": np.arange(self.n_eval),
            "arc": self.eval_arc,
            "node": [self.network.nodes[k] for k in self.eval_node],
            "n_crossings": [len(c) for c in self._bin_nodes],
            "n0": np.bincount(self.bin[(self.bin >= 0) & ~T], minlength=self.n_eval),
            "n1": np.bincount(self.bin[(self.bin >= 0) & T], minlength=self.n_eval),
            "n_branches": [len(b) for b in self._bin_branches],
            "masked": self.eval_masked,
        })
        for k, name in enumerate(self.features):
            out[name] = self.reference[:, k]
        return out


def bin_coordinates(
    coords: NetworkCoordinates,
    *,
    features=(),
    topology: bool = True,
    similarity: bool = True,
    branches: bool = True,
    h_S="rot",
    kernel_S="triangular",
) -> BinCoordinates:
    """Nearest-crossing bins for :func:`rdnet` (an alternative to the default design).

    Parameters
    ----------
    coords
        Output of :func:`network_coordinates` (its evaluation points are not used).
    features
        Edge attributes describing street segments (capacity, speed, lanes,
        baseline volume, ...). Use attributes that describe the street, not
        its piece: lengths are divided when edges are split, other attributes
        are copied. Missing values are filled with the attribute's mean.
    topology
        Add the local topology features ``n_adjacent`` (edges sharing an end),
        ``n_second`` (edges two steps away), ``degree_max`` (the larger end
        degree) and ``clustering`` (mean clustering coefficient of the ends),
        computed on the network as given (nodes inserted at the boundary and
        at locations are skipped).
    similarity
        Use the edge dissimilarity ``S``: the kernel ``K2(S / h_S)`` and the
        ``R * S`` slopes. ``False`` gives plain bins (Step 1 only).
    branches
        Give each branch its own slopes on each side.
    h_S
        Similarity bandwidth: ``"rot"`` (default, :func:`similarity_bandwidth`
        for the rank of the feature covariance), a positive number, or
        ``np.inf`` (no kernel in ``S``; the ``R * S`` slopes stay).
    kernel_S
        Kernel ``K2`` (default triangular, as ``K1``).

    Returns
    -------
    BinCoordinates
        One evaluation point per bin, ordered along the boundary. Use
        :meth:`BinCoordinates.bins_table` to inspect the bins.

    Notes
    -----
    Branches follow shortest paths over the metric's edge costs ignoring
    edge directions (exact for undirected cost metrics; with directed
    networks, transition penalties or precomputed metrics the tree is an
    approximation, and lengths are used where edge costs are not available).
    ``rdbwnet`` selects ``h_R`` per bin from the members' fit in ``R`` with the
    ``K2`` weights, without the branch and ``S`` slopes.

    There is no data-driven selector for ``h_S``: unlike ``R``, ``S`` is not a
    running variable extrapolated to a cutoff, and the MSE theory behind
    ``rdbwnet`` does not cover it. Report :func:`similarity_sensitivity`.
    """
    if isinstance(coords, BinCoordinates):
        raise ValueError("coords are already bins; pass the output of network_coordinates")
    if isinstance(features, str):
        features = [features]
    features = list(features)
    net = coords.network
    n_nodes = net.n_nodes
    E = net.edge_frame()
    eu, ev = E["u"].to_numpy(), E["v"].to_numpy()
    node_on = coords.node_on_boundary

    # one bin per crossing (boundary node) that some location reaches first
    first = coords._first
    cross = coords.crossing_node[first]  # per location
    ok_loc = (cross >= 0) & np.isfinite(coords.s[first])
    used = np.unique(cross[ok_loc])
    arc = coords.boundary.project(net.coords[used]) if used.size else np.zeros(0)
    order = np.argsort(arc, kind="stable")
    rep, arc = used[order], arc[order]
    M = len(rep)
    bin_node = np.full(n_nodes, -1)
    bin_node[rep] = np.arange(M)
    bin_loc = np.where(ok_loc, bin_node[np.maximum(cross, 0)], -1)
    bin_obs = bin_loc[coords.loc]
    bin_nodes = [rep[m:m + 1] for m in range(M)]
    med = arc

    tree = _crossing_tree(net, node_on, coords.node_treated, E)
    out = object.__new__(BinCoordinates)
    out.__dict__.update(coords.__dict__)
    out._rows = {}
    out.eval_arc = med
    out.eval_node_arc = med.copy()
    out.eval_node = rep
    out.eval_xy = net.coords[rep] if M else np.zeros((0, 2))
    out.eval_snap = Snap(rep, rep.copy(), np.zeros(M), np.zeros(M), np.zeros(M))
    out.eval_masked = coords.boundary.is_masked(med)
    out.bin = bin_obs
    out._bin_nodes = bin_nodes

    sn = coords._lsnap
    branch_loc = _location_branches(tree, sn, coords.treated[first]) if branches else np.full(len(first), -1)
    # a location whose tree path starts in another bin (ties, or a directed metric): trunk
    sc = tree["source"]
    tb = np.where(sc >= 0, bin_node[np.maximum(sc, 0)], -1)
    via_u = sn.cost_u + tree["dist"][sn.u]
    end = np.where(via_u <= sn.cost_v + tree["dist"][sn.v], sn.u, sn.v)
    other_start = (tb[end] != bin_loc) & (bin_loc >= 0)
    branch_loc = np.where(other_start, -1, branch_loc)
    out.branch = branch_loc[coords.loc]
    out._bin_branches = []
    for m in range(M):
        b = np.unique(branch_loc[(bin_loc == m) & (branch_loc >= 0)])
        out._bin_branches.append(b)

    names: list[str] = []
    if similarity and (features or topology):
        orig, Z, names = _edge_feature_matrix(net, E, list(features), topology)
        Omega = np.atleast_2d(np.cov(Z, rowvar=False)) if len(Z) > 1 else np.zeros((Z.shape[1], Z.shape[1]))
        Oinv = np.linalg.pinv(Omega, rcond=1e-10, hermitian=True)
        rank = int(np.linalg.matrix_rank(Omega, tol=1e-10 * max(np.abs(Omega).max(), 1e-300), hermitian=True)) \
            if Omega.size else 0
        # each bin's reference: mean features of the segments at its crossings
        ref = np.zeros((M, len(names)))
        inc_bin = np.r_[bin_node[eu], bin_node[ev]]
        inc_orig = np.r_[orig, orig]
        keep = inc_bin >= 0
        pairs = np.unique(np.column_stack([inc_bin[keep], inc_orig[keep]]), axis=0)
        for m in range(M):
            ids = pairs[pairs[:, 0] == m, 1]
            if ids.size:
                ref[m] = Z[ids].mean(axis=0)
        row_loc = _location_edge_rows(sn, tree)
        z_loc = np.where((row_loc >= 0)[:, None], Z[orig[np.maximum(row_loc, 0)]], np.nan)
        dz = z_loc - ref[np.maximum(bin_loc, 0)]
        dz = np.where((row_loc >= 0)[:, None], dz, 0.0)  # at the crossing itself: S = 0
        S_loc = np.sqrt(np.clip(np.einsum("ij,jk,ik->i", dz, Oinv, dz), 0.0, None))
        S_loc[bin_loc < 0] = np.nan
        out.similarity = S_loc[coords.loc]
        out.reference = ref
        out.omega = Omega
        out.feature_rank = rank
    else:
        out.similarity = np.full(coords.n_obs, np.nan)
        out.reference = np.zeros((M, 0))
        out.omega = np.zeros((0, 0))
        out.feature_rank = 0
    out.features = names
    out._use_similarity = bool(similarity)
    get_kernel(kernel_S)
    out.kernel_S = kernel_S
    if isinstance(h_S, str):
        if h_S != "rot":
            raise ValueError("h_S must be 'rot', a positive number or np.inf")
        h_S = similarity_bandwidth(out.feature_rank)
    out._set_h_S(h_S)
    if M == 0:
        warnings.warn("no observation reaches the boundary: there are no bins", stacklevel=2)
    return out


def _crossing_tree(net, node_on, node_T, E) -> dict:
    """Shortest-path tree from the boundary nodes and the branch of every node.

    A branch starts at the first junction on the path from the crossing: a
    node with three or more neighbours, or a crossing with two or more
    neighbours off the boundary on that side. Nodes before the first junction are on the
    trunk (-1). Branches are labelled by the edge-frame row of their first edge.
    """
    n = net.n_nodes
    eu, ev = E["u"].to_numpy(), E["v"].to_numpy()
    cost = E["cost"].to_numpy(dtype=float)
    length = E["length"].to_numpy(dtype=float)
    cost = np.where(np.isfinite(cost), cost, length)
    loop = eu == ev
    A = _min_csr(np.r_[eu[~loop], ev[~loop]], np.r_[ev[~loop], eu[~loop]], np.r_[cost[~loop], cost[~loop]], n)
    sources = np.flatnonzero(node_on)
    dist, pred, src = dijkstra(A, directed=True, indices=sources, min_only=True, return_predecessors=True)
    pred = np.where(pred < 0, -1, pred)
    # undirected neighbours (distinct, no self-loops)
    pairs = np.unique(np.column_stack([np.minimum(eu, ev), np.maximum(eu, ev)])[~loop], axis=0)
    deg = np.bincount(pairs.ravel(), minlength=n)
    # crossings: neighbours off the boundary on each side
    a, b = np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]]
    off = node_on[a] & ~node_on[b]
    side_T = np.bincount(a[off & node_T[b]], minlength=n)
    side_C = np.bincount(a[off & ~node_T[b]], minlength=n)
    full_key = np.minimum(eu, ev) * n + np.maximum(eu, ev)
    row_order = np.argsort(full_key, kind="stable")
    sorted_key = full_key[row_order]

    def row_of(p, x):
        """Edge-frame row of the segment between nodes ``p`` and ``x`` (arrays or scalars)."""
        p, x = np.asarray(p, dtype=np.int64), np.asarray(x, dtype=np.int64)
        k = np.minimum(p, x) * n + np.maximum(p, x)
        return row_order[np.minimum(np.searchsorted(sorted_key, k), len(sorted_key) - 1)]

    def splits(p, treated_side):
        """Whether a path leaving node ``p`` towards ``treated_side`` starts a branch (vectorised)."""
        p, ts = np.asarray(p, dtype=np.int64), np.asarray(treated_side, dtype=bool)
        return np.where(node_on[p], np.where(ts, side_T[p], side_C[p]) >= 2, deg[p] >= 3)

    # label[x]: the first edge on the path from the crossing to x that leaves a
    # junction, i.e. the topmost "start" edge among x's ancestors; found by pointer
    # jumping over the tree (each round doubles the stretch of path summarised)
    label = np.full(n, -1, dtype=np.int64)
    child = np.flatnonzero(pred >= 0)
    st = splits(pred[child], node_T[child])
    label[child[st]] = row_of(pred[child[st]], child[st])
    anc = pred.copy()
    while True:
        m = np.flatnonzero(anc >= 0)
        if m.size == 0:
            break
        up = anc[m]
        label[m] = np.where(label[up] >= 0, label[up], label[m])
        anc[m] = anc[up]
    return {"dist": dist, "pred": pred, "label": label, "splits": splits, "row_of": row_of,
            "source": np.where(src >= 0, src, -1)}


def _location_branches(tree, sn: Snap, T) -> np.ndarray:
    """Branch of each location (one per row of ``sn``)."""
    dist, label = tree["dist"], tree["label"]
    swap = sn.cost_u + dist[sn.u] > sn.cost_v + dist[sn.v]
    a = np.where(swap, sn.v, sn.u)  # the end the location is reached through
    b = np.where(swap, sn.u, sn.v)
    out = np.where(label[a] >= 0, label[a], -1)
    new = (out < 0) & (a != b) & tree["splits"](a, T)
    out[new] = tree["row_of"](a[new], b[new])
    return out


def _location_edge_rows(sn: Snap, tree) -> np.ndarray:
    """Edge-frame row of the segment each location sits on (-1 at a crossing node)."""
    pred, row_of = tree["pred"], tree["row_of"]
    out = np.full(len(sn.u), -1, dtype=np.int64)
    on_edge = sn.u != sn.v
    out[on_edge] = row_of(sn.u[on_edge], sn.v[on_edge])
    at_node = ~on_edge & (pred[sn.u] >= 0)  # on a node: the segment it was reached along
    out[at_node] = row_of(pred[sn.u[at_node]], sn.u[at_node])
    return out


def _edge_feature_matrix(net, E, attrs: list[str], topology: bool):
    """``(orig, Z, names)``: original-segment id of each edge-frame row and features per original segment.

    Nodes the package inserted (at boundary crossings and at locations) with
    two neighbours are skipped, so a segment split into pieces is one segment
    with the topology of the network as given.
    """
    n = net.n_nodes
    eu, ev = E["u"].to_numpy(), E["v"].to_numpy()
    deg_split = np.bincount(np.r_[eu[eu != ev], ev[eu != ev]], minlength=n)
    labels = net.nodes
    inserted = np.array([bool(net.graph.nodes[labels[k]].get("rdnet_crossing") or
                              net.graph.nodes[labels[k]].get("rdnet_point")) for k in range(n)], dtype=bool)
    transparent = inserted & (deg_split == 2)
    nb = np.full((n, 2), -1, dtype=np.int64)
    fill = np.zeros(n, dtype=np.int64)
    for a, b in zip(eu, ev):
        for x, y in ((a, b), (b, a)):
            if transparent[x] and fill[x] < 2:
                nb[x, fill[x]] = y
                fill[x] += 1

    def walk(a, prev):
        for _ in range(n):
            if not transparent[a]:
                return a
            nxt = nb[a, 1] if nb[a, 0] == prev else nb[a, 0]
            prev, a = a, nxt
        return a

    oa, ob = eu.copy(), ev.copy()
    for r in np.flatnonzero(transparent[eu] | transparent[ev]):
        oa[r], ob[r] = walk(eu[r], ev[r]), walk(ev[r], eu[r])
    okey = np.minimum(oa, ob) * n + np.maximum(oa, ob)
    ukeys, orig = np.unique(okey, return_inverse=True)
    m = len(ukeys)
    A_, B_ = ukeys // n, ukeys % n
    cols, names = [], []
    for attr in attrs:
        vals = np.full(len(E), np.nan)
        for r, (a, b) in enumerate(zip(eu, ev)):
            data = net.graph.get_edge_data(labels[a], labels[b])
            if data is None and net.graph.is_directed():
                data = net.graph.get_edge_data(labels[b], labels[a])
            if data is not None and net.graph.is_multigraph():
                data = next(iter(data.values()), None)
            x = None if data is None else data.get(attr)
            try:
                vals[r] = float(_scalar(x)) if x is not None else np.nan
            except (TypeError, ValueError):
                vals[r] = np.nan
        sums = np.bincount(orig, weights=np.where(np.isfinite(vals), vals, 0.0), minlength=m)
        cnt = np.bincount(orig, weights=np.isfinite(vals).astype(float), minlength=m)
        col = np.divide(sums, cnt, out=np.full(m, np.nan), where=cnt > 0)
        miss = ~np.isfinite(col)
        if miss.all():
            raise ValueError(f"edge attribute {attr!r} is missing on every segment")
        if miss.any():
            warnings.warn(f"edge attribute {attr!r} is missing on {int(miss.sum())} of {m} segments; filled with its "
                          "mean", stacklevel=3)
            col[miss] = np.nanmean(col)
        cols.append(col)
        names.append(attr)
    if topology:
        real = A_ != B_
        deg = np.bincount(np.r_[A_[real], B_[real]], minlength=n).astype(float)
        cols.append(np.where(real, deg[A_] + deg[B_] - 2, deg[A_]))
        names.append("n_adjacent")
        inc = sparse.csr_matrix((np.ones(2 * int(real.sum())), (np.r_[A_[real], B_[real]],
                                np.r_[np.flatnonzero(real), np.flatnonzero(real)])), shape=(n, m))
        adj = (inc.T @ inc).tocsr()
        adj.setdiag(0)
        adj.eliminate_zeros()
        adj.data[:] = 1.0
        two = (adj @ adj).tocsr()
        two.data[:] = 1.0
        two = two - two.multiply(adj)
        two.setdiag(0)
        two.eliminate_zeros()
        cols.append(np.asarray(two.getnnz(axis=1), dtype=float))
        names.append("n_second")
        cols.append(np.maximum(deg[A_], deg[B_]))
        names.append("degree_max")
        cols.append(_clustering(A_[real], B_[real], n)[[A_, B_]].mean(axis=0))
        names.append("clustering")
    Z = np.column_stack(cols) if cols else np.zeros((m, 0))
    return orig, Z, names


def _clustering(a: np.ndarray, b: np.ndarray, n: int) -> np.ndarray:
    """Local clustering coefficient of every node of a simple undirected graph (0 below degree 2)."""
    pairs = np.unique(np.column_stack([np.minimum(a, b), np.maximum(a, b)]), axis=0)
    A = sparse.csr_matrix((np.ones(2 * len(pairs)), (np.r_[pairs[:, 0], pairs[:, 1]],
                                                     np.r_[pairs[:, 1], pairs[:, 0]])), shape=(n, n))
    deg = np.asarray(A.sum(axis=1)).ravel()
    tri = np.asarray((A @ A).multiply(A).sum(axis=1)).ravel() / 2  # triangles through each node
    den = deg * (deg - 1) / 2
    return np.divide(tri, den, out=np.zeros(n), where=den > 0)


# --- pooled fixed-effects estimator -------------------------------------------------------


@dataclass
class RDNetFEResult:
    """Pooled boundary average treatment effect from :func:`rdnet_fe`."""

    tau: float
    se: float
    tau_bc: float
    se_rb: float
    ci_lo: float
    ci_hi: float
    pvalue: float
    level: float
    p: int
    q: int
    vce: str
    df: float
    bins: pd.DataFrame = field(repr=False)
    terms: list = field(default_factory=list, repr=False)

    def summary(self) -> str:
        used = self.bins["used"]
        lines = [
            "rd2dnet: pooled nearest-crossing BATE with bin fixed effects",
            f"  p={self.p}  q={self.q}  vce={self.vce}  bins used: {int(used.sum())} of {len(self.bins)}  "
            f"n0={int(self.bins.loc[used, 'n0'].sum())}  n1={int(self.bins.loc[used, 'n1'].sum())}",
            f"  slope terms: {', '.join(self.terms) if self.terms else 'none'}",
            "",
            pd.DataFrame([self.as_series()]).to_string(index=False, float_format=lambda v: f"{v:.3f}"),
        ]
        return "\n".join(lines)

    def as_series(self) -> pd.Series:
        return pd.Series(dict(tau=self.tau, se=self.se, tau_bc=self.tau_bc, se_rb=self.se_rb, ci_lo=self.ci_lo,
                              ci_hi=self.ci_hi, pvalue=self.pvalue), name="BATE")

    def __str__(self) -> str:
        return self.summary()


def rdnet_fe(
    y,
    bins: BinCoordinates,
    *,
    h=None,
    p: int = 1,
    q: int | None = None,
    kernel="triangular",
    vce: str = "hc3",
    cluster=None,
    cluster_vce: str = "cr2",
    conley_cutoff: float | None = None,
    conley_distance: str = "network",
    weights=None,
    covs=None,
    ipw: bool = True,
    bin_slopes: bool = True,
    level: float = 95.0,
    min_obs: int = 1,
    bwselect: str = "mserd",
    scaleregul: float = 3.0,
) -> RDNetFEResult:
    """Pooled BATE with nearest-crossing fixed effects, for many sparsely populated crossings.

    One weighted regression over the members of all bins::

        y = a_m + T_B * T + poly_p(R) + T * poly_p(R)
              + R * Z_m' th1 + T * R * Z_m' th2 + d1 * R * S + d2 * T * R * S + covs

    with ``T`` the 0/1 treatment indicator, ``R = s / h`` (each side's
    bandwidth), ``a_m`` one intercept per bin (so the jump is identified
    within bins only), ``Z_m`` the bin's crossing features (standardised
    across bins) and ``S`` the edge dissimilarity. Weights are the product
    kernel times, with ``ipw=True``, ``1 / p_m`` for treated and ``1 / (1 - p_m)``
    for control members, ``p_m`` the kernel-weighted treated share of bin
    ``m``. Fixed-effects OLS weights a bin by ``p_m (1 - p_m)`` times its mass;
    after reweighting both sides of every bin carry the bin's kernel mass, so
    ``T_B`` weights bins by their local sample (the density-weighted BATE).

    ``tau`` comes from the order-``p`` fit and the bias-corrected ``tau_bc``
    (with ``se_rb``, the interval and the p-value) from the order-``q`` fit,
    as in :func:`rdnet`.

    Parameters
    ----------
    h
        Bandwidth: a scalar, ``(h_control, h_treated)``, one value per bin, a
        ``(M, 2)`` array, or ``None`` to select per bin with :func:`rdbwnet`
        (bins whose selection fails use the median of the selected ones).
    vce
        ``"hc0"``-``"hc3"``; ``cluster`` or ``conley_cutoff`` switch to cluster
        (CR2 by default, ``cluster_vce="cr1"`` for CR1) or network Conley errors.
        Conley errors here use the pooled fit's own residuals.
    ipw
        Reweight bins to the density-weighted BATE (default). ``False`` gives
        plain fixed-effects weights ``p_m (1 - p_m)``.
    bin_slopes
        Interact the slopes with the bins' crossing features ``Z_m`` (needs
        :func:`bin_coordinates` with features or topology).
    min_obs
        Minimum distinct locations per side for a bin to enter.
    """
    if not isinstance(bins, BinCoordinates):
        raise TypeError("rdnet_fe needs bin_coordinates(...)")
    y = _check_obs(y, bins.n_obs, "y")
    p = int(p)
    q = p + 1 if q is None else int(q)
    if p < 0 or q < p:
        raise ValueError("need 0 <= p <= q")
    if not 1 < level < 100:
        raise ValueError(f"level is a percentage, e.g. 95 (got {level})")
    kfun = get_kernel(kernel)
    intensity = np.ones(bins.n_obs) if weights is None else _check_weights(weights, bins.n_obs)
    Zc = None if covs is None else _check_covs(covs, bins.n_obs)
    if cluster is not None:
        vce = "cluster"
        cluster = _check_cluster(cluster, bins.n_obs)
    if conley_cutoff is not None:
        vce = "conley"
        if conley_cutoff <= 0:
            raise ValueError("conley_cutoff must be positive")
    if vce not in ("hc0", "hc1", "hc2", "hc3", "cluster", "conley"):
        raise ValueError("vce must be 'hc0'-'hc3', 'cluster' or 'conley'")
    if vce == "cluster" and cluster is None:
        raise ValueError("vce='cluster' needs a cluster variable")
    if cluster_vce not in ("cr2", "cr1"):
        raise ValueError("cluster_vce must be 'cr2' or 'cr1'")
    M = bins.n_eval
    if h is None:
        bw = rdbwnet(y, bins, p=p, kernel=kernel, bwselect=bwselect, weights=weights, scaleregul=scaleregul,
                     vce=vce if vce in ("hc0", "hc1", "hc2", "hc3") else "hc0")
        H = bw[["h0", "h1"]].to_numpy(dtype=float)
        good = np.all(np.isfinite(H) & (H > 0), axis=1)
        if not good.any():
            raise ValueError("no bin's bandwidth could be selected; pass h=")
        H[~good] = np.median(H[good], axis=0)
    else:
        from .estimate import _bandwidth_matrix

        H = _bandwidth_matrix(h, M)
    T = bins.treated
    s = bins.s
    b = bins.bin
    kw = bins.kernel_weights if bins.kernel_weights is not None else np.ones(bins.n_obs)
    ok = (b >= 0) & np.isfinite(s) & np.isfinite(y) & (intensity > 0) & ~bins.eval_masked[np.maximum(b, 0)]
    if Zc is not None:
        ok &= np.all(np.isfinite(Zc), axis=1)
    hside = np.where(T, H[np.maximum(b, 0), 1], H[np.maximum(b, 0), 0])
    W = np.where(ok, kfun(np.abs(s) / hside) * kw * intensity, 0.0)
    ok &= W > 0
    tab = pd.DataFrame({"bin": np.arange(M), "arc": bins.eval_arc, "h0": H[:, 0], "h1": H[:, 1]})
    n1 = np.array([len(np.unique(bins.loc[ok & T & (b == m)])) for m in range(M)])
    n0 = np.array([len(np.unique(bins.loc[ok & ~T & (b == m)])) for m in range(M)])
    w1 = np.bincount(b[ok & T], weights=W[ok & T], minlength=M)
    w0 = np.bincount(b[ok & ~T], weights=W[ok & ~T], minlength=M)
    phat = np.divide(w1, w1 + w0, out=np.full(M, np.nan), where=(w1 + w0) > 0)
    use_bin = (n0 >= min_obs) & (n1 >= min_obs)
    tab["n0"], tab["n1"], tab["phat"], tab["used"] = n0, n1, phat, use_bin
    idx = np.flatnonzero(ok & use_bin[np.maximum(b, 0)])
    if use_bin.sum() == 0 or len(idx) == 0:
        raise ValueError("no bin has observations on both sides within the bandwidth")
    bi, Ti = b[idx], T[idx].astype(float)
    w = W[idx]
    if ipw:
        ph = phat[bi]
        w = w * np.where(T[idx], 1.0 / ph, 1.0 / (1.0 - ph))
    R = s[idx] / hside[idx]
    ub, binv = np.unique(bi, return_inverse=True)
    FE = np.zeros((len(idx), len(ub)))
    FE[np.arange(len(idx)), binv] = 1.0
    extra, terms = [], []
    if bin_slopes and bins.reference.shape[1]:
        Zm = bins.reference[ub]
        sd = Zm.std(axis=0)
        vary = sd > 1e-12 * np.maximum(np.abs(Zm).max(axis=0), 1.0)
        if vary.any():
            Zs = (Zm[:, vary] - Zm[:, vary].mean(axis=0)) / sd[vary]
            Zi = Zs[binv]
            extra += [R[:, None] * Zi, (Ti * R)[:, None] * Zi]
            names = [f for f, v in zip(bins.features, vary) if v]
            terms += [f"R*{f}" for f in names] + [f"T*R*{f}" for f in names]
    if bins.has_similarity:
        Si = bins.similarity[idx]
        extra += [(R * Si)[:, None], (Ti * R * Si)[:, None]]
        terms += ["R*S", "T*R*S"]
    z = stats.norm.ppf(0.5 + level / 200)
    cr2 = vce == "cluster" and cluster_vce == "cr2"
    base_vce = vce if vce in ("hc1", "hc2", "hc3") else "hc0"

    def fit(order):
        P = poly_basis(R, np.zeros_like(R), order, dim=1)[:, 1:]
        X = np.column_stack([FE, Ti, P, Ti[:, None] * P, *extra] + ([Zc[idx]] if Zc is not None else []))
        from .estimate import _identifiable

        keep = _identifiable(X, w)
        tcol = len(ub)
        if not keep[tcol]:
            raise ValueError("the treatment effect is not identified (no bin has both sides)")
        X = X[:, keep]
        c = np.zeros(X.shape[1])
        c[int(keep[:tcol].sum())] = 1.0
        f = wls(X, y[idx], w)
        tau = float(c @ f.beta)
        df = np.inf
        if vce == "conley":
            psi = f.loadings(c) * f.resid
            var = float(np.squeeze(bins.hac_quadratic(idx, psi[:, None], conley_cutoff, conley_distance)))
        elif cr2:
            Radj, df = f.cr2(c, cluster[idx])
            _, codes = np.unique(cluster[idx], return_inverse=True)
            var = float(np.sum(np.bincount(codes, weights=Radj * f.resid) ** 2))
        elif vce == "cluster":
            from .estimate import _cluster_dof

            var = variance(f.influence(c, "hc0"), "cluster", cluster=cluster[idx]) * _cluster_dof(
                len(idx), X.shape[1], cluster[idx])
        else:
            var = float(np.sum(f.influence(c, base_vce) ** 2))
        return tau, math.sqrt(max(var, 0.0)), df

    tau_p, se_p, _ = fit(p)
    tau_q, se_q, df = fit(q)
    zq = float(stats.t.ppf(stats.norm.cdf(z), df)) if np.isfinite(df) and df > 0 else z
    if se_q > 0:
        pv = 2 * (stats.t.sf(abs(tau_q / se_q), df) if np.isfinite(df) and df > 0 else stats.norm.sf(abs(tau_q / se_q)))
    else:
        pv = np.nan
    return RDNetFEResult(tau_p, se_p, tau_q, se_q, tau_q - zq * se_q, tau_q + zq * se_q, float(pv), level, p, q, vce,
                         df, tab, terms)


def similarity_sensitivity(y, bins: BinCoordinates, h_S=None, *, wbate="density", **kwargs) -> pd.DataFrame:
    """Aggregate effect over a range of similarity bandwidths ``h_S``.

    ``h_S`` defaults to half, one and two times the rule of thumb, plus
    ``np.inf`` (no similarity kernel). Each value refits :func:`rdnet` on the
    bins (``kwargs`` are passed on; with ``h`` not given, ``h_R`` is selected
    again for each ``h_S``) and aggregates with ``wbate`` (``"density"``,
    ``"gls"``, ``"place"`` or explicit weights).
    """
    from .estimate import rdnet

    if not isinstance(bins, BinCoordinates):
        raise TypeError("similarity_sensitivity needs bin_coordinates(...)")
    if h_S is None:
        base = similarity_bandwidth(bins.feature_rank)
        h_S = [0.5 * base, base, 2.0 * base, np.inf] if np.isfinite(base) else [np.inf]
    rows = []
    for hs in h_S:
        res = rdnet(y, bins.with_similarity_bandwidth(hs), **kwargs)
        rec = {"h_S": float(hs), "bins_ok": int(res._ok().sum())}
        try:
            rec.update(res.wbate(wbate).to_dict())
        except ValueError:
            rec.update(tau=np.nan, tau_bc=np.nan, se_rb=np.nan, ci_lo=np.nan, ci_hi=np.nan, pvalue=np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)
