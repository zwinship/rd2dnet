"""rd2d-style entry point: ``rd2dnet(Y, X, t, b, network)``."""

from __future__ import annotations

import networkx as nx
import numpy as np

from .boundary import Boundary
from .coordinates import _points_array, network_coordinates
from .estimate import RDNetResult, rdnet
from .network import NetworkGraph, project_lonlat


def as_network(network, metric="length", **kwargs) -> NetworkGraph:
    """Turn common network inputs into a :class:`NetworkGraph`.

    Accepts a :class:`NetworkGraph` (returned as is), a networkx graph whose
    nodes have ``x``/``y`` attributes, a python-igraph graph, or line
    geometries (a GeoDataFrame/GeoSeries of LineStrings, or a DataFrame of
    segments with ``x1, y1, x2, y2``), passed to :meth:`NetworkGraph.from_lines`.
    """
    if isinstance(network, NetworkGraph):
        return network
    if isinstance(network, nx.Graph):
        return NetworkGraph(network, metric, **kwargs)
    if type(network).__module__.startswith("igraph"):
        return NetworkGraph.from_igraph(network, metric, **kwargs)
    return NetworkGraph.from_lines(network, metric, **kwargs)


def rd2dnet(
    Y,
    X,
    t,
    b,
    network,
    *,
    boundary: Boundary | None = None,
    treated_side: str | None = None,
    metric="length",
    snap: str = "edge",
    lonlat: bool = False,
    **kwargs,
) -> RDNetResult:
    """Boundary discontinuity estimates with the data entry of rd2d.

    ``rd2dnet(Y, X, t, b, network)`` mirrors ``rd2d(Y, X, t, b)``: the same
    outcome, locations, treatment and boundary points, plus the network on
    which distance is measured.

    Parameters
    ----------
    Y
        Outcome, length ``n``.
    X
        Locations, ``(n, 2)`` (array, DataFrame with ``x``/``y``, or a
        GeoSeries of points), in the network's projected coordinates.
    t
        Treatment indicator (0/1 or boolean), length ``n``.
    b
        Evaluation points along the boundary, ``(J, 2)``, in order along it.
    network
        A :class:`NetworkGraph`, a networkx graph with ``x``/``y`` node
        attributes, an igraph graph, or line geometries (see
        :func:`as_network`).
    boundary
        The boundary line. By default it is the line through ``b``, with the
        treated side taken from where the treated observations lie. Pass a
        :class:`Boundary` when the line extends beyond the first and last
        evaluation points or when a treated region is known.
    treated_side
        ``"left"`` or ``"right"`` of the line through ``b`` (instead of
        inferring it from ``X`` and ``t``).
    metric
        Distance on the network when ``network`` is not already a
        :class:`NetworkGraph`.
    snap
        ``"edge"`` or ``"node"``, passed to :func:`network_coordinates`.
    lonlat
        ``X`` and ``b`` are longitude/latitude; they are projected into the
        network's CRS.
    **kwargs
        Passed to :func:`rdnet`: ``h``, ``p``, ``kernel``, ``bwselect``, ``vce``,
        ``cluster``, ``fuzzy``, ``covs``, ``level``, ...

    Returns
    -------
    RDNetResult
        Estimates at each point of ``b``; ``result.coords`` holds the network
        coordinates for :func:`density` and :func:`rdbwnet`.

    Examples
    --------
    >>> import rd2dnet
    >>> from rd2dnet.datasets import load_mock_midtown
    >>> city = load_mock_midtown()
    >>> b = city.boundary.interpolate(city.boundary.grid(9))
    >>> obs = city.obs
    >>> res = rd2dnet.rd2dnet(obs["outcome"], obs[["x", "y"]], obs["treated"], b, city.network(),
    ...                       boundary=city.boundary, seed=0)
    """
    net = as_network(network, metric)
    xy = _points_array(X, net.crs)
    b = np.asarray(b, dtype=float).reshape(-1, 2)
    if lonlat:
        if net.crs is None:
            raise ValueError("network has no CRS to project lon/lat into")
        xy = project_lonlat(xy[:, 0], xy[:, 1], net.crs)
        b = project_lonlat(b[:, 0], b[:, 1], net.crs)
    t = np.asarray(t).astype(float).ravel()
    if np.isnan(t).any() or not np.isin(t, (0, 1)).all():
        raise ValueError("t must be a 0/1 (or boolean) treatment indicator without missing values")
    t = t.astype(bool)
    if len(t) != len(xy):
        raise ValueError("X and t must have the same number of observations")
    if boundary is None:
        boundary = Boundary.from_points(b, xy, t, treated_side=treated_side)
    coords = network_coordinates(net, boundary, xy, treated=t, eval_points=b, snap=snap)
    return rdnet(Y, coords, **kwargs)
