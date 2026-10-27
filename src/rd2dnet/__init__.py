"""rd2dnet: boundary discontinuity designs on street networks.

A direct extension of the location-based boundary discontinuity framework and
software of rd2d (Cattaneo, Titiunik and Yu; https://rdpackages.github.io/rd2d/)
from Euclidean coordinates to network-native coordinates and user-defined
network metrics. Please cite the rd2d papers when using this package.

Licensed under GPL-3.0-only, like rd2d and rddensity, from which parts are ported.

With the data entry of rd2d (outcome, locations, treatment, boundary points)::

    res = rd2dnet.rd2dnet(Y, X, t, b, network)

Step by step::

    import rd2dnet
    net = rd2dnet.NetworkGraph(G, metric="length")          # or travel time, transitions, ...
    bnd = rd2dnet.Boundary(line, treated_region=zone)
    coords = rd2dnet.network_coordinates(net, bnd, points)
    res = rd2dnet.rdnet(y, coords)
    print(res)
"""

from .bandwidth import rdbwnet
from .bins import BinCoordinates, RDNetFEResult, bin_coordinates, rdnet_fe, similarity_bandwidth, similarity_sensitivity
from .boundary import Boundary
from .coordinates import NetworkCoordinates, network_coordinates
from .density import DensityResult, density
from .estimate import RDNetResult, rdnet
from .interface import as_network, rd2dnet
from .manipulation import ManipulationTestResult, manipulation_test
from .kernels import get_kernel
from .metrics import EdgeCost, NetworkMetric, Precomputed, TransitionCost, node_contrast, travel_time
from .network import NetworkGraph, project_lonlat

__version__ = "0.2.0"

__all__ = [
    "BinCoordinates",
    "Boundary",
    "DensityResult",
    "EdgeCost",
    "ManipulationTestResult",
    "NetworkCoordinates",
    "NetworkGraph",
    "NetworkMetric",
    "Precomputed",
    "RDNetFEResult",
    "RDNetResult",
    "TransitionCost",
    "as_network",
    "bin_coordinates",
    "density",
    "get_kernel",
    "manipulation_test",
    "network_coordinates",
    "node_contrast",
    "project_lonlat",
    "rd2dnet",
    "rdbwnet",
    "rdnet",
    "rdnet_fe",
    "similarity_bandwidth",
    "similarity_sensitivity",
    "travel_time",
]
