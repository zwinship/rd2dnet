"""rdnet: boundary discontinuity designs on street networks.

Extends the location-based boundary discontinuity framework of rd2d
(Cattaneo, Titiunik and Yu) from Euclidean coordinates to network-native
coordinates and user-defined network metrics.

Typical workflow::

    import rdnet
    net = rdnet.NetworkGraph(G, metric="length")          # or travel time, transitions, ...
    bnd = rdnet.Boundary(line, treated_region=zone)
    coords = rdnet.network_coordinates(net, bnd, points)
    res = rdnet.rdnet(y, coords)
    print(res)
"""

from .bandwidth import rdbwnet
from .boundary import Boundary
from .coordinates import NetworkCoordinates, network_coordinates
from .density import DensityResult, density
from .estimate import RDNetResult, rdnet
from .kernels import get_kernel
from .metrics import EdgeCost, NetworkMetric, Precomputed, TransitionCost, node_contrast, travel_time
from .network import NetworkGraph, project_lonlat

__version__ = "0.1.0"

__all__ = [
    "Boundary",
    "DensityResult",
    "EdgeCost",
    "NetworkCoordinates",
    "NetworkGraph",
    "NetworkMetric",
    "Precomputed",
    "RDNetResult",
    "TransitionCost",
    "density",
    "get_kernel",
    "network_coordinates",
    "node_contrast",
    "project_lonlat",
    "rdbwnet",
    "rdnet",
    "travel_time",
]
