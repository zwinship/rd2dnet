# rdnet

Boundary discontinuity designs on street networks. `rdnet` adapts the location-based
boundary discontinuity framework of [`rd2d`](https://github.com/rdpackages/rd2d)
(Cattaneo, Titiunik & Yu) from Euclidean coordinates to **network-native coordinates**
and **user-defined network metrics**. The motivating application is a congestion pricing
boundary like Manhattan's 60th Street.

## Method

`rd2d` fits a local polynomial in `(X1 - b1, X2 - b2)` on each side of the boundary at
every evaluation point `b`, using kernel weights based on Euclidean distance. `rdnet`
replaces both pieces:

| | rd2d | rdnet |
|---|---|---|
| coordinates | `(X1 - b1, X2 - b2)` | `s_i`: signed network distance to the boundary (negative on the control side)<br>`o_i(b)`: signed offset along the boundary from `b` to where `i` reaches it |
| kernel distance | `‖X_i - b‖` | `d_i(b)`: shortest-path distance from `b` to `i` |
| estimate | side-specific intercept difference | same, with polynomial order `p` in `(s, o)` |
| inference | robust bias correction with `q = p + 1` | same, plus uniform bands |

Points snap to the nearest street segment on their own side. The off-network
"depth" is charged through `snap_scale`.

## Distance is a metric you choose

Each `NetworkGraph` carries the metric that defines distance on it. Distance here means
**dissimilarity**. It does not have to be metres.

```python
import rdnet
from rdnet import EdgeCost, TransitionCost, Precomputed, travel_time, node_contrast

rdnet.NetworkGraph(G, "length")                      # metres along streets
rdnet.NetworkGraph(G, travel_time())                 # seconds, from speed_kph
rdnet.NetworkGraph(G, node_contrast("far", 50.0))    # length + 50 * |land-use change|
rdnet.NetworkGraph(G, lambda u, v, d: ...)           # any non-negative edge cost
rdnet.NetworkGraph(G, Precomputed(D, nodes))         # your own node-by-node matrix

# Switching road class costs extra. Ramps are transparent connectors, so
# highway -> ramp -> highway is close while highway -> ramp -> local street is far.
rdnet.NetworkGraph(G, TransitionCost(
    "length", class_attr="highway",
    penalty={("motorway", "residential"): 2000.0},   # or default_penalty=..., or a callable
    connectors={"motorway_link", "trunk_link"},
))
```

The coordinates `s` and `o`, the kernel distances and the bandwidths are all in the
metric's units. Treatment side is always geographic.

## Quick start

```python
import rdnet
from rdnet.datasets import load_mock_midtown   # bundled, works offline

city = load_mock_midtown()
net = city.network("length")
coords = rdnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]], n_eval=21)

res = rdnet.rdnet(city.obs["outcome"], coords,
                  weights=city.obs["far"],          # intensity / floor-area scaling
                  cluster=city.obs["block"])        # or conley_cutoff=150.0
print(res)          # pointwise BATEC table
res.plot()

rdnet.rdbwnet(city.obs["outcome"], coords)          # bandwidths
print(rdnet.density(coords, normalize="network"))   # sorting test per unit street length
```

For real data, use `rdnet.NetworkGraph.from_place("Manhattan, New York")`, which needs
`pip install -e .[osm]`. It falls back to the bundled mock network with a warning when
offline. You can also pass your own graph:
`NetworkGraph(G)`, `NetworkGraph.from_adjacency(A, coords)` or `NetworkGraph.from_igraph(g)`.

### Features

- `Boundary(line, treated_region=... | treated_side="left"/"right", mask=[park_polygon, (a, b)])`:
  masked stretches with no mass, such as Central Park, are skipped.
- Edges that cross the boundary are split automatically, so every cross-side path
  passes through a boundary node.
- `vce`: `hc0` to `hc3`, `cluster`, or `conley` (Bartlett spatial HAC over **network**
  or Euclidean distance).
- `weights`: parcel intensity (floor-area ratio (FAR), building height) multiplies the
  kernel weights.
- `pooled=True` gives the common-slope model `y ~ 1 + T + s + o`. The default is
  side-specific fits, as in rd2d.
- `rdbwnet`: direct plug-in MSE (`mserd`) or common IMSE (`imserd`) bandwidths.
- `density`: a McCrary-style test on `s`. With `normalize="network"` it is per unit of
  street length, so the test is not confounded by uneven street supply.

## Development

```bash
git clone https://github.com/zwinship/rdnet.git && cd rdnet
python -m venv .venv
.venv\Scripts\activate          # Windows   (source .venv/bin/activate on macOS/Linux)
pip install -e ".[dev]"
pytest                          # offline; `pytest -m online` runs the live OSM test
```

On Windows, check that `python` is your real interpreter and not one bundled with
another app (such as Inkscape). `py -0p` lists the installed interpreters.

## References

- Cattaneo, Titiunik & Yu. *rd2d: Causal Inference in Boundary Discontinuity Designs.* arXiv:2505.07989.
- McCrary (2008). *Manipulation of the running variable in the regression discontinuity design.*
- Conley (1999). *GMM estimation with cross sectional dependence.*
