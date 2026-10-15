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
print(rdnet.density(coords))                         # sorting test per metre of street
rdnet.rdnet(city.obs["outcome"], coords, fuzzy=take_up)   # fuzzy design (take_up = your 0/1 or share)
```

For real data, use `rdnet.NetworkGraph.from_place("Manhattan, New York")`, which needs
`pip install -e .[osm]`. It falls back to the bundled mock network with a warning when
offline. You can also pass your own graph:
`NetworkGraph(G)`, `NetworkGraph.from_adjacency(A, coords)` or `NetworkGraph.from_igraph(g)`.

### Features

- `Boundary(line, treated_region=... | treated_side="left"/"right", mask=[park_polygon, (a, b)])`:
  masked stretches with no mass, such as Central Park, are skipped by both the
  estimator and the density test.
- Edges that cross the boundary are split automatically, so every cross-side path
  passes through a boundary node. Points on a street segment that runs along the
  boundary count as already being on it.
- **Bandwidths** (`rdbwnet`): a port of rd2d's two-stage direct plug-in. It uses a
  rule-of-thumb pilot, a separate pilot bandwidth for the bias, regularisation
  (`scaleregul`) and a higher-order bias term (`scalebiascrct`). `bwselect` takes
  `mserd`, `msetwo` (one bandwidth per side), `imserd` or `imsetwo`, plus the
  coverage-error-optimal `cerrd`, `certwo`, `icerrd` and `icertwo`. You can also
  pass `h` as a scalar, `(h_control, h_treated)`, or one value per evaluation point.
- **Covariates** (`covs=`): the coefficients are common to both sides, as in rd2d's
  `covs_eff`, and estimated in one joint local regression whose influence function
  accounts for estimating them. The point estimate matches rd2d's plug-in exactly.
  Covariates that are constant in a local sample are dropped there.
- **Derivative targets**: `deriv=(1, 0)` gives the jump in the slope *across* the
  boundary (a kink) and `deriv=(0, 1)` the jump in the slope *along* it.
  `tangvec=(t_s, t_o)` gives any directional derivative. Bandwidth selection adapts,
  as in rd2d.
- **Fuzzy designs**: `rdnet(y, coords, fuzzy=D)` estimates the ratio of the outcome
  jump to the take-up jump at each boundary point. It reports the first stage `fs`
  and `fs_se`, and selects the bandwidth from the linearised outcome, as rd2d does.
- **Standard errors** (`vce`): `hc0` to `hc3`, `cluster`, or `conley` (Bartlett
  spatial HAC over **network** or Euclidean distance). Conley is computed block by
  block, so memory stays bounded, and it also drives the uniform band's joint
  covariance across boundary points. Conley residuals come from a smoother
  auxiliary fit (`conley_resid_scale`, default 3 × h, and `conley_resid_order`,
  default q). Residuals from the estimation fit itself absorb spatially correlated
  shocks and make the HAC variance far too small. The default was chosen by
  simulation across five error designs (see Validation).
- `weights`: parcel intensity (floor-area ratio (FAR), building height) multiplies the
  kernel weights.
- `pooled=True` gives the common-slope model `y ~ 1 + T + s + o`. The default is
  side-specific fits, as in rd2d.
- **Density test** (`density`): the local polynomial density estimator of
  Cattaneo, Jansson & Ma (as in `rddensity`), with jackknife standard errors and
  bias-corrected inference. The default `normalize="network"` measures position in
  *cumulative street length* from the boundary, so the test asks whether there are
  more parcels per metre of street on one side. On a street grid the raw density
  of `s` is jagged, and street supply differs across sides (for example a park
  frontage), so the unnormalised test (`normalize="none"`) can flag geometry as
  sorting. The default `bwselect="common"` uses one bandwidth for both sides.

## Validation (Monte Carlo on the bundled mock city)

| check | result |
|---|---|
| robust 95% CI, pointwise coverage | 0.925 (n = 4,000), 0.940 (n = 16,000) |
| uniform 95% band coverage | 0.873 (n = 4,000), 0.940 (n = 16,000) |
| first stage, fuzzy design | unbiased; HC standard error matches Monte Carlo sd within 5% |
| kink (`deriv=(1, 0)`) | exact without noise; unbiased, standard error matches Monte Carlo sd, coverage 0.94–0.955 |
| density test size at 5%, `network` / `none` (no-park city) | 0.030 / 0.065 |
| density test size at 5%, `network` / `none` (park city, masked) | 0.025 / 0.085 |

Spatially correlated errors (n = 8,000, h = 500 m, Conley cutoff = 2 × correlation scale;
each cell is pointwise / uniform-band coverage):

| error design | HC1 | Conley (default) |
|---|---|---|
| independent | 0.917 / 0.842 | 0.962 / 0.958 |
| smooth field, 100 m | 0.825 / 0.617 | 0.927 / 0.808 |
| smooth field, 200 m | 0.895 / 0.825 | 0.955 / 0.933 |
| block shocks, 150 m | 0.900 / 0.800 | 0.958 / 0.950 |
| block shocks, 300 m | 0.885 / 0.858 | 0.958 / 0.883 |

**Known limitations.**
- Coverage is below nominal when local samples are small (about 150 observations
  per side) and improves with n. rd2d notes the same.
- Conley standard errors are conservative (up to about 1.5 × the true sd) when the
  spatial correlation scale is close to the bandwidth. A residual-based HAC cannot
  fully separate correlated noise from the regression function at that scale, and
  the default errs toward wider intervals.
- The uniform band is weakest for short-range smooth correlation (0.81). Pointwise
  coverage there is 0.93.

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

See [ROADMAP.md](ROADMAP.md) for planned next steps.

## References

- Cattaneo, Titiunik & Yu. *rd2d: Causal Inference in Boundary Discontinuity Designs.* arXiv:2505.07989.
- McCrary (2008). *Manipulation of the running variable in the regression discontinuity design.*
- Conley (1999). *GMM estimation with cross sectional dependence.*
