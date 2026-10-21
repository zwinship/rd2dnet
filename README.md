# rd2dnet

Boundary discontinuity designs on street networks. `rd2dnet` adapts the location-based
boundary discontinuity framework of [`rd2d`](https://github.com/rdpackages/rd2d)
(Cattaneo, Titiunik & Yu) from Euclidean coordinates to **network-native coordinates**
and **user-defined network metrics**.

> **Built on rd2d.** `rd2dnet` is a direct extension of the boundary discontinuity
> methodology and software of **Matias D. Cattaneo, Rocío Titiunik and Ruiqi Rae Yu**
> ([rd2d](https://github.com/rdpackages/rd2d), part of the
> [rdpackages](https://rdpackages.github.io/) project). The estimands, side-specific local
> polynomial estimation, robust bias-corrected inference, uniform bands and the
> bandwidth selector (a port of `rdbw2d`'s direct plug-in) follow their work. The
> density test follows Cattaneo, Jansson and Ma's `rddensity`. What is new here is
> carrying that framework over to network coordinates and user-defined network
> metrics. If you use `rd2dnet`, please cite the rd2d papers listed under
> [References](#references), as well as this package. `rd2dnet` is an independent
> project and is not affiliated with or endorsed by the rd2d authors.

## Installation

```bash
pip install rd2dnet              # core
pip install "rd2dnet[osm]"       # plus OpenStreetMap download (osmnx, pyproj)
pip install "rd2dnet[plot]"      # plus plotting (matplotlib)
```

Requires Python 3.10 or newer. Releases are on [PyPI](https://pypi.org/project/rd2dnet/).

**R package.** The same methods are available in R, in [`R-package/`](R-package/)
(same functions in snake case, same numbers on the bundled data):

```r
remotes::install_github("zwinship/rd2dnet", subdir = "R-package")
```

See [`R-package/README.md`](R-package/README.md) and `vignette("rd2dnet")`.

**Stata (draft).** A Stata port lives in [`stata/`](stata/) (`rdnetset`, `rdnet`, `rdbwnet`,
`rdnetdensity`, `estat wbate`/`lbate`). It is still being tested;
[`stata/README.md`](stata/README.md) explains its status and how to run the test suite.

## Method

`rd2d` fits a local polynomial in `(X1 - b1, X2 - b2)` on each side of the boundary at
every evaluation point `b`, using kernel weights based on Euclidean distance. `rd2dnet`
replaces both pieces:

| | rd2d | rd2dnet |
|---|---|---|
| coordinates | `(X1 - b1, X2 - b2)` | `s_i`: signed network distance to the boundary (negative on the control side)<br>`o_i(b)`: signed offset along the boundary from `b` to where `i` reaches it |
| kernel distance | `‖X_i - b‖` | `d_i(b)`: shortest-path distance from `b` to `i` |
| estimate | side-specific intercept difference | same, with polynomial order `p` in `(s, o)` |
| inference | robust bias correction with `q = p + 1` | same, plus uniform bands |

Points snap to the nearest street segment on their own side. The off-network
"depth" is charged at the snapped street's own cost per metre by default (1 for
length, the street's seconds per metre for travel time), or at a fixed `snap_scale`.

## Distance is a metric you choose

Each `NetworkGraph` carries the metric that defines distance on it. Distance here means
**dissimilarity**. It does not have to be metres.

```python
import rd2dnet
from rd2dnet import EdgeCost, TransitionCost, Precomputed, travel_time, node_contrast

rd2dnet.NetworkGraph(G, "length")                      # metres along streets
rd2dnet.NetworkGraph(G, travel_time())                 # seconds, from speed_kph
rd2dnet.NetworkGraph(G, node_contrast("far", 50.0))    # length + 50 * |land-use change|
rd2dnet.NetworkGraph(G, lambda u, v, d: ...)           # any non-negative edge cost
rd2dnet.NetworkGraph(G, Precomputed(D, list(G.nodes))) # your own node-by-node matrix (rows in that order)

# Switching road class costs extra. Ramps are transparent connectors, so
# highway -> ramp -> highway is close while highway -> ramp -> local street is far.
rd2dnet.NetworkGraph(G, TransitionCost(
    "length", class_attr="highway",
    penalty={("motorway", "residential"): 2000.0},   # or default_penalty=..., or a callable
    connectors={"motorway_link", "trunk_link"},
))
```

The coordinates `s` and `o`, the kernel distances and the bandwidths are all in the
metric's units. Treatment side is always geographic.

## Quick start

```python
import rd2dnet
from rd2dnet.datasets import load_mock_midtown   # bundled, works offline

city = load_mock_midtown()
net = city.network("length")
coords = rd2dnet.network_coordinates(net, city.boundary, city.obs[["x", "y"]], n_eval=21)

res = rd2dnet.rdnet(city.obs["outcome"], coords,
                  weights=city.obs["far"],          # intensity / floor-area scaling
                  cluster=city.obs["block"])        # or conley_cutoff=150.0
print(res)          # pointwise BATEC table
res.plot()

res.bandwidths                                         # the bandwidths rdnet selected
print(rd2dnet.density(coords))                         # sorting test (density of the distance to the boundary)
print(rd2dnet.manipulation_test(coords))               # sorting test (balance near the boundary)
take_up = city.obs["treated"]                          # stand-in: your 0/1 (or share) take-up
rd2dnet.rdnet(city.obs["outcome"], coords, fuzzy=take_up)   # fuzzy design
```

### Getting your data in

**Same data entry as rd2d.** `rd2d(Y, X, t, b)` takes an outcome, locations, a
treatment indicator and points along the boundary. `rd2dnet` takes the same four, plus
the network:

```python
res = rd2dnet.rd2dnet(Y, X, t, b, network)   # X: (n, 2), t: 0/1, b: (J, 2) ordered along the boundary
res.wbate()                                  # equal weights over the estimated points
rd2dnet.density(res.coords)                  # coordinates are kept on the result
```

The boundary is the line through `b`, and the treated side is the side where the
`t = 1` observations lie. Pass `boundary=rd2dnet.Boundary(...)` when the line extends
beyond the first and last points, or when you have the treated zone as a polygon.

**Networks, in whatever form you have them:**

| you have | use |
|---|---|
| a GIS layer of roads, rivers, rail or pipes (GeoDataFrame of LineStrings, any CRS) | `NetworkGraph.from_lines(gdf)`: line ends that meet become nodes, other columns become edge attributes, lon/lat is projected to UTM, and `split_at_intersections=True` handles layers that are not split at junctions |
| a table of segments `x1, y1, x2, y2` (no node ids) | `NetworkGraph.from_lines(df)` |
| a node table and an edge table | `NetworkGraph.from_frames(nodes, edges, node_id="id", source="from", target="to")` |
| a networkx graph (nodes with `x`, `y`) | `NetworkGraph(G)` |
| an OSMnx graph, or just a place name | `NetworkGraph.from_osmnx(G)`, `NetworkGraph.from_place("...")` |
| an igraph graph or a weighted adjacency matrix | `NetworkGraph.from_igraph(g)`, `NetworkGraph.from_adjacency(A, coords)` |
| a travel-time matrix from another router | `NetworkGraph(G, Precomputed(D, nodes))`, or `Precomputed.from_od(table, "origin", "destination", "cost")` for a long table (the boundary must pass through nodes) |

Directed networks (one-way streets, rivers drawn downstream) take `directed=True` in
`from_lines`/`from_frames` and `direction="out"` (distances run from the boundary to
each location) or `"in"` (from each location to the boundary).

**Boundaries:** `Boundary(line, treated_side="left")`, `Boundary(line, treated_region=zone)`,
`Boundary.from_region(zone)` for a zone's outline (a school attendance zone, a city, a
congestion zone; `within=` keeps only the stretch shared with a neighbour), and
`Boundary.from_points(b, X, t)` as in rd2d. Closed outlines are handled: offsets go the
shorter way round. If the boundary follows a street, set `tol` (in the network's units,
for example 0.5 m) so the street's nodes count as lying on it.

**Coordinates:** points can be an array, a DataFrame with `x`/`y`, or a GeoSeries in
any CRS (it is reprojected to the network's). Use `lonlat=True` for plain lon/lat arrays.

**Sparse crossings.** Where the network meets the boundary at only a few points (rivers,
rail lines, highways crossing a border), there is no variation along the boundary near
each crossing. Those evaluation points are fitted in `s` alone, and the `basis` column of
the results reads `"s"`. Put one evaluation point at each crossing.

**Repeated measurements** (panels, repeat sales, sensor readings): each location is
snapped once. Use `density(coords, unique_locations=True)` for the sorting test. For
standard errors, clustering by unit can understate the variance of spatially correlated
shocks; `conley_cutoff=...` is the safer choice.

The tutorial notebook [`examples/tutorial.ipynb`](examples/tutorial.ipynb) walks through
every feature on the mock city (metrics, coordinates, bandwidths, aggregate effects,
covariates, fuzzy designs, derivatives, cluster and Conley errors, the density test) and
the recipes for donut designs, covariate balance, placebo boundaries and
difference-in-discontinuities.

For real data, use `rd2dnet.NetworkGraph.from_place("Manhattan, New York")`, which needs
`pip install "rd2dnet[osm]"`. The query's top geocoding match is used: its outline when it
is an area of at least `min_area` (250,000 m² by default), otherwise (a neighbourhood
mapped as a point, a single building) the streets within `dist` metres of it (1,000 by
default). If the download fails it raises the error; `offline_fallback=True` returns the
bundled mock network instead (for demos only).

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
  `rdnet` selects with `scaleregul=3`, as rd2d does for estimation; the bandwidths it
  used are in `res.bandwidths`.
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
- **Standard errors** (`vce`): `hc0` to `hc3` (default `hc3`; rd2d uses `hc1`, but on
  networks local samples often contain high-leverage points, and HC3 kept pointwise
  coverage at 0.945-0.951 across our Monte Carlo designs where HC1 gave 0.92-0.94),
  `cluster` (CR1 with rd2d's
  `(n - 1)/(n - k) * g/(g - 1)` factor), or `conley` (Bartlett
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
- **Aggregate effects**: `res.wbate(weights)` gives the weighted boundary average
  treatment effect `w' tau` (weights normalised to sum to one) with
  `se = sqrt(w' Sigma w)`, and `res.lbate()` the largest effect `max_j tau(b_j)` with
  the interval `[max_j(tau_bc - c se), max_j(tau_bc + c se)]`, `c` the uniform critical
  value. Both follow rd2d's `summary(WBATE=, LBATE=)`; `res.summary(wbate=w, lbate=True)`
  prints them. `res.cov` always holds the J x J covariance `Sigma` of the bias-corrected
  estimates across evaluation points (NaN for points not estimated; its diagonal is
  `se_rb**2`), so aggregates account for the correlation between nearby points.
  Masked or unestimated points need weight 0.
- **Sorting tests.** Sorting is about where units locate, so both tests work in the plane,
  not on the network.
  - `density(coords)`: the local polynomial density estimator of Cattaneo, Jansson & Ma
    (as in `rddensity`), with jackknife standard errors and bias-corrected inference, on
    the signed straight-line distance from each location to the boundary
    (`normalize="euclidean"`, the default). In Monte Carlo runs it held its size
    (0.03-0.06 at 5%) on straight, round, wavy and irregular boundaries, and had the most
    power against sorting in one direction. The alternatives `normalize="none"` (network
    distance `s`) and `"network"` (observations per metre of street) are oversized at
    curved boundaries (0.08-0.33 and 0.34-1.00 without any sorting).
  - `manipulation_test(coords)`: the test of Bugni, Crippa & Restrepo (2026): among the
    `q` locations nearest the boundary (default 200), the treated share should be 1/2
    in every stretch of the boundary. It held its size on every design and is the only
    test with power against sorting that offsets along the boundary; report a few `q`,
    e.g. `q=[50, 100, 200]`. Python, R and Stata (`rdnetmanip`) give identical answers.
  - Neither is valid when the observations lie on the network itself (points on
    streets); a warning says so.

## Validation (Monte Carlo on the bundled mock city)

| check | result |
|---|---|
| robust 95% CI, pointwise coverage | 0.925 (n = 4,000), 0.940 (n = 16,000) |
| uniform 95% band coverage | 0.873 (n = 4,000), 0.940 (n = 16,000) |
| first stage, fuzzy design | unbiased; HC standard error matches Monte Carlo sd within 5% |
| kink (`deriv=(1, 0)`) | exact without noise; unbiased, standard error matches Monte Carlo sd, coverage 0.94–0.955 |
| WBATE (equal weights), 95% CI coverage | 0.910 (n = 4,000), 0.947 (n = 16,000); ignoring the cross-point covariance gives 0.717 / 0.774 |
| WBATE standard error vs Monte Carlo sd | 0.128 vs 0.129 (n = 16,000); median 0.237 vs robust sd 0.262 (n = 4,000) |
| LBATE 95% CI coverage | 0.982 (n = 4,000), 0.997 (n = 16,000): conservative by construction |
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
- At n = 4,000 the WBATE error is heavy-tailed (a few replications with tiny local
  samples, max error about 8 with matching large standard errors). The plain Monte Carlo
  sd (0.42) then exceeds the mean SE (0.25); on the robust scale the SE is about 10%
  too small, and coverage is 0.91.
- **Clustering with few clusters per side within the bandwidth under-covers badly.** In a
  Monte Carlo with 300 m block shocks and clustering by block (5-6 clusters per side in the
  bandwidth), pointwise coverage was 0.70-0.74 and uniform-band coverage about 0.2. Use
  Conley errors for spatially correlated shocks. For panels, clustering by unit needs many
  sites per side (coverage 0.88 with about 48 sites per side, 0.94 with about 190).
- Cluster-robust errors use the local fit's own residuals, as rd2d does. When clusters
  are small relative to the bandwidth they are too small: clustering by city block
  (n = 4,000, selected bandwidths, 400 replications) gives 95% coverage of 0.78 with
  independent errors (HC1: 0.93), and neither HC1 nor cluster errors capture block-level
  shocks (coverage 0.69 for both; their SEs are about half the Monte Carlo sd). The
  local polynomial absorbs within-cluster variation, the same effect that motivates
  the Conley default's wider auxiliary fit. For spatially correlated shocks use Conley.
- Observations at exactly the same location are at network distance 0 in the Conley
  kernel (as with `conley_distance="euclidean"`). Before the sparse kernel they were
  treated as `2 x snapping depth` apart, which understated the correlation between
  repeated readings at one site.
- The LBATE interval inherits the uniform band's joint coverage and is conservative
  for the maximum (0.98-1.00). Its point estimate `max_j tau_bc` is biased upward in
  small samples (+0.42 at n = 4,000, +0.05 at n = 16,000), as any maximum of noisy
  estimates is.

## Scale

Timings in seconds on one desktop (AMD Ryzen 7 3700X, 8 cores, 32 GB), 25 evaluation
points, bandwidths selected by `rdbwnet`, Conley cutoff 150 m (city), 200 m (suburb) or
1 km (river). The networks are synthetic (`benchmarks/synthetic.py`): a jittered grid
city with arterials, a motorway with ramps and a park on the boundary; an irregular
Delaunay suburb with cul-de-sacs and a curved boundary that cuts edges; and a dendritic
river basin (a tree with edges directed downstream). Observations are uniform
parcels, GPS-like pings near the streets, or a sensor panel (many readings at a few
thousand fixed sites).

| network (nodes) | observations | coords | `rdbwnet` | `rdnet` HC1 | cluster | Conley | `density` |
|---|---|---|---|---|---|---|---|
| city (50k) | 100k parcels | 2.1 | 0.9 | 0.2 | 0.2 | 1.1 | 0.7 |
| city (50k) | 1M parcels | 4.5 | 7.1 | 2.9 | 3.3 | 11.8 | 3.9 |
| city (50k) | 1M GPS pings | 4.0 | 6.5 | 2.7 | 2.7 | 12.9 | 4.1 |
| city (50k) | 5M GPS pings | 18.1 | 35.9 | 12.7 | 14.0 | 280 | 26.2 |
| city (50k) | 1M readings, 5k sensors | 1.0 | 5.4 | 1.2 | 1.4 | 4.2 | 3.0 |
| metro (200k) | 1M parcels | 7.3 | 9.0 | 2.6 | 3.2 | 8.7 | 6.8 |
| suburb (50k) | 200k parcels | 5.2 | 1.8 | 0.6 | 0.6 | 4.0 | 1.3 |
| river (50k) | 200k GPS pings | 2.1 | 1.3 | 0.5 | 0.5 | 2.8 | 0.5 |
| river (50k) | 10M readings, 10k sensors | 2.4 | 62.0 | 10.3 | 10.9 | 24.9 | 8.9 |

Building a 50k-node `TransitionCost` engine takes about 0.6 s (2.4 s for 200k nodes).

How the cost grows:
- **Coordinates**: snapping and boundary distances are computed once per distinct
  location (x, y, side), so a panel costs as much as its sites. Snapping and the
  geometry tests run on threads.
- **`rdbwnet` and `rdnet`**: linear in n per evaluation point. rd2d's bias pilot uses
  the median distance, so those fits see about half of each side's sample; this is
  the largest cost at 10M observations.
- **Conley**: proportional to the number of location pairs closer than the cutoff.
  That grows with the square of local density: 5x the observations in the same city
  is about 25x the pairs (5M pings, roughly 10^9 pairs, takes under 5 minutes). Node
  pairs come from Dijkstra searches bounded at the cutoff on the local subgraph,
  and the kernel is assembled in sparse row blocks of bounded memory. With very dense
  data, a smaller cutoff or clustered standard errors are far cheaper.
- **Density test**: O(n log n) plus a pass over the edges.

`RD2DNET_NUM_THREADS` caps the number of threads (default: up to 16 for geometry,
8 for the Conley kernel, 4 for evaluation points in `rdbwnet`; `1` disables
threading). Threaded results are identical to serial ones.

Reproduce with `python benchmarks/profile_scale.py --net grid_city --n 1000000
--design gps` (add `--profile STAGE` for a cProfile breakdown).

## Development

```bash
git clone https://github.com/zwinship/rd2dnet.git && cd rd2dnet
python -m venv .venv
.venv\Scripts\activate          # Windows   (source .venv/bin/activate on macOS/Linux)
pip install -e ".[dev]"
pytest                          # offline, includes doctests; `-m online` runs the live OSM test
pytest -m notebook              # executes examples/tutorial.ipynb
```

On Windows, check that `python` is your real interpreter and not one bundled with
another app (such as Inkscape). `py -0p` lists the installed interpreters.

GitHub Actions (`.github/workflows/ci.yml`) runs the tests on Python 3.10-3.13 (Linux) and
3.12 (Windows), executes the tutorial, and builds and checks the distribution.
`examples/make_tutorial.py` regenerates the notebook.

The R package is developed in `R-package/` (C++ core via Rcpp, roxygen2 documentation,
testthat tests including parity checks against this package's numbers):

```bash
cd R-package
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'
R CMD build . && R CMD check --as-cran rd2dnet_*.tar.gz
```

`.github/workflows/R-CMD-check.yml` runs `R CMD check --as-cran` on Linux (R release,
devel and oldrel), macOS and Windows. `R-package/data-raw/` regenerates the bundled data
from this package's files.

## References

`rd2dnet` builds directly on the following work. Please cite it.

**rd2d (Cattaneo, Titiunik and Yu)**, see https://rdpackages.github.io/rd2d/

- Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). rd2d: Causal Inference in Boundary
  Discontinuity Designs. Working paper. arXiv:2505.07989.
- Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). Boundary Discontinuity Designs:
  Theory and Practice. *Advances in Economics and Econometrics: Thirteenth World
  Congress*, Cambridge University Press, Vol. 1, Ch. 2, to appear.
- Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). Estimation and Inference in Boundary
  Discontinuity Designs: Location-Based Methods. *Journal of the American Statistical
  Association*, revise and resubmit.
- Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). Estimation and Inference in Boundary
  Discontinuity Designs: Distance-Based Methods. *Journal of Econometrics* 256(A): 106266.
- Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). Estimation and Inference in Boundary
  Discontinuity Designs: Pooling-Based Methods. Working paper.

**Sorting tests**

- Bugni, F. A., Crippa, F. and Restrepo, D. (2026). Manipulation Testing in Boundary
  Discontinuity Designs. arXiv:2609.00350.

- Cattaneo, M. D., Jansson, M. and Ma, X. (2020). Simple Local Polynomial Density
  Estimators. *Journal of the American Statistical Association* 115(531): 1449–1455.
- Cattaneo, M. D., Jansson, M. and Ma, X. (2018). Manipulation Testing Based on Density
  Discontinuity. *Stata Journal* 18(1): 234–261.

**Robust bias-corrected inference and spatial HAC**

- Calonico, S., Cattaneo, M. D. and Titiunik, R. (2014). Robust Nonparametric Confidence
  Intervals for Regression-Discontinuity Designs. *Econometrica* 82(6): 2295–2326.
- Conley, T. G. (1999). GMM Estimation with Cross Sectional Dependence. *Journal of
  Econometrics* 92(1): 1–45.

## License

GPL-3.0-only, the same license as rd2d and rddensity, from which parts of the
implementation are ported. See [LICENSE](LICENSE).
