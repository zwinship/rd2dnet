# Changelog

## 0.2.0

Easier data entry, from testing the package on many kinds of applications (school zones,
minimum-wage borders, police precincts, river monitoring stations, subway fare zones,
highway exits, sensor panels).

New:

- `rd2dnet(Y, X, t, b, network)`: the data entry of rd2d (outcome, locations, treatment,
  boundary points) plus a network. `as_network()` turns graphs, igraph objects and line
  layers into a `NetworkGraph`.
- `NetworkGraph.from_lines()`: networks from GIS line layers (GeoDataFrames in any CRS,
  projected to UTM when in lon/lat), lists of lines, or segment tables with endpoint
  coordinates; joins line ends within `tol`, optionally splits lines at crossings, keeps
  attributes, and supports directed lines.
- `NetworkGraph.from_frames()`: networks from a node table and an edge table.
- `Boundary.from_points()` (treated side inferred from `t`) and `Boundary.from_region()`
  (a zone's outline, optionally restricted with `within=`).
- Distance-only local fits where the network meets the boundary at few points (rivers,
  rail, highways): the `basis` column reports `"s"` there.
- `density(unique_locations=True)` for repeated measurements.
- `RDNetResult.wbate()` defaults to equal weights over the estimated points.
- `network_coordinates(extensive=...)`.

Fixes:

- Closed-ring boundaries: offsets wrap the shorter way round, and the default evaluation
  grid no longer repeats the start point.
- The density test kept control observations at `s = 0` on the treated side.
- Splitting an edge at the boundary now divides the metric's own cost attribute (such
  as `weight`), not only `length` and `travel_time`.
- Parallel edges: snapping uses the cheapest edge, as shortest paths do.
- The bandwidth floor (the distance to the `bwcheck`-th nearest observation) is nudged
  just past that distance. The triangular kernel gives weight 0 at exactly `d = h`, so a
  bandwidth at the floor used to leave one observation fewer than `min_obs` and the point
  was reported as `insufficient`. Fixed in Python, R and Stata.
- Points with fewer than two clusters on a side get status `few clusters` instead of a
  zero standard error.
- `rdnet` reports the selector's reason (`insufficient`, `singular`) instead of
  `no bandwidth`; `h=(h0, h1)` is always read as control/treated bandwidths.
- GeoSeries points in another CRS are reprojected to the network's; networks in
  lon/lat give a warning; a missing metric attribute is reported when the network is
  built.
- Input checks: `level` must be a percentage, weights finite, cluster labels complete
  and of the right length, `treated` boolean or 0/1, `vce` valid in `rdbwnet`.
- Warnings for evaluation points beyond the boundary's network span, evaluation points
  that coincide, boundary nodes just outside `tol`, and both `cluster` and
  `conley_cutoff` given.
- `travel_time(default_kph=...)` also fills missing (NaN) speeds.
- `from_igraph` uses `length` when there is no `weight` attribute.

Fixes from cross-language and invariance testing:

- Directed networks: a location partway along a one-way edge was reached through
  either end node, as if the edge ran both ways, so `s` and kernel distances were too
  short. Transition costs: the penalty for switching onto the class of the edge a
  location sits on was never charged. With `direction="out"`/`"in"` on a directed
  network, or with `TransitionCost`, every location and evaluation point is now placed
  on a node of the network (its edge is split there), so paths respect edge direction
  and class changes. Fixed in Python, R and Stata.
- Nodes inserted where edges cross the boundary had no node attributes, so
  `node_contrast()` costs on split edges counted the contrast twice (against 0). Numeric
  node attributes are now interpolated along the edge.
- Covariates (and `pooled=True`) with different bandwidths on the two sides: each side's
  kernel weights are now divided by its `h**2`, as in rd2d, which sets the relative weight
  of the sides in the common coefficients. Estimates with one bandwidth are unchanged.
- HC2/HC3: `1 - h_ii` is clipped at 1e-8 (rd2d's value) instead of 1e-12. Observations
  with leverage 1 have zero residuals, and the tighter clip turned rounding error in
  those residuals into visible standard-error noise (about 1e-3).
- `parity/`: random designs run in Python and R (CI job `parity`) and Stata
  (`parity/run_stata.do`), compared to 1e-7 (differences at rounding level in
  automatically selected bandwidths reach a few times 1e-8).

Stata (first run in Stata 19; the test suite passes):

- `rdnetdensity` now returns its results in `r()` (they were cleared when the command ended).
- `net install` now installs the mock-city datasets, so `rdnetmock` works after installation
  (ancillary files are listed with `F` in `rd2dnet.pkg`).
- `stata/test/install_check.do` installs the package at a given commit and runs every
  help-file example from an empty directory; `stata/test/fetch.do` downloads the files
  for the test suite.

## 0.1.1

- `NetworkGraph.from_place` now uses Nominatim's top match for the query. osmnx skips
  matches without an outline, so a query for a neighbourhood mapped as a point (such as
  "Midtown Manhattan, New York") used to download the streets inside an unrelated
  building further down the list (Grand Central Terminal: 32 nodes). Matches that are
  points, or areas smaller than `min_area` (default 250,000 m²), now download the
  streets within `dist` metres (default 1,000) of the match, with a warning.
- The live OpenStreetMap test (`pytest -m online`) runs in a new GitHub Actions
  workflow (`.github/workflows/online.yml`).

## 0.1.0 (2026-10-16)

First release.

- Network boundary discontinuity estimation (`rdnet`): side-specific local polynomials
  in network coordinates, robust bias-corrected inference, uniform confidence bands,
  fuzzy designs, covariates, derivative and directional-derivative targets, HC0-HC3,
  cluster and network/Euclidean Conley standard errors.
- Aggregate effects as in rd2d: `RDNetResult.wbate()` and `RDNetResult.lbate()`, with the
  joint covariance across evaluation points stored in `RDNetResult.cov`.
- Bandwidth selection (`rdbwnet`), ported from rd2d's direct plug-in, including the
  coverage-error-optimal variants.
- Density (sorting) test (`density`), ported from rddensity, with a street-length
  normalisation.
- User-defined network metrics: edge costs, travel time, node contrasts, transition
  penalties between road classes (`TransitionCost`) and precomputed matrices.
- Scales to city-sized networks and millions of observations (see the README's Scale section).
- Tutorial notebook with donut, covariate-balance, placebo-boundary and
  difference-in-discontinuities recipes.
