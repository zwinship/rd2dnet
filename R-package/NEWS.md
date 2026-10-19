# rd2dnet 0.2.0

## Easier data and network entry

* `rd2dnet(Y, X, t, b, network)`: the data entry of rd2d's `rd2d(Y, X, t, b)`, with the
  boundary drawn through the evaluation points `b` and the treated side inferred from `t`.
* `network_graph_from_lines()`: networks from GIS line layers (sf `LINESTRING`s, lists of
  vertex matrices, or segment tables with `x1, y1, x2, y2`). Ends within `tol` are joined,
  attributes are kept for metrics, layers in longitude/latitude are projected to their UTM
  zone, and `split_at_intersections = TRUE` nodes layers that are not split at junctions.
* `boundary_from_points()` and `boundary_from_region()` (the outline of a treated zone,
  optionally restricted to the stretch `within` a polygon).
* `network_graph()` computes `length` when it is missing, and errors early when the
  metric's column is missing. sf points in another CRS are transformed to the network's,
  data frames with two unnamed coordinate columns are accepted, and networks in
  longitude/latitude give a warning.
* `network_graph_from_igraph()` falls back to `length` when there is no `weight`, accepts
  `directed`, and reads sfnetworks-style vertex geometries.

## Estimation

* Where the network meets the boundary at only a few points (rivers, rail lines, highways
  crossing a border), the local fit uses a polynomial in `s` alone instead of failing as
  singular; the new `basis` column of `rdnet()` and `rdbwnet()` results says which fit was used.
* Closed boundaries (zone outlines) are treated as rings: offsets wrap around and the
  default evaluation grid does not repeat the starting point.
* `"few clusters"` status when a side has fewer than two clusters in the local sample
  (clustered standard errors were zero there).
* The bandwidth floor now gives the `bwcheck`-th observation positive kernel weight, so
  points selected at the floor are estimated instead of reported as `"insufficient"`.
* `rdnet()` reports the selector's reason (`"insufficient"`, `"singular"`) instead of
  `"no bandwidth"`; `wbate()` defaults to equal weights over the estimated points.
* `density()` keeps control observations exactly on the boundary on the control side and
  gains `unique_locations`, with a warning when locations repeat (panels).
* Splitting an edge at the boundary divides the metric's own cost column, not only
  `length` and `travel_time` (`network_coordinates(extensive = )` adds others).
* `edge_frame()` keeps the cheapest of parallel edges, matching shortest paths.

## Fixes from cross-language testing

* Directed networks: locations partway along one-way edges were reached through either
  end node, as if the edge were two-way. Transition costs: the penalty for switching onto
  the class of a location's own edge was not charged. For directed networks (with
  `direction` other than `"both"`) and `transition_cost()`, locations and evaluation
  points are now placed on nodes (their edges are split there).
* Covariates (and `pooled = TRUE`) with different bandwidths on the two sides: each side's
  kernel weights are divided by its `h^2`, as in rd2d. Estimates with one bandwidth are
  unchanged.
* Splitting a network that was split before (such as `coords$network`, reused for a
  placebo boundary) no longer reuses the labels of nodes inserted the first time.
* Panels and other repeated measurements: the bandwidth floor (`bwcheck`) now counts distinct
  locations instead of observations.
* Input checks: matrices with several columns for `y`, `fuzzy` or `weights`, boundaries with
  coinciding vertices, non-finite `eval_points`, `n_eval < 1` and negative `snap_scale` are
  errors.
* With one evaluation point the uniform critical value is `qnorm()` itself, not simulated.
* `min_obs` counts distinct locations rather than rows.
* `boundary_grid()` on a closed ring no longer returns both 0 and the ring's length.
* Side inference in `boundary_from_points()` and `rd2dnet()` ignores observations beyond the
  ends of the line through `b`.
* `network_coordinates(treated = )` warns when treatment labels contradict the side of the
  boundary an observation lies on.
* `network_coordinates()` warns when the points lie far from the network (lon/lat points
  without `lonlat = TRUE`, or a CRS mismatch).
* Cluster-robust standard errors warn when some point has fewer than 10 clusters on a side.
* `plot(true_tau = )` accepts a single value or a function of arc position.
* HC2/HC3: `1 - h_ii` is clipped at 1e-8 (rd2d's value) instead of 1e-12, so observations
  with leverage 1 no longer add rounding noise to the standard errors.
* Nodes inserted where edges cross the boundary now get interpolated numeric node
  attributes, so `node_contrast()` costs on split edges are right.

## Checks and messages

* Errors for `level` given as a proportion, non-finite weights, cluster labels of the wrong
  length or with missing values, non-positive bandwidths, invalid `vce` in `rdbwnet()`, and
  treatment indicators that are not logical or 0/1. A length-2 `h` always means
  `c(h_control, h_treated)`.
* Warnings when both `cluster` and `conley_cutoff` are given, when evaluation points lie
  beyond the network's crossings or pile up on few crossings, and when nodes lie just
  outside the boundary tolerance. Unreachable observations are counted by rows and locations.

# rd2dnet 0.1.0

First release of the R package, a port of the Python package `rd2dnet` 0.1.0 with the
same methods and, on the bundled mock city, the same numbers.

* Network boundary discontinuity estimation (`rdnet()`): side-specific local polynomials
  in network coordinates, robust bias-corrected inference, uniform confidence bands,
  fuzzy designs, covariates, derivative and directional-derivative targets, HC0-HC3,
  cluster and network/Euclidean Conley standard errors.
* Aggregate effects as in rd2d: `wbate()` and `lbate()`, with the joint covariance across
  evaluation points stored in the result (`vcov()`).
* Bandwidth selection (`rdbwnet()`), ported from rd2d's direct plug-in, including the
  coverage-error-optimal variants.
* Density (sorting) test (`density()` method for network coordinates), ported from
  rddensity, with a street-length normalisation.
* User-defined network metrics: edge costs, travel time, node contrasts, transition
  penalties between road classes (`transition_cost()`) and precomputed matrices.
* Networks from data frames, adjacency matrices, igraph graphs or OpenStreetMap (via
  dodgr), with sf geometries accepted for points, boundaries, regions and masks.
* Shortest paths, nearest-segment snapping and the sparse Conley kernel are implemented
  in C++.
* Vignette with donut, covariate-balance, placebo-boundary and
  difference-in-discontinuities recipes.
