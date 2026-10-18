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
