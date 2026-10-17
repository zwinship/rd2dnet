# Changelog

## 0.1.0 (unreleased)

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
