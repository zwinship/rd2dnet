# Changelog

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
