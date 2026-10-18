## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Test environments

* GitHub Actions: ubuntu-latest (R release, devel, oldrel-1), macos-latest (R release),
  windows-latest (R release), all `R CMD check --as-cran`: Status OK.

## Notes for CRAN

* The package ports the bandwidth selector of 'rd2d' and the density estimator of
  'rddensity' (both GPL-3); their authors are listed as contributors and copyright
  holders in Authors@R.
* `network_graph_from_place()` downloads OpenStreetMap data and is only run in
  `\donttest{}`; offline it falls back to bundled data with a warning.
