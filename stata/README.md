# rd2dnet for Stata

Boundary discontinuity designs on street networks: the Stata version of the
Python and R packages in this repository. It carries the location-based framework
of [rd2d](https://rdpackages.github.io/rd2d/) (Cattaneo, Titiunik and Yu) over to
network coordinates and user-defined network distances. Please cite the rd2d papers
when you use it (see `help rd2dnet`).

> **Status:** tested in Stata 19 (SE). `test/rd2dnet_tests.do` checks every feature
> against reference numbers from Python and R and passes (`ALL TESTS PASSED`), and
> `test/install_check.do` installs the package with `net install` and runs every
> help-file example. The code uses only Stata 16 features, but has not been run in
> Stata 16-18. The Stata version mirrors the R port function by function (the R port
> matches the Python package to about 1e-12).

## Installation

Requires Stata 16 or newer (frames).

```stata
net install rd2dnet, from("https://raw.githubusercontent.com/zwinship/rd2dnet/main/stata") replace
```

Once it is on SSC: `ssc install rd2dnet`.

## Quick start

```stata
rdnetmock, clear                      // bundled mock city: parcels + frames rdn_nodes, rdn_edges
rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) boundary(rdn_line) ///
    treatedregion(rdn_region) maskpolygons(rdn_park) neval(21)
rdnet outcome [aw=far], vce(cluster block)   // BATEC with bandwidths from rdbwnet
estat wbate, weights(equal)                  // weighted average effect
estat lbate                                  // largest effect
rdbwnet outcome, bwselect(msetwo)            // bandwidths
rdnetdensity                                 // sorting test per metre of street
rdnet outcome, h(700) fuzzy(takeup)          // fuzzy design
```

Your own data, with the observations in memory:

```stata
* data entry as in rd2d: points b along the boundary (in order) and a treatment variable
matrix b = (100, 1600 \ 1000, 1600 \ 2000, 1600 \ 2900, 1600)
rdnetset x y, nodes(rdn_nodes) edges(rdn_edges) evalpoints(b) treated(treated)

* a GIS road layer exported as segments with end coordinates: no node table needed
frame create roads
frame roads: import delimited roads.csv        // x1 y1 x2 y2 length speed ...
rdnetset x y, edges(roads) jointol(0.5) boundary(line) treatedside(left)
```

Otherwise put nodes (`id x y`) and edges (`from to length ...`) in two frames, the
boundary vertices in a k x 2 matrix (first = last vertex for a zone outline), and the
treated zone in a polygon matrix (or use `treatedside(left|right)`). Where the network
meets the boundary at only a few points (rivers, rail lines, highways), those points are
fitted in the distance to the boundary alone (the `basis` column of `e(results)` is 2).
Use `rdnetdensity, unique` with repeated measurements.

| | Python / R | Stata |
|---|---|---|
| network | `NetworkGraph` / `network_graph()` | node and edge frames in `rdnetset` |
| network from line segments | `NetworkGraph.from_lines()` / `network_graph_from_lines()` | `rdnetset ..., edges()` with `x1 y1 x2 y2`, no `nodes()` |
| rd2d-style entry | `rd2dnet(Y, X, t, b, network)` | `rdnetset x y, evalpoints(b) treated(t)`, then `rdnet y` |
| metric | `EdgeCost`, `travel_time()`, `node_contrast()` | `cost(varname)`: build any cost variable in Stata |
| transition penalties | `TransitionCost` | `classvar()`, `penalty()`, `penaltymatrix()`, `connectors()` |
| coordinates | `network_coordinates()` | `rdnetset` (creates `_rdn_s`, `_rdn_treated`, `_rdn_arc`) |
| estimation | `rdnet()` | `rdnet` (results in `e()`; `e(b)`, `e(V)` work with `test`/`lincom`) |
| aggregates | `wbate()`, `lbate()` | `estat wbate`, `estat lbate` |
| nearest-crossing bins | `bin_coordinates()`, `rdnet_fe()`, `wbate("density"/"gls")` | `rdnetset ..., features()`, `rdnetbins`, `rdnetfe`, `estat wbate, weights(density\|gls)` |
| bandwidths | `rdbwnet()` | `rdbwnet` |
| density test | `density()` | `rdnetdensity` |
| mock data | `load_mock_midtown()` | `rdnetmock` |

Not in the Stata version yet: precomputed distance matrices, curved edge geometries
(edges are straight segments between nodes), and OpenStreetMap download. Coordinates
live in memory, so rerun `rdnetset` after reloading the data.

## Files

| file | contents |
|---|---|
| `lrd2dnet.mata` | the Mata library: shortest paths, geometry, coordinates, estimation, bandwidths, density test, Conley HAC. Compiled at run time by `rd2dnet_mata.ado` (no `.mlib` to build). |
| `rdnetset.ado`, `rdnet.ado`, `rdnet_estat.ado`, `rdbwnet.ado`, `rdnetdensity.ado`, `rdnetmock.ado` | commands (thin wrappers that parse options and call Mata entry points `rdn_*_entry()`) |
| `*.sthlp` | help files |
| `rdnet_mock_*.dta` | the bundled mock city (same data as the Python and R packages) |
| `rd2dnet.pkg`, `stata.toc` | `net install` manifest |
| `test/rd2dnet_tests.do` | test suite with reference values |

## Development notes

* **Running the tests:** `cd stata/test`, then `do rd2dnet_tests.do`. Without a
  local clone, `do https://raw.githubusercontent.com/zwinship/rd2dnet/<sha>/stata/test/fetch.do <sha>`
  downloads the files at a commit into `~/rd2dnet_test`. To check an installation,
  `do https://raw.githubusercontent.com/zwinship/rd2dnet/<sha>/stata/test/install_check.do <sha>`.
  Use commit SHAs: raw.githubusercontent.com caches branch URLs for several minutes.
* **Test suite details:** the suite adds the
  parent directory to the adopath, drops any previously compiled `rdn_*()` Mata
  functions, runs every block even when one fails, and prints a summary. Fix
  compile errors in `lrd2dnet.mata` first (the first command that runs compiles
  it), then work through the FAIL lines.
* **Design:** each ado-file parses its options into locals, then calls a Mata
  entry point, which reads them with `st_local()`. The network coordinates are
  kept in the Mata external `rdn_C` (a `struct rdn_coords`) between commands, and
  observations are matched to them through the `_rdn_id` variable.
* **Mirror of the R port:** the Mata functions follow `R-package/R/*.R` and
  `R-package/src/core.cpp` closely (same names without the dot, same steps), so
  when a number disagrees, compare with the R function of the same name.
  Internally, infinite distances are `1e300`; they become missing values before
  results reach Stata.
* **Reference values** in the tests come from Python where results are identical,
  and from R for automatic bandwidths, because Python breaks exact ties between
  equidistant boundary nodes differently (R and this port send ties to the first
  boundary node).
* **Recompiling during development:** after editing `lrd2dnet.mata`, run
  `mata: mata drop rdn_*()` (or `clear all`) so `rd2dnet_mata` recompiles it.

## Publishing to SSC

SSC packages are submitted by email to Kit Baum (kit.baum@bc.edu) with the files
attached (a zip of this directory without `test/`), after the test suite passes.
Update `Distribution-Date` in `rd2dnet.pkg` and the version lines in the ado-files
and help files for each release.

## License

GPL-3, the same license as rd2d and rddensity, from which parts of the implementation
are ported.
