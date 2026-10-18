# rd2dnet (R)

Boundary discontinuity designs on street networks. `rd2dnet` adapts the location-based
boundary discontinuity framework of [`rd2d`](https://github.com/rdpackages/rd2d)
(Cattaneo, Titiunik & Yu) from Euclidean coordinates to **network-native coordinates**
and **user-defined network metrics**. This is the R package; the Python package lives at
the root of the same repository and gives the same numbers on the bundled data.

> **Built on rd2d.** `rd2dnet` is a direct extension of the boundary discontinuity
> methodology and software of **Matias D. Cattaneo, Rocío Titiunik and Ruiqi Rae Yu**
> ([rd2d](https://github.com/rdpackages/rd2d), part of the
> [rdpackages](https://rdpackages.github.io/) project). The estimands, side-specific local
> polynomial estimation, robust bias-corrected inference, uniform bands and the
> bandwidth selector (a port of `rdbw2d`'s direct plug-in) follow their work. The
> density test follows Cattaneo, Jansson and Ma's `rddensity`. If you use `rd2dnet`,
> please cite the rd2d papers (`citation("rd2dnet")`), as well as this package.
> `rd2dnet` is an independent project and is not affiliated with or endorsed by the rd2d
> authors.

## Installation

```r
# development version, from GitHub
# install.packages("remotes")
remotes::install_github("zwinship/rd2dnet", subdir = "R-package")
```

The package needs a C++ compiler when installed from source (Rtools on Windows, Xcode
command line tools on macOS). Optional: `sf` (sf geometries, projecting lon/lat),
`igraph` (networks from igraph graphs), `dodgr` (OpenStreetMap download).

## Quick start

```r
library(rd2dnet)
city <- load_mock_midtown()                    # bundled, works offline
net <- city$network("length")                  # = network_graph(city$edges, city$nodes, "length")
coords <- network_coordinates(net, city$boundary, city$obs, n_eval = 21)

res <- rdnet(city$obs$outcome, coords,
             weights = city$obs$far,           # intensity / floor-area scaling
             cluster = city$obs$block)         # or conley_cutoff = 150
res                                            # pointwise BATEC table
plot(res)
wbate(res, as.numeric(res$table$status == "ok"))   # weighted average effect
lbate(res)                                         # largest effect

rdbwnet(city$obs$outcome, coords)              # bandwidths
density(coords)                                # sorting test per metre of street
rdnet(city$obs$outcome, coords, fuzzy = take_up)   # fuzzy design (take_up = your 0/1 or share)
```

`vignette("rd2dnet")` walks through every feature on the mock city (metrics, coordinates,
bandwidths, aggregate effects, covariates, fuzzy designs, derivatives, cluster and Conley
errors, the density test) and the recipes for donut designs, covariate balance, placebo
boundaries and difference-in-discontinuities.

## Your own data

A network is a data frame of edges (`from`, `to`, attributes such as `length`, optional
`geometry`) and a data frame of nodes (`name`, projected `x`, `y`):

```r
net <- network_graph(edges, nodes, metric = "length")
net <- network_graph_from_igraph(g)                    # igraph with x/y vertex attributes
net <- network_graph_from_adjacency(A, coords)
net <- network_graph_from_place("Manhattan, New York") # OpenStreetMap via dodgr + sf

bnd <- boundary(line, treated_region = zone,           # matrices or sf geometries
                mask = list(park_polygon, c(a, b)))    # masked stretches have no mass
```

## Distance is a metric you choose

```r
network_graph(edges, nodes, "length")                       # metres along streets
network_graph(edges, nodes, travel_time())                  # seconds, from speed_kph
network_graph(edges, nodes, node_contrast("far", 50))       # length + 50 * |land-use change|
network_graph(edges, nodes, function(edges) 2 * edges$length)   # any vectorised edge cost
network_graph(edges, nodes, precomputed(D))                 # your own node-by-node matrix

# Switching road class costs extra. Ramps are transparent connectors, so
# highway -> ramp -> highway is close while highway -> ramp -> local street is far.
network_graph(edges, nodes, transition_cost(
  "length", class_attr = "highway",
  penalty = data.frame(from = "motorway", to = "residential", penalty = 2000),
  connectors = c("motorway_link", "trunk_link")))
```

The coordinates `s` and `o`, the kernel distances and the bandwidths are all in the
metric's units. Treatment side is always geographic.

## Method

| | rd2d | rd2dnet |
|---|---|---|
| coordinates | `(X1 - b1, X2 - b2)` | `s_i`: signed network distance to the boundary (negative on the control side)<br>`o_i(b)`: signed offset along the boundary from `b` to where `i` reaches it |
| kernel distance | `‖X_i - b‖` | `d_i(b)`: shortest-path distance from `b` to `i` |
| estimate | side-specific intercept difference | same, with polynomial order `p` in `(s, o)` |
| inference | robust bias correction with `q = p + 1` | same, plus uniform bands |

See the README at the repository root for the full feature list, the Monte Carlo
validation, known limitations and references; all of it applies to the R package.

## Python and R

Function names follow the Python package in snake case: `NetworkGraph` is
`network_graph()`, `Boundary` is `boundary()`, `EdgeCost`/`TransitionCost`/`Precomputed`
are `edge_cost()`/`transition_cost()`/`precomputed()`, `RDNetResult.wbate()` is `wbate()`,
and `rd2dnet.density()` is the `density()` method for network coordinates. Indices
(evaluation points, nodes) are 1-based. Cost functions are vectorised over the edge table.
On the bundled mock city the two packages agree to about 1e-12; they can differ only in
which boundary node an observation is assigned to when two are exactly equally far.

## Development

```bash
Rscript -e 'devtools::test()'                          # offline test suite
RD2DNET_ONLINE=true Rscript -e 'devtools::test(filter = "online")'   # live OpenStreetMap download
```

The online test needs internet access and `dodgr`, `osmdata` and `sf`, like the Python
package's `pytest -m online`. `.github/workflows/online.yml` runs both on demand.

## License

GPL-3, the same license as rd2d and rddensity, from which parts of the implementation
are ported.
