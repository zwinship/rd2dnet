#' rd2dnet: Boundary Discontinuity Designs on Networks
#'
#' A network-native extension of the location-based boundary discontinuity
#' framework and software of rd2d (Cattaneo, Titiunik and Yu;
#' <https://rdpackages.github.io/rd2d/>), from Euclidean coordinates to
#' network coordinates and user-defined network metrics. Please cite the rd2d
#' papers when using this package (see `citation("rd2dnet")`).
#'
#' Typical workflow:
#' \preformatted{
#' net    <- network_graph(edges, nodes, metric = "length")  # or travel_time(), ...
#' bnd    <- boundary(line, treated_region = zone)
#' coords <- network_coordinates(net, bnd, points)
#' res    <- rdnet(y, coords)
#' }
#'
#' Main functions: [network_graph()] and the [metrics], [boundary()],
#' [network_coordinates()], [rdnet()] (with [wbate()] and [lbate()]),
#' [rdbwnet()] and the density test [density.rd2dnet_coords()].
#'
#' @keywords internal
#' @useDynLib rd2dnet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef vcov confint
"_PACKAGE"
