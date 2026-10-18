# Bundled offline mock city.
#
# Layout (metres): 11 north-south avenues 250 m apart ("primary"), 41 east-west
# streets 80 m apart numbered 40th-80th ("residential"), a west-side highway
# ("motorway") joined to the grid by ramps ("motorway_link"), and a park north
# of the boundary between avenues 3 and 7 with no interior streets or
# observations. The boundary is 60th Street; the treated zone lies south of it.

.MOCK <- list(
  avenue_spacing = 250, street_spacing = 80, n_avenues = 11, first_street = 40,
  last_street = 80, boundary_street = 60, highway_x = -300,
  ramp_streets = c(42, 50, 57, 66, 72, 79),
  speed_kph = c(residential = 25, primary = 40, motorway = 88, motorway_link = 40)
)
.MOCK$boundary_y <- (.MOCK$boundary_street - .MOCK$first_street) * .MOCK$street_spacing
.MOCK$park <- local({
  x <- c(3, 7, 7, 3, 3) * .MOCK$avenue_spacing
  y <- .MOCK$boundary_y + c(0, 0, 15, 15, 0) * .MOCK$street_spacing
  cbind(x, y, deparse.level = 0)
})
.MOCK$line <- rbind(c(.MOCK$highway_x - 100, .MOCK$boundary_y),
                    c((.MOCK$n_avenues - 1) * .MOCK$avenue_spacing + 100, .MOCK$boundary_y))

#' Mock city treatment effect
#'
#' The treatment effect used to simulate the mock city's outcome, by arc
#' length along its boundary: `1.5 + 0.75 sin(2 pi x / 2500)` with `x` the
#' easting.
#'
#' @param arc Arc-length positions along the mock boundary.
#' @return Numeric vector of true effects.
#' @examples
#' mock_true_tau(c(0, 500, 1000))
#' @export
mock_true_tau <- function(arc) {
  x <- as.numeric(arc) + .MOCK$line[1, 1]
  1.5 + 0.75 * sin(2 * pi * x / ((.MOCK$n_avenues - 1) * .MOCK$avenue_spacing))
}

.mock_boundary <- function(mask_park = TRUE) {
  boundary(.MOCK$line, treated_region = rect_polygon(-1e5, -1e5, 1e5, .MOCK$boundary_y),
           mask = if (mask_park) list(.MOCK$park) else NULL)
}

.mock_city <- function(nodes, edges, obs, park = TRUE) {
  force(nodes)
  force(edges)
  structure(list(
    nodes = nodes,
    edges = edges,
    obs = obs,
    boundary = .mock_boundary(park),
    park = if (park) .MOCK$park else NULL,
    boundary_y = .MOCK$boundary_y,
    true_tau = mock_true_tau,
    network = function(metric = "length", ...) {
      network_graph(edges, nodes, metric = metric, crs = "local-metres", ...)
    }
  ), class = "rd2dnet_mock_city")
}

#' @export
print.rd2dnet_mock_city <- function(x, ...) {
  cat(sprintf("<rd2dnet mock city> %d nodes, %d edges, %d observations; boundary y = %g%s\n",
              nrow(x$nodes), nrow(x$edges), nrow(x$obs), x$boundary_y,
              if (is.null(x$park)) "" else ", park masked"))
  cat("  $nodes $edges $obs $boundary $park $true_tau() $network(metric)\n")
  invisible(x)
}

#' Mock Midtown city
#'
#' `load_mock_midtown()` loads the bundled synthetic Midtown-Manhattan-style
#' street grid and 3,000 simulated parcels (no internet needed); it is the
#' same data as the Python package's bundled city. `make_mock_midtown()`
#' simulates a new city of the same design with R's random number generator.
#'
#' Layout (metres): 11 north-south avenues 250 m apart (`"primary"`), 41
#' east-west streets 80 m apart numbered 40th-80th (`"residential"`), a
#' west-side highway (`"motorway"`) joined to the grid by ramps
#' (`"motorway_link"`), and a park north of the boundary between avenues 3
#' and 7 with no interior streets or parcels. The boundary is 60th Street
#' (`y = 1600`); the treated zone lies south of it. The outcome is
#' `0.0008 (y - 1600) + 0.6 cos(x / 700) + 0.2 log(far) + tau(arc) T + N(0, 1)`
#' with `tau` given by [mock_true_tau()].
#'
#' @param n_obs Number of simulated parcels.
#' @param seed Random seed (the global random number stream is left untouched).
#' @param park Include the park (and mask its boundary frontage).
#' @return A `rd2dnet_mock_city` list with `nodes`, `edges`, `obs`, the
#'   [boundary()] `boundary`, the park polygon `park`, `boundary_y`,
#'   `true_tau` (a function of arc length) and `network(metric, ...)`, a
#'   function returning the [network_graph()] with the given metric.
#' @seealso [mock_midtown_obs], [mock_midtown_nodes], [mock_midtown_edges].
#' @examples
#' city <- load_mock_midtown()
#' city
#' head(city$obs)
#' city$network("length")
#' small <- make_mock_midtown(n_obs = 500, seed = 1)
#' @export
load_mock_midtown <- function() {
  env <- new.env()
  utils::data(list = c("mock_midtown_nodes", "mock_midtown_edges", "mock_midtown_obs"),
              package = "rd2dnet", envir = env)
  .mock_city(env$mock_midtown_nodes, env$mock_midtown_edges, env$mock_midtown_obs)
}

#' @rdname load_mock_midtown
#' @export
make_mock_midtown <- function(n_obs = 3000, seed = 0, park = TRUE) {
  .with_seed(seed, {
    g <- .make_mock_graph(park)
    xmax <- (.MOCK$n_avenues - 1) * .MOCK$avenue_spacing
    ymax <- (.MOCK$last_street - .MOCK$first_street) * .MOCK$street_spacing
    pb <- c(range(.MOCK$park[, 1]), range(.MOCK$park[, 2]))
    pts <- matrix(numeric(0), 0, 2)
    while (nrow(pts) < n_obs) {
      cand <- cbind(stats::runif(n_obs, 0, xmax), stats::runif(n_obs, 0, ymax))
      inside <- cand[, 1] > pb[1] & cand[, 1] < pb[2] & cand[, 2] > pb[3] & cand[, 2] < pb[4]
      pts <- rbind(pts, if (park) cand[!inside, , drop = FALSE] else cand)
    }
    x <- pts[seq_len(n_obs), 1]
    y <- pts[seq_len(n_obs), 2]
    treated <- y < .MOCK$boundary_y
    arc <- x - .MOCK$line[1, 1]
    far <- round(exp(stats::rnorm(n_obs, 1.2, 0.5)), 3)
    f <- 0.0008 * (y - .MOCK$boundary_y) + 0.6 * cos(x / 700) + 0.2 * log(far)
    outcome <- f + mock_true_tau(arc) * treated + stats::rnorm(n_obs)
    obs <- data.frame(
      x = round(x, 3), y = round(y, 3), treated = treated, far = far,
      block = as.integer(floor(x / .MOCK$avenue_spacing) * 1000 + floor(y / .MOCK$street_spacing)),
      outcome = round(outcome, 6)
    )
    .mock_city(g$nodes, g$edges, obs, park)
  })
}

.make_mock_graph <- function(park = TRUE) {
  M <- .MOCK
  streets <- M$first_street:M$last_street
  in_park <- function(i, j) park && i > 3 && i < 7 && j > M$boundary_street && j < M$boundary_street + 15
  grid <- expand.grid(j = streets, i = 0:(M$n_avenues - 1))[, c("i", "j")]
  grid <- grid[!mapply(in_park, grid$i, grid$j), ]
  mu <- 1.5 - 0.02 * abs(grid$j - 52)
  far <- round(exp(stats::rnorm(nrow(grid), mu, 0.4)), 3)
  nodes <- data.frame(name = grid$i * 100 + grid$j, x = grid$i * M$avenue_spacing,
                      y = (grid$j - M$first_street) * M$street_spacing, far = far)
  hw <- data.frame(name = 9000 + streets, x = M$highway_x,
                   y = (streets - M$first_street) * M$street_spacing, far = 0)
  nodes <- rbind(nodes, hw)
  ids <- nodes$name
  from <- c()
  to <- c()
  cls <- c()
  add <- function(a, b, k) {
    from <<- c(from, a)
    to <<- c(to, b)
    cls <<- c(cls, k)
  }
  for (r in seq_len(nrow(grid))) {
    id <- grid$i[r] * 100 + grid$j[r]
    east <- id + 100
    north <- id + 1
    if (grid$i[r] + 1 < M$n_avenues && east %in% ids) add(id, east, "residential")
    if (grid$j[r] + 1 <= M$last_street && north %in% ids) add(id, north, "primary")
  }
  for (j in streets[-1]) add(9000 + j - 1, 9000 + j, "motorway")
  for (j in M$ramp_streets) add(9000 + j, j, "motorway_link")
  xa <- nodes$x[match(from, ids)]
  ya <- nodes$y[match(from, ids)]
  xb <- nodes$x[match(to, ids)]
  yb <- nodes$y[match(to, ids)]
  len <- sqrt((xa - xb)^2 + (ya - yb)^2)
  speed <- unname(M$speed_kph[cls])
  edges <- data.frame(from = from, to = to, length = len, highway = cls, speed_kph = speed,
                      travel_time = len / (speed / 3.6), stringsAsFactors = FALSE)
  list(nodes = nodes, edges = edges)
}

#' Mock Midtown data
#'
#' The bundled synthetic city used in examples, tests and the vignette (see
#' [load_mock_midtown()], which assembles them with the boundary).
#'
#' @format
#' `mock_midtown_obs`: 3,000 simulated parcels with columns
#' \describe{
#'   \item{x, y}{projected coordinates (metres)}
#'   \item{treated}{south of the boundary (60th Street, `y < 1600`)}
#'   \item{far}{floor-area ratio, used as intensity weight}
#'   \item{block}{city block id, for clustering}
#'   \item{outcome}{simulated outcome}
#' }
#'
#' `mock_midtown_nodes`: 450 street intersections with columns `name`, `x`,
#' `y` and a node-level floor-area ratio `far`.
#'
#' `mock_midtown_edges`: street segments with columns `from`, `to`, `length`
#' (metres), `highway` (road class), `speed_kph` and `travel_time` (seconds).
#' @source Simulated; identical to the data bundled with the Python package
#'   `rd2dnet`.
#' @name mock_midtown
#' @aliases mock_midtown_obs mock_midtown_nodes mock_midtown_edges
#' @keywords datasets
NULL
