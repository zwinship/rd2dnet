# Shared fixtures (computed once per test run).

.fixture_env <- new.env()

mock_city <- function() {
  if (is.null(.fixture_env$city)) .fixture_env$city <- load_mock_midtown()
  .fixture_env$city
}

mock_coords <- function() {
  if (is.null(.fixture_env$coords)) {
    city <- mock_city()
    .fixture_env$coords <- network_coordinates(city$network("length"), city$boundary, city$obs, n_eval = 21)
  }
  .fixture_env$coords
}

big_coords <- function() {
  if (is.null(.fixture_env$big)) {
    city <- make_mock_midtown(n_obs = 8000, seed = 21)
    coords <- network_coordinates(city$network(), city$boundary, city$obs,
                                  eval_points = c(600, 900, 2300, 2600))
    .fixture_env$big <- list(city = city, coords = coords)
  }
  .fixture_env$big
}

# 5 x 5 grid (100 m spacing) with a boundary along y = 200, treated below.
grid_fixture <- function() {
  ij <- expand.grid(j = 0:4, i = 0:4)
  nodes <- data.frame(name = seq_len(25) - 1L, x = 100 * ij$i, y = 100 * ij$j)
  id <- function(i, j) i * 5L + j
  e <- rbind(
    expand.grid(i = 0:3, j = 0:4)[, c("i", "j")] |> transform(a = id(i, j), b = id(i + 1L, j)),
    expand.grid(i = 0:4, j = 0:3)[, c("i", "j")] |> transform(a = id(i, j), b = id(i, j + 1L))
  )
  edges <- data.frame(from = e$a, to = e$b, length = 100)
  bnd <- boundary(rbind(c(-50, 200), c(450, 200)), treated_region = rect_polygon(-1e3, -1e3, 1e3, 200))
  list(nodes = nodes, edges = edges, boundary = bnd)
}

# Highway A -- ramp -- highway B, with a local street hanging off the ramp.
#
#   a0 --hw-- a1 --ramp-- r --ramp-- b1 --hw-- b0
#                         |
#                       local -- l1
highway_fixture <- function() {
  nodes <- data.frame(name = c("a0", "a1", "r", "b1", "b0", "l1"),
                      x = c(0, 100, 200, 300, 400, 200), y = c(0, 0, 0, 0, 0, 100))
  edges <- data.frame(from = c("a0", "a1", "r", "b1", "r"), to = c("a1", "r", "b1", "b0", "l1"),
                      length = 100,
                      highway = c("motorway", "motorway_link", "motorway_link", "motorway", "residential"))
  list(nodes = nodes, edges = edges)
}

# Reference: dense Bartlett kernel; co-located observations are at distance 0.
dense_hac <- function(coords, idx, Psi, cutoff) {
  Dm <- pairwise_distances(coords, idx, cutoff = cutoff)
  Dm[outer(coords$loc[idx], coords$loc[idx], "==")] <- 0
  K <- pmax(1 - Dm / cutoff, 0)
  out <- t(Psi) %*% K %*% Psi
  (out + t(out)) / 2
}
