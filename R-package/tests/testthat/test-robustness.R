# Edge cases found by trying the package on many kinds of networks and data.

path_network <- function(n = 41, spacing = 250) {
  nodes <- data.frame(name = seq_len(n), x = (seq_len(n) - 1 - n %/% 2) * spacing, y = 0)
  edges <- data.frame(from = seq_len(n - 1), to = 2:n, length = spacing)
  network_graph(edges, nodes)
}

# Three parallel rivers crossing a border at x = 0: o has no variation near each crossing.
rivers <- function(seed = 0, n = 3000) {
  set.seed(seed)
  xs <- seq(-6000, 6000, length.out = 25)
  lines <- list()
  for (yy in c(-3000, 0, 3000)) {
    for (i in 1:24) lines[[length(lines) + 1]] <- rbind(c(xs[i], yy), c(xs[i + 1], yy))
  }
  net <- network_graph_from_lines(lines)
  r <- sample(1:3, n, replace = TRUE)
  x <- stats::runif(n, -6000, 6000)
  y <- c(-3000, 0, 3000)[r] + stats::rnorm(n, 0, 30)
  t <- x > 0
  Y <- 1 + 0.0002 * x - 1.5 * t + stats::rnorm(n, 0, 0.5)
  bd <- boundary(rbind(c(0, -4000), c(0, 4000)), treated_side = "right")
  list(net = net, bd = bd, X = cbind(x, y), t = t, Y = Y)
}

test_that("few crossings fall back to a distance-only fit", {
  rv <- rivers()
  b <- rbind(c(0, -3000), c(0, 0), c(0, 3000))
  res <- rd2dnet(rv$Y, rv$X, rv$t, b, rv$net, boundary = rv$bd, seed = 0)
  expect_true(all(res$table$status == "ok"))
  expect_true(all(res$table$basis == "s"))
  expect_true(all(res$bandwidths$basis == "s"))
  expect_lt(abs(wbate(res)[["tau_bc"]] + 1.5), 0.5)
  r <- rdnet(mock_city()$obs$outcome, mock_coords(), h = 900, seed = 0)
  expect_equal(unique(r$table$basis[r$table$status == "ok"]), "s,o")
})

test_that("many evaluation points on few crossings warn", {
  rv <- rivers(n = 600)
  expect_warning(network_coordinates(rv$net, rv$bd, rv$X, n_eval = 15), "distinct places")
})

test_that("a single cluster per side is flagged", {
  co <- mock_coords()
  cl <- ifelse(co$treated, "A", "B")
  r <- rdnet(mock_city()$obs$outcome, co, h = 900, cluster = cl, seed = 0)
  expect_true("few clusters" %in% r$table$status)
  expect_true(all(r$table$status %in% c("few clusters", "masked", "insufficient")))
})

test_that("the bandwidth floor keeps points estimable", {
  co <- mock_coords()
  y <- mock_city()$obs$outcome
  bw <- rdbwnet(y, co)
  ok <- bw$status == "ok"
  h <- ifelse(ok, pmax(bw$hmin0, bw$hmin1), 1000)
  r <- rdnet(y, co, h = h, seed = 0)
  expect_false(any(r$table$status[ok] == "insufficient"))
})

test_that("closed rings: offsets wrap and the grid does not repeat", {
  ij <- expand.grid(j = 0:20, i = 0:20)
  nodes <- data.frame(name = seq_len(nrow(ij)), x = 100 * ij$i, y = 100 * ij$j)
  id <- function(i, j) i * 21L + j + 1L
  e <- rbind(expand.grid(i = 0:19, j = 0:20) |> transform(a = id(i, j), b = id(i + 1L, j)),
             expand.grid(i = 0:20, j = 0:19) |> transform(a = id(i, j), b = id(i, j + 1L)))
  net <- network_graph(data.frame(from = e$a, to = e$b, length = 100), nodes)
  bnd <- boundary_from_region(rect_polygon(550, 550, 1450, 1450))
  expect_true(bnd$closed)
  set.seed(0)
  pts <- matrix(stats::runif(4000, 0, 2000), ncol = 2)
  co <- network_coordinates(net, bnd, pts, n_eval = 12)
  expect_gt(min(diff(sort(co$eval_arc))), 0.5 * bnd$length / 12)
  j <- which.min(co$eval_arc)
  at <- coords_at(co, j)
  arc <- co$crossing_arc
  L <- bnd$length
  near <- is.finite(at$o) & at$d < 400
  before <- near & arc > L / 2
  after <- near & arc < L / 2 & arc > co$eval_node_arc[j]
  expect_true(all(at$o[before] <= 0))
  expect_true(all(at$o[after] >= 0))
})

test_that("density keeps boundary controls and flags repeated locations", {
  city <- mock_city()
  idx <- rep(seq_len(nrow(city$obs)), each = 3)
  c3 <- network_coordinates(city$network(), city$boundary, city$obs[idx, c("x", "y")], n_eval = 5)
  expect_warning(density(c3), "share a location")
  d1 <- density(mock_coords())
  d3 <- density(c3, unique_locations = TRUE)
  expect_equal(d3$pvalue, d1$pvalue, tolerance = 1e-6)
})

test_that("split edges divide the metric's own column", {
  net <- network_graph_from_adjacency(rbind(c(0, 200), c(200, 0)), rbind(c(-100, 0), c(100, 0)))
  co <- network_coordinates(net, boundary(rbind(c(0, -50), c(0, 50)), treated_side = "right"),
                            rbind(c(-100, 0), c(100, 0)), n_eval = 1)
  expect_equal(co$s, c(-100, 100))
})

test_that("parallel edges use the cheaper one", {
  out <- lapply(list(c(300, 100), c(100, 300)), function(len) {
    net <- network_graph(data.frame(from = c("a", "a"), to = c("b", "b"), length = len),
                         data.frame(name = c("a", "b"), x = 0, y = c(-100, 100)))
    network_coordinates(net, boundary(rbind(c(-50, 0), c(50, 0)), treated_side = "left"),
                        rbind(c(0, 50), c(0, -50)), eval_points = 50)$s
  })
  expect_equal(out[[1]], out[[2]])
})

test_that("input validation", {
  co <- mock_coords()
  y <- mock_city()$obs$outcome
  expect_error(rdnet(y, co, h = 900, level = 0.95), "percentage")
  w <- rep(1, length(y))
  w[3] <- Inf
  expect_error(rdnet(y, co, h = 900, weights = w), "finite")
  expect_error(rdnet(y, co, h = 900, cluster = 1:5), "one value per observation")
  expect_error(rdnet(y, co, h = 900, cluster = ifelse(seq_along(y) < 5, NA, 1)), "missing")
  expect_error(rdnet(y, co, h = -1), "positive")
  expect_error(rdbwnet(y, co, vce = "bogus"), "vce")
  city <- mock_city()
  expect_error(network_coordinates(city$network(), city$boundary, city$obs[, c("x", "y")],
                                   treated = ifelse(city$obs$treated, "yes", "no")), "logical or 0/1")
  expect_warning(rdnet(y, co, h = 900, cluster = seq_along(y) %% 50, conley_cutoff = 100, seed = 0), "Conley")
})

test_that("h = c(h0, h1) with two evaluation points", {
  city <- mock_city()
  c2 <- network_coordinates(city$network(), city$boundary, city$obs, eval_points = c(600, 2600))
  r <- rdnet(city$obs$outcome, c2, h = c(800, 1000), seed = 0)
  expect_true(all(r$table$h0 == 800) && all(r$table$h1 == 1000))
  r2 <- rdnet(city$obs$outcome, c2, h = cbind(c(800, 1000), c(800, 1000)), seed = 0)
  expect_equal(r2$table$h0, c(800, 1000))
})

test_that("status says why selection failed", {
  net <- path_network()
  bd <- boundary(rbind(c(0, -500), c(0, 500)), treated_side = "right")
  set.seed(0)
  x <- c(stats::runif(300, -5000, 0), stats::runif(10, 0, 5000))
  co <- network_coordinates(net, bd, cbind(x, 5), eval_points = 500)
  r <- rdnet(stats::rnorm(length(x)), co)
  expect_equal(r$table$status[1], "insufficient")
})

test_that("wbate defaults to equal weights", {
  r <- rdnet(mock_city()$obs$outcome, mock_coords(), h = 900, seed = 0)
  expect_equal(wbate(r), wbate(r, as.numeric(r$table$status == "ok")))
})

test_that("travel_time default covers missing speeds", {
  seg <- data.frame(x1 = c(0, 100), y1 = 0, x2 = c(100, 200), y2 = 0, speed_kph = c(36, NA))
  net <- network_graph_from_lines(seg, travel_time(default_kph = 18))
  expect_equal(node_distance(net, 1, 3), 10 + 20)
})

test_that("boundary() merges connected multilines and rejects polygons", {
  skip_if_not_installed("sf")
  ml <- sf::st_multilinestring(list(rbind(c(0, 0), c(1, 0)), rbind(c(1, 0), c(2, 0))))
  expect_equal(boundary(ml, treated_side = "left")$length, 2)
  poly <- sf::st_polygon(list(rect_polygon(0, 0, 1, 1)))
  expect_error(boundary(poly, treated_region = poly), "boundary_from_region")
})

test_that("a boundary just off the street warns about tol", {
  nodes <- data.frame(name = 1:7, x = c(-100, 0, 100, 0, 0, -100, 100), y = c(0, 0, 0, 100, -100, 100, 100))
  edges <- data.frame(from = c(1, 2, 4, 5, 6, 7), to = c(2, 3, 2, 2, 1, 3), length = 100)
  bd <- boundary(rbind(c(-500, 0.01), c(500, 0.01)), treated_side = "left")
  expect_warning(network_coordinates(network_graph(edges, nodes), bd, rbind(c(0, 50), c(0, -50)), n_eval = 3),
                 "increase tol")
})

test_that("igraph input: weight falls back to length, directed override, vertex geometry", {
  skip_if_not_installed("igraph")
  g <- igraph::make_graph(c(1, 2, 2, 3), directed = TRUE)
  igraph::V(g)$x <- c(0, 100, 200)
  igraph::V(g)$y <- 0
  igraph::E(g)$length <- c(100, 100)
  net <- network_graph_from_igraph(g)
  expect_equal(node_distance(net, 1, 3), 200)
  net2 <- network_graph_from_igraph(g, directed = FALSE)
  expect_false(net2$directed)
  skip_if_not_installed("sf")
  g2 <- igraph::make_graph(c(1, 2), directed = FALSE)
  igraph::V(g2)$geometry <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(30, 40)), crs = 32618)
  net3 <- network_graph_from_igraph(g2, metric = "length")
  expect_equal(node_distance(net3, 1, 2), 50)
  expect_equal(net3$crs, "EPSG:32618")
})

test_that("inputs that used to pass silently are rejected", {
  cc <- mock_coords()
  y <- mock_city()$obs$outcome
  expect_error(rdnet(cbind(y, y)[seq_len(cc$n_obs / 2), ], cc, h = 600), "must be a vector")
  expect_error(density(cc, weights = rep(-1, cc$n_obs)), "finite and non-negative")
  expect_error(boundary(rbind(c(500, 1000), c(500, 1000)), treated_side = "left"), "zero length")
  city <- mock_city()
  net <- city$network("length")
  expect_error(network_coordinates(net, city$boundary, city$obs, eval_points = c(NA, 500)), "finite")
  expect_error(network_coordinates(net, city$boundary, city$obs, n_eval = 0), "positive integer")
  nodes <- data.frame(name = 1:2, x = c(0, 100), y = 0)
  expect_error(network_graph(data.frame(from = 1, to = 2, length = 100), nodes, snap_scale = -5), "snap_scale")
})
