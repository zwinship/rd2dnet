segments_of <- function(net) {
  E <- edge_frame(net)
  P <- net$coords
  data.frame(x1 = P[E$u, 1], y1 = P[E$u, 2], x2 = P[E$v, 1], y2 = P[E$v, 2])
}

eval_b <- function(city) boundary_interpolate(city$boundary, boundary_grid(city$boundary, 9))

# --- rd2d-style entry ------------------------------------------------------------

test_that("rd2dnet() matches the step-by-step interface", {
  city <- mock_city()
  obs <- city$obs
  b <- eval_b(city)
  net <- city$network()
  res <- suppressWarnings(rd2dnet(obs$outcome, obs[, c("x", "y")], obs$treated, b, net,
                                  boundary = city$boundary, h = 900, seed = 0))
  co <- suppressWarnings(network_coordinates(net, city$boundary, obs[, c("x", "y")], eval_points = b))
  ref <- rdnet(obs$outcome, co, h = 900, seed = 0)
  expect_equal(res$table, ref$table)
  expect_true(is.finite(density(res$coords)$pvalue))
})

test_that("rd2dnet() infers the treated side", {
  city <- mock_city()
  obs <- city$obs
  b <- eval_b(city)
  X <- as.matrix(obs[, c("x", "y")])
  side <- boundary_from_points(b, X, obs$treated)$treated_side
  expect_false(identical(boundary_from_points(b[nrow(b):1, ], X, obs$treated)$treated_side, side))
  a <- suppressWarnings(rd2dnet(obs$outcome, X, obs$treated, b, city$network(), h = 900, seed = 0))
  a2 <- suppressWarnings(rd2dnet(obs$outcome, X, obs$treated, b, city$network(), treated_side = side,
                                 h = 900, seed = 0))
  expect_equal(a$table, a2$table)
  expect_true(all(a$coords$s[obs$treated] >= 0, na.rm = TRUE))
})

test_that("rd2dnet() accepts segment tables and igraph", {
  city <- mock_city()
  obs <- city$obs
  b <- eval_b(city)
  net <- city$network()
  ref <- suppressWarnings(rd2dnet(obs$outcome, obs[, c("x", "y")], obs$treated, b, net,
                                  boundary = city$boundary, h = 900, seed = 0))$table
  got <- suppressWarnings(rd2dnet(obs$outcome, obs[, c("x", "y")], obs$treated, b, segments_of(net),
                                  boundary = city$boundary, h = 900, seed = 0))$table
  expect_equal(got$tau_bc, ref$tau_bc, tolerance = 1e-10)
  expect_error(rd2dnet(obs$outcome, obs[, c("x", "y")], obs$treated, b, city$network), "call it")
})

test_that("rd2dnet() input errors", {
  city <- mock_city()
  obs <- city$obs
  b <- eval_b(city)
  net <- city$network()
  expect_error(rd2dnet(obs$outcome, obs[, c("x", "y")], rep(2, nrow(obs)), b, net), "0/1")
  expect_error(rd2dnet(obs$outcome, obs[, c("x", "y")], obs$treated[1:10], b, net), "same number")
  expect_error(boundary_from_points(b, obs[, c("x", "y")], rep(1, nrow(obs))), "both treated and control")
  expect_error(boundary_from_points(rbind(c(0, 0), c(0, 0)), treated_side = "left"), "two distinct")
})

test_that("side inference warns when t disagrees with the geometry", {
  city <- mock_city()
  set.seed(1)
  t <- sample(0:1, nrow(city$obs), replace = TRUE)
  expect_warning(boundary_from_points(eval_b(city), city$obs[, c("x", "y")], t), "side of the boundary")
})

# --- networks from lines -------------------------------------------------------------

test_that("a segment table reproduces the network", {
  city <- mock_city()
  net <- city$network()
  net2 <- network_graph_from_lines(segments_of(net))
  expect_equal(net2$n_nodes, net$n_nodes)
  expect_equal(nrow(net2$edges), nrow(edge_frame(net)))
  pts <- city$obs[1:200, c("x", "y")]
  c1 <- suppressWarnings(network_coordinates(net, city$boundary, pts))
  c2 <- suppressWarnings(network_coordinates(net2, city$boundary, pts))
  expect_equal(c1$s, c2$s, tolerance = 1e-8)
})

test_that("segment attributes are available to metrics", {
  seg <- data.frame(x1 = c(0, 100), y1 = 0, x2 = c(100, 200), y2 = 0, speed_kph = c(36, 3.6))
  net <- network_graph_from_lines(seg, travel_time())
  expect_equal(node_distance(net, 1, 3), 10 + 100)
})

test_that("nearby ends are joined and crossings split", {
  gap <- list(rbind(c(0, 0), c(100, 0)), rbind(c(100.4, 0), c(200, 0)))
  expect_equal(network_graph_from_lines(gap)$n_nodes, 4)
  expect_equal(network_graph_from_lines(gap, tol = 1)$n_nodes, 3)
  expect_error(network_graph_from_lines(gap, split_at_intersections = TRUE), "sf")
  skip_if_not_installed("sf")
  cross <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, -500), c(0, 500))),
                                                     sf::st_linestring(rbind(c(-500, 0), c(500, 0)))))
  expect_equal(network_graph_from_lines(cross)$n_nodes, 4)
  noded <- network_graph_from_lines(cross, split_at_intersections = TRUE)
  expect_equal(noded$n_nodes, 5)
  centre <- which.min(rowSums(abs(noded$coords)))
  expect_equal(sum(noded$edges$from == centre | noded$edges$to == centre), 4)
  expect_true(all(c("id") %in% names(noded$edges)))
})

test_that("curved lines and rings", {
  arc <- rbind(c(0, 0), c(50, 50), c(100, 0))
  net <- network_graph_from_lines(list(arc))
  expect_equal(node_distance(net, 1, 2), 2 * sqrt(2) * 50)
  ring <- network_graph_from_lines(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0))))
  expect_equal(ring$n_nodes, 2)
  expect_equal(nrow(ring$edges), 2)
  expect_warning(network_graph_from_lines(list(rbind(c(0, 0), c(0, 0)), rbind(c(0, 0), c(1, 0)))),
                 "zero-length")
})

test_that("directed lines follow the flow", {
  river <- list(rbind(c(0, 0), c(100, 0)), rbind(c(100, 0), c(200, 0)))
  net <- network_graph_from_lines(river, directed = TRUE, direction = "out")
  up <- which.min(net$coords[, 1])
  down <- which.max(net$coords[, 1])
  expect_equal(node_distance(net, up, down), 200)
  expect_equal(node_distance(net, down, up), Inf)
})

test_that("sf layers in lon/lat are projected and points follow", {
  skip_if_not_installed("sf")
  ll <- sf::st_sf(name = c("a", "b"), geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(-73.99, 40.75), c(-73.98, 40.75))),
    sf::st_linestring(rbind(c(-73.98, 40.75), c(-73.98, 40.76))), crs = 4326))
  net <- expect_no_warning(network_graph_from_lines(ll))
  expect_equal(net$crs, "EPSG:32618")
  total <- node_distance(net, 1, 3)
  expect_gt(total, 1800)
  expect_lt(total, 2000)
  expect_equal(net$edges$name, c("a", "b"))
  pts <- sf::st_sfc(sf::st_point(c(-73.985, 40.7501)), crs = 4326)
  P <- net$coords
  bnd <- boundary(rbind(c(P[2, 1] - 1, P[2, 2] - 500), c(P[2, 1] - 1, P[2, 2] + 2000)), treated_side = "left")
  co <- suppressWarnings(network_coordinates(net, bnd, pts))
  expect_true(all(is.finite(co$s)))
})

test_that("network_graph fills length and checks the metric column", {
  nodes <- data.frame(id = c("a", "b", "c"), x = c(0, 100, 100), y = c(0, 0, 50))
  net <- network_graph(data.frame(from = c("a", "b"), to = c("b", "c")), nodes)
  expect_equal(node_distance(net, "a", "c"), 150)
  expect_error(network_graph(data.frame(from = "a", to = "b"), nodes, metric = "minutes"), "available")
})

test_that("boundary_from_region uses the outline", {
  zone <- rect_polygon(1000, -1000, 3000, 1000)
  whole <- boundary_from_region(zone)
  expect_equal(whole$length, 8000)
  expect_true(whole$closed)
  part <- boundary_from_region(zone, within = rect_polygon(-1e6, -1e6, 1500, 1e6))
  expect_equal(part$length, 3000)
  expect_false(part$closed)
  expect_equal(boundary_is_treated(part, rbind(c(2000, 0), c(0, 0))), c(TRUE, FALSE))
})

test_that("lon/lat networks warn", {
  nodes <- data.frame(name = 1:2, x = c(-73.98, -73.97), y = c(40.75, 40.76))
  expect_warning(network_graph(data.frame(from = 1, to = 2, length = 1200), nodes), "longitude/latitude")
  expect_no_warning(network_graph(data.frame(from = 1, to = 2, length = 1200), nodes, crs = "EPSG:32618"))
})

test_that("two-column data frames are points", {
  city <- mock_city()
  net <- city$network()
  pts <- city$obs[1:50, c("x", "y")]
  names(pts) <- c("east", "north")
  c1 <- suppressWarnings(network_coordinates(net, city$boundary, pts))
  c0 <- suppressWarnings(network_coordinates(net, city$boundary, city$obs[1:50, c("x", "y")]))
  expect_equal(c1$s, c0$s)
})
