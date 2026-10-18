test_that("edges crossing the boundary are split", {
  g <- grid_fixture()
  b <- boundary(rbind(c(-50, 150), c(450, 150)), treated_region = rect_polygon(-1e3, -1e3, 1e3, 150))
  net <- with_boundary_crossings(network_graph(g$edges, g$nodes), b)
  expect_equal(attr(net, "n_crossings_inserted"), 5)
  cr <- which(net$nodes$rdnet_crossing)
  expect_length(cr, 5)
  expect_equal(net$nodes$y[cr], rep(150, 5))
  for (k in cr) {
    lens <- sort(net$edges$length[net$edges$from == k | net$edges$to == k])
    expect_equal(lens, c(50, 50))
  }
  # the split network is cached for the same boundary
  net0 <- network_graph(g$edges, g$nodes)
  expect_identical(with_boundary_crossings(net0, b), with_boundary_crossings(net0, b))
})

test_that("splitting keeps edge geometries and coordinates work across the split", {
  nodes <- data.frame(name = 1:3, x = c(0, 0, 100), y = c(-100, 100, 100))
  geom <- list(rbind(c(0, -100), c(0, 0), c(0, 100)), NULL)
  edges <- data.frame(from = c(1, 2), to = c(2, 3), length = c(200, 100))
  edges$geometry <- geom
  b <- boundary(rbind(c(-50, 50), c(150, 50)), treated_side = "right")
  net <- with_boundary_crossings(network_graph(edges, nodes), b)
  expect_equal(attr(net, "n_crossings_inserted"), 1)
  expect_equal(sort(net$edges$length[1:3]), c(50, 100, 150))
  cc <- network_coordinates(network_graph(edges, nodes), b, rbind(c(0, -50), c(50, 100)))
  expect_equal(cc$s, c(100, -100))
})

test_that("snapping to an edge splits its cost", {
  g <- grid_fixture()
  sn <- snap_points(network_graph(g$edges, g$nodes), rbind(c(30, 10)))
  expect_equal(sn$depth, 10)
  expect_equal(sort(c(sn$cost_u, sn$cost_v)), c(30, 70))
})

test_that("snapping respects allowed nodes", {
  g <- grid_fixture()
  net <- network_graph(g$edges, g$nodes)
  sn <- snap_points(net, rbind(c(0, 0)), allowed = net$coords[, 2] >= 300, mode = "node")
  expect_equal(net$coords[sn$u, 2], 300)
  expect_error(snap_points(net, rbind(c(0, 0)), mode = "bogus"), "mode")
})

test_that("networks from adjacency matrices", {
  A <- rbind(c(0, 2, 0), c(2, 0, 3), c(0, 3, 0))
  net <- network_graph_from_adjacency(A, rbind(c(0, 0), c(1, 0), c(2, 0)))
  expect_equal(node_distance(net, 1, 3), 5)
  expect_equal(nrow(net$edges), 2)
  skip_if_not_installed("Matrix")
  net2 <- network_graph_from_adjacency(Matrix::Matrix(A, sparse = TRUE), rbind(c(0, 0), c(1, 0), c(2, 0)))
  expect_equal(node_distance(net2, 1, 3), 5)
})

test_that("networks from igraph", {
  skip_if_not_installed("igraph")
  g <- igraph::make_ring(4)
  igraph::V(g)$x <- c(0, 1, 1, 0)
  igraph::V(g)$y <- c(0, 0, 1, 1)
  igraph::E(g)$weight <- c(1, 2, 3, 4)
  net <- network_graph_from_igraph(g)
  expect_equal(node_distance(net, 1, 3), 3)
  expect_output(print(net), "nodes=4")
})

test_that("from_place falls back to the mock network offline", {
  local_mocked_bindings(.osm_download = function(...) stop("no internet"))
  expect_warning(net <- network_graph_from_place("Manhattan, New York"), "bundled mock")
  expect_equal(net$n_nodes, 450)
  expect_error(network_graph_from_place("Manhattan", offline_fallback = FALSE), "no internet")
  expect_error(network_graph_from_place(), "query or a bbox")
})

test_that("from_place builds a projected network from a dodgr table", {
  skip_if_not_installed("sf")
  fake <- data.frame(from_id = c("a", "b"), from_lon = c(-73.98, -73.979), from_lat = c(40.76, 40.76),
                     to_id = c("b", "c"), to_lon = c(-73.979, -73.979), to_lat = c(40.76, 40.761),
                     d = c(84.3, 111.2), highway = "residential", time = c(60, 80))
  local_mocked_bindings(.osm_download = function(...) fake)
  net <- network_graph_from_place("somewhere")
  expect_equal(net$n_nodes, 3)
  expect_equal(net$crs, "EPSG:32618")
  expect_equal(node_distance(net, "a", "c"), 84.3 + 111.2)
})

test_that("missing coordinates are an error", {
  expect_error(network_graph(data.frame(from = 1, to = 2), data.frame(name = 1:2)), "coordinates")
  expect_error(network_graph(data.frame(from = 1, to = 3), data.frame(name = 1:2, x = 0, y = 0)), "node id")
})

test_that("edge_frame keeps one row per street and orients geometries", {
  nodes <- data.frame(name = 1:2, x = c(0, 10), y = 0)
  edges <- data.frame(from = c(1, 2), to = c(2, 1), length = c(10, 12))
  edges$geometry <- list(rbind(c(10, 0), c(5, 1), c(0, 0)), NULL)
  E <- edge_frame(network_graph(edges, nodes))
  expect_equal(nrow(E), 1)
  expect_equal(E$geometry[[1]][1, ], c(0, 0))
  expect_equal(E$length, 10)
})
