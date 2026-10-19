hw_net <- function(metric = "length", ...) {
  h <- highway_fixture()
  network_graph(h$edges, h$nodes, metric = metric, ...)
}

test_that("edge length is the sum of lengths", {
  net <- hw_net()
  expect_equal(node_distance(net, "a0", "b0"), 400)
  expect_equal(node_distance(net, "a0", "l1"), 300)
})

test_that("a cost function is applied to all edges", {
  expect_equal(node_distance(hw_net(function(edges) 2 * edges$length), "a0", "b0"), 800)
  expect_error(node_distance(hw_net(function(edges) 1), "a0", "b0"), "one value per edge")
})

test_that("travel time metric", {
  nodes <- data.frame(name = 0:1, x = c(0, 1000), y = 0)
  net <- network_graph(data.frame(from = 0, to = 1, length = 1000, speed_kph = 36), nodes, travel_time())
  expect_equal(node_distance(net, 0, 1), 100) # 1 km at 10 m/s
  net2 <- network_graph(data.frame(from = 0, to = 1, length = 1000), nodes, travel_time(default_kph = 3.6))
  expect_equal(node_distance(net2, 0, 1), 1000)
  expect_error(node_distance(network_graph(data.frame(from = 0, to = 1, length = 1), nodes, travel_time()), 0, 1),
               "no 'speed_kph'")
})

test_that("transition costs: highway -> ramp -> highway stays close", {
  m <- transition_cost("length", class_attr = "highway", default_penalty = 5000, connectors = "motorway_link")
  net <- hw_net(m)
  expect_equal(node_distance(net, "a0", "b0"), 400) # no class change
  expect_equal(node_distance(net, "a0", "l1"), 300 + 5000) # motorway -> residential
  expect_output(print(m), "transition_cost")
})

test_that("transition penalties from a data frame, a matrix or a function; asymmetric", {
  pen <- data.frame(from = c("motorway", "residential"), to = c("residential", "motorway"), penalty = c(1000, 10))
  net <- hw_net(transition_cost("length", penalty = pen, connectors = "motorway_link"))
  expect_equal(node_distance(net, "a0", "l1"), 1300)
  expect_equal(node_distance(net, "l1", "a0"), 310)
  M <- matrix(c(0, 10, 1000, 0), 2, dimnames = list(c("motorway", "residential"), c("motorway", "residential")))
  net_m <- hw_net(transition_cost("length", penalty = M, connectors = "motorway_link"))
  expect_equal(node_distance(net_m, "a0", "l1"), 1300)
  f <- function(a, b) if (a == b) 0 else 7
  expect_equal(node_distance(hw_net(transition_cost(penalty = f, connectors = "motorway_link")), "a0", "l1"), 307)
  bad <- hw_net(transition_cost(penalty = function(a, b) -1, connectors = "motorway_link"))
  expect_error(node_distance(bad, "a0", "l1"), "invalid cost")
})

test_that("transition costs without penalties match edge costs", {
  net_t <- hw_net(transition_cost("length"))
  net_e <- hw_net("length")
  expect_equal(distances_from(net_t, 1:6), distances_from(net_e, 1:6))
})

test_that("nearest source matches brute force", {
  net <- hw_net(transition_cost("length", default_penalty = 50, connectors = "motorway_link"))
  src <- match(c("a0", "l1"), net$nodes$name)
  ns <- nearest_source(net, src)
  full <- distances_from(net, src)
  expect_equal(ns$dist, apply(full, 2, min))
  expect_equal(ns$src, src[apply(full, 2, which.min)])
})

test_that("node contrast metric", {
  nodes <- data.frame(name = 0:1, x = c(0, 100), y = 0, landuse = c(0, 3))
  net <- network_graph(data.frame(from = 0, to = 1, length = 100), nodes, node_contrast("landuse", weight = 10))
  expect_equal(node_distance(net, 0, 1), 130)
})

test_that("precomputed metric", {
  h <- highway_fixture()
  D <- abs(outer(1:6, 1:6, "-"))
  net <- network_graph(h$edges, h$nodes, precomputed(D, nodes = h$nodes$name))
  expect_equal(node_distance(net, "a0", "b1"), 3)
  expect_error(precomputed(matrix(-1, 2, 2)), "non-negative")
  expect_error(precomputed(matrix(1, 2, 3)), "square")
  expect_error(node_distance(network_graph(h$edges, h$nodes, precomputed(diag(2))), "a0", "b1"), "nodes")
})

test_that("negative costs are rejected", {
  expect_error(node_distance(hw_net(function(edges) rep(-1, nrow(edges))), "a0", "b0"), "invalid cost")
})

test_that("infinite costs block an edge", {
  net <- hw_net(function(edges) ifelse(edges$highway == "residential", Inf, edges$length))
  expect_equal(node_distance(net, "a0", "l1"), Inf)
})

test_that("zero-cost edges are kept", {
  nodes <- data.frame(name = 0:2, x = 0:2, y = 0)
  net <- network_graph(data.frame(from = 0:1, to = 1:2, length = 0), nodes)
  expect_equal(node_distance(net, 0, 2), 0)
})

test_that("direction of directed graphs", {
  nodes <- data.frame(name = 0:1, x = 0:1, y = 0)
  e <- data.frame(from = 0, to = 1, length = 5)
  expect_equal(node_distance(network_graph(e, nodes, directed = TRUE, direction = "out"), 1, 0), Inf)
  expect_equal(node_distance(network_graph(e, nodes, directed = TRUE, direction = "in"), 1, 0), 5)
  expect_equal(node_distance(network_graph(e, nodes, directed = TRUE, direction = "both"), 1, 0), 5)
  expect_error(network_graph(e, nodes, direction = "sideways"), "direction")
})

test_that("metric specifications are coerced", {
  expect_s3_class(rd2dnet:::as_metric("length"), "rd2dnet_edge_cost")
  expect_s3_class(rd2dnet:::as_metric(matrix(0, 2, 2)), "rd2dnet_precomputed")
  expect_error(rd2dnet:::as_metric(3), "network metric")
})

test_that("bounded searches find exactly the pairs within the limit", {
  city <- mock_city()
  nets <- list(
    edge = city$network("length"),
    transition = city$network(transition_cost("length", default_penalty = 120, connectors = "motorway_link"))
  )
  nets$precomputed <- city$network(precomputed(distances_from(nets$edge, seq_len(450))))
  for (net in nets) {
    src <- seq(1, net$n_nodes, by = 7)
    full <- distances_from(net, src)
    pr <- rd2dnet:::.engine_pairs(rd2dnet:::.engine(net), src, 400)
    got <- matrix(FALSE, nrow(full), ncol(full))
    got[cbind(pr$i, pr$node)] <- TRUE
    expect_equal(got, full < 400)
    expect_equal(pr$dist, full[cbind(pr$i, pr$node)])
  }
})

small_net <- function(pos, pairs, attrs = NULL, ...) {
  nodes <- data.frame(name = names(pos), x = vapply(pos, `[`, 0, 1), y = vapply(pos, `[`, 0, 2))
  from <- vapply(pairs, `[`, "", 1)
  to <- vapply(pairs, `[`, "", 2)
  len <- sqrt((nodes$x[match(from, nodes$name)] - nodes$x[match(to, nodes$name)])^2 +
                (nodes$y[match(from, nodes$name)] - nodes$y[match(to, nodes$name)])^2)
  edges <- data.frame(from = from, to = to, length = len)
  if (!is.null(attrs)) edges <- cbind(edges, attrs)
  network_graph(edges, nodes, ...)
}

test_that("points on one-way edges are reached in the edge's direction", {
  # X is on the boundary; the one-way edge N -> M can only be entered at N (via X -> R -> N)
  pos <- list(X = c(0, 0), M = c(0, -100), N = c(0, -200), R = c(300, -100), T = c(0, 100),
              L = c(-100, 0), Rb = c(100, 0))
  pairs <- list(c("X", "M"), c("N", "M"), c("X", "R"), c("R", "N"), c("X", "T"), c("L", "X"), c("X", "Rb"))
  net <- small_net(pos, pairs, directed = TRUE, direction = "out")
  bnd <- boundary(rbind(c(-100, 0), c(100, 0)), treated_side = "left")
  cc <- network_coordinates(net, bnd, rbind(c(0, -110), c(0, -150)), eval_points = 100)
  want <- 2 * sqrt(300^2 + 100^2) + c(90, 50)
  expect_equal(cc$s, -want, tolerance = 1e-12)
  expect_equal(coords_at(cc, 1)$d, want, tolerance = 1e-12)
})

test_that("the transition penalty applies to the edge a point sits on", {
  pos <- list(P = c(0, -200), Q = c(0, -100), X = c(0, 0), Y = c(0, 100), S = c(300, -100),
              Z = c(300, 0), W = c(300, 100))
  pairs <- list(c("P", "Q"), c("Q", "X"), c("X", "Y"), c("Q", "S"), c("S", "Z"), c("Z", "W"), c("X", "Z"))
  cls <- data.frame(highway = c("primary", "residential", "residential", "primary", "primary", "primary",
                                "primary"))
  m <- transition_cost("length", class_attr = "highway", default_penalty = 1000)
  net <- small_net(pos, pairs, cls, metric = m)
  bnd <- boundary(rbind(c(-50, 0), c(350, 0)), treated_side = "left")
  cc <- network_coordinates(net, bnd, rbind(c(0, -150), c(0, -199.9), c(0, -200)), eval_points = 50)
  expect_equal(cc$s, c(-450, -499.9, -500), tolerance = 1e-12)
})

test_that("node contrasts add up across an edge split at the boundary", {
  nodes <- data.frame(name = c("a", "b", "c", "d"), x = c(0, 0, 200, 200), y = c(-100, 100, -100, 100),
                      far = c(1, 5, 1, 5))
  edges <- data.frame(from = c("a", "c", "a", "b"), to = c("b", "d", "c", "d"), length = 200)
  net <- network_graph(edges, nodes, node_contrast("far", 10))
  bnd <- boundary(rbind(c(-50, 0), c(250, 0)), treated_side = "right")
  split <- with_boundary_crossings(net, bnd)
  expect_equal(node_distance(split, "a", "b"), 240)
})

test_that("reusing a split network gives the same coordinates", {
  # coords$network carries inserted nodes; splitting it again must not reuse their labels
  g <- expand.grid(i = 0:5, j = 0:5)
  nodes <- data.frame(name = paste(g$i, g$j), x = 100 * g$i, y = 100 * g$j)
  e <- list()
  for (r in seq_len(nrow(g))) {
    i <- g$i[r]
    j <- g$j[r]
    for (dd in list(c(1, 0), c(0, 1))) {
      if (i + dd[1] <= 5 && j + dd[2] <= 5) {
        a <- paste(i, j)
        b <- paste(i + dd[1], j + dd[2])
        e[[length(e) + 1]] <- data.frame(from = a, to = b, length = 100)
        if ((i + j) %% 3 != 0) e[[length(e) + 1]] <- data.frame(from = b, to = a, length = 100)
      }
    }
  }
  net <- network_graph(do.call(rbind, e), nodes, directed = TRUE, direction = "out")
  bnd <- boundary(rbind(c(-50, 250), c(550, 250)), treated_side = "right")
  set.seed(0)
  pts <- matrix(runif(600, 0, 500), ncol = 2)
  cc <- network_coordinates(net, bnd, pts, n_eval = 5)
  cc2 <- network_coordinates(cc$network, bnd, pts, n_eval = 5)
  expect_equal(cc2$s, cc$s, tolerance = 1e-12)
})
