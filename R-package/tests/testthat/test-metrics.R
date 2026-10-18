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
