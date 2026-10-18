test_that("treated side left and right", {
  line <- rbind(c(0, 0), c(10, 0))
  pts <- rbind(c(5, 1), c(5, -1))
  expect_equal(boundary_is_treated(boundary(line, treated_side = "left"), pts), c(TRUE, FALSE))
  expect_equal(boundary_is_treated(boundary(line, treated_side = "right"), pts), c(FALSE, TRUE))
})

test_that("boundary requires exactly one side specification", {
  expect_error(boundary(rbind(c(0, 0), c(1, 0))), "exactly one")
  expect_error(boundary(rbind(c(0, 0), c(1, 0)), treated_side = "up"), "left")
})

test_that("masks by geometry and by arc interval", {
  b <- boundary(rbind(c(0, 0), c(100, 0)), treated_side = "left",
                mask = list(rect_polygon(20, 0, 40, 10), c(80, 90)))
  expect_equal(boundary_is_masked(b, c(10, 30, 50, 85)), c(FALSE, TRUE, FALSE, TRUE))
  expect_false(boundary_is_masked(b, NA))
})

test_that("boundary geometry helpers", {
  b <- boundary(rbind(c(0, 0), c(100, 0), c(100, 50)), treated_side = "left")
  expect_equal(b$length, 150)
  expect_equal(boundary_project(b, rbind(c(10, 5), c(120, 20))), c(10, 120))
  expect_equal(boundary_distance(b, rbind(c(10, 5), c(120, 20))), c(5, 20))
  expect_equal(boundary_interpolate(b, c(50, 125, 999)), rbind(c(50, 0), c(100, 25), c(100, 50)))
  expect_equal(boundary_grid(b, n = 4), c(0, 50, 100, 150))
  expect_equal(boundary_grid(b, spacing = 60), c(0, 60, 120))
  expect_output(print(b), "length=150")
})

test_that("treated regions accept holes and several rings", {
  outer <- rect_polygon(0, 0, 10, 10)
  hole <- rect_polygon(4, 4, 6, 6)
  b <- boundary(rbind(c(0, 0), c(10, 0)), treated_region = list(outer, hole))
  expect_equal(boundary_is_treated(b, rbind(c(1, 1), c(5, 5), c(20, 20))), c(TRUE, FALSE, FALSE))
})

test_that("coordinates on a grid have the right signs and values", {
  g <- grid_fixture()
  net <- network_graph(g$edges, g$nodes)
  pts <- rbind(c(100, 50), c(100, 350), c(300, 120))
  cc <- network_coordinates(net, g$boundary, pts, eval_points = 150) # eval at x = 100
  expect_equal(cc$treated, c(TRUE, FALSE, TRUE))
  # Points on the avenue x = 100: straight up/down the avenue to y = 200.
  expect_equal(cc$s[1:2], c(150, -150))
  # (300, 120) snaps onto the x = 300 avenue: 80 m to the boundary, no depth.
  expect_equal(cc$s[3], 80)
  at <- coords_at(cc, 1)
  expect_equal(at$o[1:2], c(0, 0))
  expect_equal(at$o[3], 200) # crossing at x = 300 is 200 m along the boundary
  expect_equal(at$d[3], 280)
  expect_output(print(cc), "eval_points=1")
})

test_that("evaluation points between nodes", {
  g <- grid_fixture()
  cc <- network_coordinates(network_graph(g$edges, g$nodes), g$boundary, rbind(c(100, 50)), eval_points = 200)
  at <- coords_at(cc, 1)
  expect_equal(cc$eval_node_arc, 200)
  expect_equal(at$d, 50 + 150)
  expect_equal(at$o, -50)
  et <- eval_table(cc)
  expect_named(et, c("arc", "node", "node_arc", "x", "y", "masked"))
  expect_equal(c(et$x, et$y), c(150, 200))
})

test_that("evaluation points can be given as coordinates", {
  g <- grid_fixture()
  net <- network_graph(g$edges, g$nodes)
  cc <- network_coordinates(net, g$boundary, rbind(c(100, 50)), eval_points = rbind(c(150, 210)))
  expect_equal(cc$eval_arc, 200)
})

test_that("the metric changes the coordinates", {
  g <- grid_fixture()
  g$edges$travel_time <- g$edges$length / 10
  cc <- network_coordinates(network_graph(g$edges, g$nodes, "travel_time"), g$boundary, rbind(c(100, 0)))
  expect_equal(cc$s, 20) # 200 m at 10 m/s, snapping depth costs 0
})

test_that("unreachable points are dropped with a warning", {
  g <- grid_fixture()
  nodes <- rbind(g$nodes, data.frame(name = c(100L, 101L), x = c(1000, 1000), y = c(0, 50)))
  edges <- rbind(g$edges, data.frame(from = 100L, to = 101L, length = 50))
  expect_warning(
    cc <- network_coordinates(network_graph(edges, nodes), g$boundary, rbind(c(1000, 10), c(0, 0))),
    "cannot reach"
  )
  expect_true(is.na(cc$s[1]))
  expect_true(is.finite(cc$s[2]))
  expect_true(is.na(coords_at(cc, 1)$o[1]))
})

test_that("mock city coordinates", {
  city <- mock_city()
  cc <- mock_coords()
  expect_equal(cc$n_obs, nrow(city$obs))
  expect_true(all((cc$s >= 0) == cc$treated))
  # On a grid with Euclidean edge snapping, |s| is at least the straight-line distance.
  expect_true(all(abs(cc$s) >= abs(city$obs$y - 1600) - 1e-6))
  expect_gt(sum(cc$eval_masked), 0)
  df <- as.data.frame(cc, j = 3)
  expect_named(df, c("x", "y", "treated", "s", "crossing_arc", "depth", "d", "o"))
})

test_that("pairwise distances are symmetric", {
  D <- pairwise_distances(mock_coords(), 1:25)
  expect_equal(D, t(D))
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0))
})

test_that("distance along a shared segment is direct", {
  g <- grid_fixture()
  cc <- network_coordinates(network_graph(g$edges, g$nodes), g$boundary,
                            rbind(c(120, 210), c(120, 190)), eval_points = 200)
  at <- coords_at(cc, 1) # evaluation point at x = 150 on the boundary street
  expect_equal(at$d, c(10 + 30, 10 + 30))
  expect_equal(at$o, c(-30, -30))
  expect_equal(abs(cc$s), c(10, 10))
})

test_that("repeated locations share coordinates", {
  city <- mock_city()
  xy <- as.matrix(city$obs[1:400, c("x", "y")])
  set.seed(1)
  rep_xy <- xy[rep(1:400, each = 3), ][sample(1200), ]
  a <- network_coordinates(city$network(), city$boundary, xy, n_eval = 9)
  b <- network_coordinates(city$network(), city$boundary, rep_xy, n_eval = 9)
  expect_equal(b$n_locations, 400)
  expect_equal(b$n_obs, 1200)
  m <- match(paste(rep_xy[, 1], rep_xy[, 2]), paste(xy[, 1], xy[, 2]))
  expect_identical(b$s, a$s[m])
  for (j in 1:9) {
    expect_identical(coords_at(b, j)$d, coords_at(a, j)$d[m])
    expect_identical(coords_at(b, j)$o, coords_at(a, j)$o[m])
  }
})

test_that("sparse network HAC matches the dense kernel, with repeats", {
  city <- mock_city()
  xy <- as.matrix(city$obs[rep(1:250, each = 2), c("x", "y")])
  cc <- network_coordinates(city$network(), city$boundary, xy, n_eval = 5)
  set.seed(1)
  Psi <- matrix(rnorm(500 * 3), 500, 3)
  for (cutoff in c(120, 400)) {
    expect_equal(rd2dnet:::.hac_quadratic(cc, 1:500, Psi, cutoff), dense_hac(cc, 1:500, Psi, cutoff),
                 tolerance = 1e-10)
  }
})

test_that("sparse HAC on a directed network and with Euclidean distance", {
  g <- grid_fixture()
  drop <- g$edges$from %% 3 == 0 & g$edges$to %% 2 == 0 # some one-way streets
  edges <- rbind(g$edges, g$edges[!drop, ] |> transform(from = to, to = from))
  net <- network_graph(edges, g$nodes, directed = TRUE, direction = "out")
  set.seed(2)
  pts <- matrix(runif(240, 0, 400), 120, 2)
  cc <- network_coordinates(net, g$boundary, pts, n_eval = 3)
  idx <- which(is.finite(cc$s))
  Psi <- matrix(rnorm(length(idx) * 2), length(idx), 2)
  expect_equal(rd2dnet:::.hac_quadratic(cc, idx, Psi, 150), dense_hac(cc, idx, Psi, 150), tolerance = 1e-10)
  D <- as.matrix(stats::dist(cc$xy[idx, ]))
  K <- pmax(1 - D / 150, 0)
  expect_equal(rd2dnet:::.hac_quadratic(cc, idx, Psi, 150, "euclidean"), t(Psi) %*% K %*% Psi,
               tolerance = 1e-10)
})

test_that("points can be given as sf objects", {
  skip_if_not_installed("sf")
  g <- grid_fixture()
  net <- network_graph(g$edges, g$nodes)
  pts <- sf::st_as_sf(data.frame(x = c(100, 300), y = c(50, 120)), coords = c("x", "y"))
  line <- sf::st_linestring(rbind(c(-50, 200), c(450, 200)))
  region <- sf::st_polygon(list(rect_polygon(-1e3, -1e3, 1e3, 200)))
  b <- boundary(line, treated_region = region)
  cc <- network_coordinates(net, b, pts, eval_points = 150)
  expect_equal(cc$s, c(150, 80))
})
