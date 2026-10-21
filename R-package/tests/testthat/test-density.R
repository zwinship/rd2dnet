test_that("no sorting in the mock city", {
  r <- density(mock_coords())
  expect_s3_class(r, "rdnet_density")
  expect_equal(r$normalize, "euclidean")
  expect_true(is.finite(r$diff_bc) && r$se > 0)
  expect_gt(r$pvalue, 0.01)
  expect_output(print(r), "density test")
})

test_that("street-length normalisation", {
  r <- density(mock_coords(), normalize = "network")
  expect_true(is.finite(r$diff_bc))
  expect_gt(r$pvalue, 0.01)
  expect_lt(abs(r$theta), 0.5) # uniform parcels per metre of street on both sides
  expect_error(density(mock_coords(), normalize = "bogus"), "normalize")
})

test_that("the jackknife standard error tracks the Monte Carlo spread", {
  skip_on_cran()
  set.seed(0)
  tri <- get_kernel("triangular")
  sims <- replicate(200, {
    x <- runif(2000, -1, 1)
    f <- rd2dnet:::.cjm(x, rep(1, 2000), 0.4, 0.4, 2, tri)
    c(f$b[4], sqrt(f$V[4, 4]))
  })
  expect_equal(mean(sims[1, ]), 0.5, tolerance = 0.04)
  expect_equal(mean(sims[2, ]), sd(sims[1, ]), tolerance = 0.2)
})

test_that("the test detects bunching", {
  city <- make_mock_midtown(n_obs = 3000, seed = 5)
  obs <- city$obs
  set.seed(1)
  extra <- as.matrix(obs[sample(nrow(obs), 1500), c("x", "y")])
  extra[, 2] <- 1600 - runif(1500, 1, 60) # pile up just inside the treated zone
  cc <- network_coordinates(city$network(), city$boundary, rbind(as.matrix(obs[, c("x", "y")]), extra), n_eval = 5)
  for (norm in c("network", "none")) {
    r <- density(cc, normalize = norm)
    expect_gt(r$theta, 0)
    expect_lt(r$pvalue, 0.01)
  }
})

test_that("masked stretches are excluded", {
  cc <- mock_coords()
  r <- density(cc, normalize = "network")
  masked <- boundary_is_masked(cc$boundary, cc$crossing_arc)
  expect_true(any(masked))
  expect_equal(r$n_left + r$n_right, sum(!masked & is.finite(cc$s)))
})

test_that("bandwidth choices", {
  a <- density(mock_coords(), bwselect = "common")
  b <- density(mock_coords(), bwselect = "each")
  expect_equal(a$h_left, a$h_right)
  expect_false(b$h_left == b$h_right)
  expect_error(density(mock_coords(), bwselect = "bogus"), "bwselect")
  expect_error(density(mock_coords(), p = 0), "p <= q")
})

test_that("a local window", {
  cc <- mock_coords()
  j <- which(!cc$eval_masked)[3]
  r <- density(cc, eval_index = j, window = 600, h = 400)
  expect_lt(r$n_left + r$n_right, cc$n_obs)
  expect_error(density(cc, eval_index = j), "window")
})

grid_net <- function(n = 11, step = 100) {
  g <- expand.grid(i = 0:(n - 1), j = 0:(n - 1))
  nodes <- data.frame(name = paste(g$i, g$j), x = step * g$i, y = step * g$j)
  e <- rbind(data.frame(from = paste(g$i, g$j), to = paste(g$i + 1, g$j))[g$i < n - 1, ],
             data.frame(from = paste(g$i, g$j), to = paste(g$i, g$j + 1))[g$j < n - 1, ])
  e$length <- step
  network_graph(e, nodes)
}

test_that("the density window wraps around a closed ring", {
  a <- seq(0, 2 * pi, length.out = 257)
  ring <- boundary_from_region(cbind(500 + 300 * cos(a), 500 + 300 * sin(a)))
  set.seed(8)
  X <- matrix(stats::runif(16000, 0, 1000), ncol = 2)
  cc <- network_coordinates(grid_net(), ring, X, n_eval = 8)
  n <- vapply(1:8, function(j) density(cc, eval_index = j, window = 150, h = 40)$n_left, numeric(1))
  expect_gt(n[1], 0.7 * stats::median(n))
})

test_that("sorting tests warn for points on or near streets", {
  set.seed(3)
  t <- stats::runif(3000, 0, 1000)
  k <- sample(0:10, 3000, replace = TRUE) * 100
  horiz <- stats::runif(3000) < 0.5
  on <- cbind(ifelse(horiz, t, k), ifelse(horiz, k, t))
  bnd <- boundary(rbind(c(550, -10), c(550, 1010)), treated_side = "right")
  for (noise in c(0, 5)) {
    cc <- network_coordinates(grid_net(), bnd, on + stats::rnorm(6000, 0, noise), n_eval = 5)
    expect_warning(density(cc), "next to the network|network itself")
  }
  cc <- network_coordinates(grid_net(), bnd, matrix(stats::runif(6000, 0, 1000), ncol = 2), n_eval = 5)
  expect_no_warning(density(cc))
})
