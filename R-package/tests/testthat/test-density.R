test_that("no sorting in the mock city", {
  r <- density(mock_coords())
  expect_s3_class(r, "rdnet_density")
  expect_equal(r$normalize, "network")
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
  r <- density(cc)
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
