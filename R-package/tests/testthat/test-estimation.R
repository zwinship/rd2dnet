y_mock <- function() mock_city()$obs$outcome

test_that("polynomial basis terms", {
  e <- rd2dnet:::.poly_exponents(1)
  expect_equal(unname(e), cbind(c(0, 1, 0), c(0, 0, 1)))
  expect_equal(nrow(rd2dnet:::.poly_exponents(2)), 6)
  expect_equal(c(rd2dnet:::.poly_basis(2, 3, 2)), c(1, 2, 3, 4, 6, 9))
})

test_that("weighted least squares recovers a linear function", {
  set.seed(0)
  s <- rnorm(200)
  o <- rnorm(200)
  f <- rd2dnet:::.wls(rd2dnet:::.poly_basis(s, o, 1), 1 + 2 * s - 3 * o, runif(200, 0.1, 1))
  expect_equal(f$beta, c(1, 2, -3), tolerance = 1e-10)
})

test_that("kernels", {
  for (name in c("triangular", "epanechnikov", "uniform", "gaussian", "tri", "epa", "uni", "gau")) {
    k <- get_kernel(name)
    expect_gt(k(0), 0)
    expect_equal(k(Inf), 0)
    if (!name %in% c("gaussian", "gau")) expect_equal(k(1.5), 0)
  }
  expect_identical(get_kernel(dnorm), dnorm)
  expect_error(get_kernel("bogus"), "unknown kernel")
})

test_that("the generator builds the bundled city's design", {
  a <- load_mock_midtown()
  b <- make_mock_midtown(n_obs = 500, seed = 4)
  expect_equal(nrow(a$edges), nrow(b$edges))
  expect_equal(a$nodes[, c("name", "x", "y")], b$nodes[, c("name", "x", "y")])
  expect_equal(nrow(b$obs), 500)
  expect_named(b$obs, names(a$obs))
  expect_equal(b$obs$treated, b$obs$y < 1600)
  # the global random number stream is left untouched
  set.seed(9)
  r1 <- runif(1)
  set.seed(9)
  invisible(make_mock_midtown(n_obs = 10, seed = 1))
  expect_equal(runif(1), r1)
  expect_output(print(a), "450 nodes")
  expect_null(make_mock_midtown(n_obs = 10, park = FALSE)$park)
})

test_that("rdnet runs and masks", {
  cc <- mock_coords()
  res <- rdnet(y_mock(), cc, seed = 0)
  t <- res$table
  expect_equal(nrow(t), cc$n_eval)
  expect_true(all(t$status[cc$eval_masked] == "masked"))
  ok <- t$status == "ok"
  expect_gte(sum(ok), 8)
  expect_true(all(t$se[ok] > 0))
  expect_true(all(t$cb_lo[ok] <= t$ci_lo[ok] + 1e-12))
  expect_output(print(res), "rd2dnet")
  expect_s3_class(res$bandwidths, "data.frame")
})

test_that("rdnet recovers a large constant effect", {
  city <- make_mock_midtown(n_obs = 6000, seed = 3)
  obs <- city$obs
  set.seed(1)
  y <- 0.0008 * (obs$y - 1600) + 3 * obs$treated + rnorm(nrow(obs), 0, 0.1)
  cc <- network_coordinates(city$network(), city$boundary, obs, n_eval = 15)
  res <- rdnet(y, cc, h = 500, cbands = FALSE)
  ok <- res$table$status == "ok"
  expect_gte(sum(ok), 5)
  expect_equal(res$table$tau[ok], rep(3, sum(ok)), tolerance = 0.05)
})

test_that("pooled option runs", {
  res <- rdnet(y_mock(), mock_coords(), h = 600, pooled = TRUE, cbands = FALSE)
  expect_gte(sum(res$table$status == "ok"), 8)
})

test_that("vce variants and their ordering", {
  se <- list()
  for (v in c("hc0", "hc1", "hc2", "hc3")) {
    t <- rdnet(y_mock(), mock_coords(), h = 600, vce = v, cbands = FALSE)$table
    expect_true(all(is.finite(t$se[t$status == "ok"])))
    se[[v]] <- t$se
  }
  ok <- !is.na(se$hc0)
  expect_true(all(se$hc1[ok] >= se$hc0[ok]))
  expect_true(all(se$hc3[ok] >= se$hc0[ok]))
})

test_that("cluster and Conley standard errors", {
  y <- y_mock()
  cc <- mock_coords()
  rc <- rdnet(y, cc, h = 600, cluster = mock_city()$obs$block, seed = 0)
  expect_equal(rc$vce, "cluster")
  expect_gte(sum(rc$table$status == "ok"), 8)
  rn <- rdnet(y, cc, h = 600, conley_cutoff = 150, cbands = FALSE)
  re <- rdnet(y, cc, h = 600, conley_cutoff = 150, conley_distance = "euclidean", cbands = FALSE)
  expect_equal(rn$vce, "conley")
  ok <- rn$table$status == "ok"
  expect_true(all(is.finite(rn$table$se[ok])) && all(is.finite(re$table$se[ok])))
})

test_that("intensity weights change the estimates", {
  y <- y_mock()
  a <- rdnet(y, mock_coords(), h = 600, cbands = FALSE)$table$tau
  b <- rdnet(y, mock_coords(), h = 600, weights = mock_city()$obs$far, cbands = FALSE)$table$tau
  ok <- !is.na(a) & !is.na(b)
  expect_false(isTRUE(all.equal(a[ok], b[ok])))
})

test_that("bandwidth selection", {
  bw <- rdbwnet(y_mock(), mock_coords())
  ok <- bw$status == "ok"
  expect_gte(sum(ok), 8)
  expect_equal(bw$h0[ok], bw$h1[ok])
  expect_true(all(bw$h0[ok] >= pmax(bw$hmin0[ok], bw$hmin1[ok]) - 1e-9))
  expect_true(all(is.finite(as.matrix(bw[ok, c("V0", "V1", "B0", "B1", "Reg10", "Reg20", "bn0", "bn1", "dn")]))))
  expect_error(rdbwnet(y_mock(), mock_coords(), bwselect = "bogus"), "bwselect")
})

test_that("bandwidth selector variants", {
  for (sel in c("msetwo", "imserd", "imsetwo")) {
    bw <- rdbwnet(y_mock(), mock_coords(), bwselect = sel)
    ok <- bw$status == "ok"
    expect_true(all(bw$h0[ok] > 0 & bw$h1[ok] > 0))
    if (sel == "msetwo") expect_false(isTRUE(all.equal(bw$h0[ok], bw$h1[ok])))
    if (startsWith(sel, "i")) {
      inner <- ok & bw$h0 > bw$hmin0 & bw$h1 > bw$hmin1
      expect_lte(length(unique(round(bw$h0[inner], 8))), 2)
    }
  }
})

test_that("side-specific bandwidths in estimation", {
  res <- rdnet(y_mock(), mock_coords(), bwselect = "msetwo", cbands = FALSE)
  ok <- res$table$status == "ok"
  expect_gte(sum(ok), 6)
  expect_true(any(is.na(res$table$h[ok])))
  res2 <- rdnet(y_mock(), mock_coords(), h = c(400, 700), cbands = FALSE)
  expect_true(all(res2$table$h1[res2$table$status == "ok"] == 700))
  expect_error(rdnet(y_mock(), mock_coords(), h = 1:3), "h must be")
})

test_that("bandwidths with clusters and a fuzzy design", {
  y <- y_mock()
  cc <- mock_coords()
  set.seed(0)
  D <- as.numeric(ifelse(cc$treated, runif(length(y)) < 0.8, runif(length(y)) < 0.2))
  bw <- rdbwnet(y, cc, cluster = mock_city()$obs$block, fuzzy = D)
  expect_true(all(bw$h0[bw$status == "ok"] > 0))
})

test_that("fuzzy design recovers the local effect", {
  city <- make_mock_midtown(n_obs = 6000, seed = 11)
  obs <- city$obs
  set.seed(2)
  T <- obs$treated
  D <- as.numeric(ifelse(T, runif(nrow(obs)) < 0.75, runif(nrow(obs)) < 0.25))
  y <- 0.0008 * (obs$y - 1600) + 2 * D + rnorm(nrow(obs), 0, 0.1)
  cc <- network_coordinates(city$network(), city$boundary, obs, n_eval = 9)
  res <- rdnet(y, cc, fuzzy = D, h = 600, seed = 0)
  t <- res$table[res$table$status == "ok", ]
  expect_true(res$fuzzy)
  expect_gte(nrow(t), 4)
  expect_true(all(abs(t$fs - 0.5) <= 4 * t$fs_se))
  expect_true(all(abs(t$tau - 2) <= 4 * t$se))
  expect_gte(mean(t$ci_lo <= 2 & 2 <= t$ci_hi), 0.75)
  expect_error(rdnet(y, cc, fuzzy = D[-1]), "fuzzy")
  expect_output(print(res), "fuzzy")
})

test_that("a vanishing Conley cutoff reproduces HC0, band included", {
  y <- y_mock()
  hc <- rdnet(y, mock_coords(), h = 600, vce = "hc0", seed = 0)
  cl <- rdnet(y, mock_coords(), h = 600, conley_cutoff = 1e-9, conley_resid_scale = 1, seed = 0)
  expect_equal(cl$table$se_rb, hc$table$se_rb, tolerance = 1e-8)
  expect_equal(cl$cov, hc$cov, tolerance = 1e-8)
  wide <- rdnet(y, mock_coords(), h = 600, conley_cutoff = 300, seed = 0)
  expect_false(is.null(wide$crit))
  expect_true(any(is.finite(wide$table$cb_lo)))
})

test_that("sparse Conley kernel matches the dense one on the mock city", {
  set.seed(0)
  Psi <- matrix(rnorm(900), 300, 3)
  K <- pmax(1 - pairwise_distances(mock_coords(), 1:300, cutoff = 200) / 200, 0)
  expect_equal(rd2dnet:::.hac_quadratic(mock_coords(), 1:300, Psi, 200), t(Psi) %*% K %*% Psi, tolerance = 1e-8)
})

test_that("input validation", {
  y <- y_mock()
  cc <- mock_coords()
  expect_error(rdnet(y[-1], cc), "one value per observation")
  expect_error(rdnet(y, cc, vce = "cluster"), "cluster variable")
  expect_error(rdnet(y, cc, weights = -rep(1, length(y))), "non-negative")
  expect_error(rdnet(y, cc, vce = "bogus"), "vce")
  expect_error(rdnet(y, cc, vce = "conley"), "conley_cutoff")
  expect_error(rdnet(y, cc, conley_cutoff = 100, conley_distance = "bogus"), "conley_distance")
  expect_error(rdnet(y, "not coords"), "network_coordinates")
})

test_that("travel-time metric end to end", {
  city <- mock_city()
  cc <- network_coordinates(city$network(travel_time()), city$boundary, city$obs, n_eval = 11)
  res <- rdnet(city$obs$outcome, cc, cbands = FALSE)
  expect_gte(sum(res$table$status == "ok"), 3)
  # bandwidths are now in seconds, much smaller than metres
  expect_lt(max(res$table$h, na.rm = TRUE), 300)
})

test_that("plot and methods", {
  res <- rdnet(y_mock(), mock_coords(), h = 600, seed = 0)
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(res, true_tau = mock_true_tau(res$table$node_arc)))
  expect_length(coef(res), nrow(res$table))
  expect_equal(unname(coef(res, "tau")), res$table$tau)
  expect_equal(unname(diag(vcov(res))), res$table$se_rb^2)
  ci <- confint(res)
  ok <- res$table$status == "ok"
  expect_equal(unname(ci[ok, 1]), res$table$ci_lo[ok])
  expect_identical(as.data.frame(res), res$table)
})

test_that("covariate adjustment matches rd2d's plug-in and shrinks standard errors", {
  bc <- big_coords()
  cc <- bc$coords
  set.seed(0)
  X <- rnorm(cc$n_obs)
  y <- 0.0005 * cc$s + 2 * cc$treated + 3 * X + rnorm(cc$n_obs)
  plain <- rdnet(y, cc, h = 500, cbands = FALSE)$table
  adj <- rdnet(y, cc, h = 500, covs = X, cbands = FALSE)$table
  expect_true(all(adj$se < 0.6 * plain$se))
  # rd2d: common gamma from partialled-out side fits, then side-specific fits on y - X gamma
  at <- coords_at(cc, 1)
  w <- get_kernel("triangular")(at$d / 500)
  use <- w > 0 & is.finite(at$d)
  sides <- list(which(use & cc$treated), which(use & !cc$treated))
  A <- 0
  B <- 0
  for (i in sides) {
    R <- rd2dnet:::.poly_basis(cc$s[i] / 500, at$o[i] / 500, 1)
    G <- solve(crossprod(R * w[i], R))
    uz <- crossprod(R * w[i], X[i])
    A <- A + sum(X[i] * w[i] * X[i]) - drop(t(uz) %*% G %*% uz)
    B <- B + sum(X[i] * w[i] * y[i]) - drop(t(uz) %*% G %*% crossprod(R * w[i], y[i]))
  }
  ya <- y - (B / A) * X
  fit <- function(i) rd2dnet:::.wls(rd2dnet:::.poly_basis(cc$s[i] / 500, at$o[i] / 500, 1), ya[i], w[i])$beta[1]
  expect_equal(adj$tau[1], fit(sides[[1]]) - fit(sides[[2]]), tolerance = 1e-8)
})

test_that("constant covariates are dropped", {
  y <- y_mock()
  a <- rdnet(y, mock_coords(), h = 600, cbands = FALSE)$table
  b <- rdnet(y, mock_coords(), h = 600, covs = rep(1, length(y)), cbands = FALSE)$table
  expect_equal(a$tau, b$tau)
})

test_that("derivative targets are exact without noise", {
  cc <- big_coords()$coords
  s <- cc$s
  T <- as.numeric(cc$treated)
  set.seed(0)
  y <- 0.001 * s + 0.003 * s * T + 1e-9 * rnorm(cc$n_obs)
  k <- rdnet(y, cc, h = 500, p = 2, deriv = c(1, 0), cbands = FALSE)$table
  expect_equal(k$tau, rep(0.003, 4), tolerance = 1e-4)
  expect_equal(k$tau_bc, rep(0.003, 4), tolerance = 1e-4)
  tv <- rdnet(y, cc, h = 500, p = 2, tangvec = c(1, 0), cbands = FALSE)$table
  expect_equal(tv$tau, k$tau)
  # a jump that grows along the boundary: d tau / d o = 0.002. Parcels that reach
  # the boundary at the highway node (via a ramp, not along 60th St) are left out.
  arc <- cc$crossing_arc
  y2 <- ifelse(arc >= 399, (0.5 + 0.002 * arc) * T, NA)
  set.seed(1)
  g <- rdnet(y2 + 1e-9 * rnorm(cc$n_obs), cc, h = 600, p = 1, deriv = c(0, 1), cbands = FALSE)$table
  expect_equal(g$tau, rep(0.002, 4), tolerance = 1e-3)
})

test_that("derivative input validation", {
  y <- y_mock()
  cc <- mock_coords()
  expect_error(rdnet(y, cc, p = 1, deriv = c(2, 0)), "deriv")
  expect_error(rdnet(y, cc, deriv = c(1, 0), tangvec = c(1, 0)), "either deriv or tangvec")
  expect_error(rdnet(y, cc, pooled = TRUE, deriv = c(1, 0)), "level jumps")
  expect_error(rdnet(y, cc, p = 0, tangvec = c(1, 0)), "p >= 1")
})

test_that("coverage-error-optimal bandwidths shrink the MSE ones", {
  y <- y_mock()
  cc <- mock_coords()
  for (pair in list(c("cerrd", "mserd"), c("certwo", "msetwo"))) {
    m <- rdbwnet(y, cc, bwselect = pair[2], bwcheck = 1)
    c2 <- rdbwnet(y, cc, bwselect = pair[1], bwcheck = 1)
    ok <- m$status == "ok" & c2$h0 > c2$hmin0 + 1e-9
    expect_true(any(ok))
    ratio <- c2$h0[ok] / m$h0[ok]
    expect_true(all(ratio < 1))
    if (pair[1] == "cerrd") {
      n <- nrow(unique(cc$xy))
      expect_equal(ratio, rep(n^(1 / 6 - 1 / 5), length(ratio)), tolerance = 1e-8)
    }
  }
})

test_that("bandwidths with derivatives and covariates", {
  bc <- big_coords()
  y <- bc$city$obs$outcome
  set.seed(0)
  X <- rnorm(bc$coords$n_obs)
  for (kw in list(list(p = 2, deriv = c(1, 0)), list(p = 2, tangvec = c(0, 1)), list(covs = X))) {
    bw <- do.call(rdbwnet, c(list(y, bc$coords), kw))
    expect_true(all(bw$h0[bw$status == "ok"] > 0))
  }
  res <- rdnet(y, bc$coords, covs = X, bwselect = "certwo", conley_cutoff = 200, seed = 0)
  expect_gte(sum(res$table$status == "ok"), 3)
})

test_that("Conley uses a wider auxiliary fit by default", {
  y <- y_mock()
  cc <- mock_coords()
  a <- rdnet(y, cc, h = 600, conley_cutoff = 150, cbands = FALSE)$table
  b <- rdnet(y, cc, h = 600, conley_cutoff = 150, conley_resid_scale = 3, cbands = FALSE)$table
  c3 <- rdnet(y, cc, h = 600, conley_cutoff = 150, conley_resid_order = 3, cbands = FALSE)$table
  expect_equal(a$se_rb, b$se_rb)
  expect_false(isTRUE(all.equal(a$se_rb, c3$se_rb)))
})

test_that("the joint covariance is full and matches the pointwise errors", {
  y <- y_mock()
  cc <- mock_coords()
  for (vce in c("hc1", "cluster", "conley")) {
    args <- list(y, cc, h = 600, vce = vce, cbands = FALSE, seed = 0)
    if (vce == "cluster") args$cluster <- mock_city()$obs$block
    if (vce == "conley") args$conley_cutoff <- 150
    res <- do.call(rdnet, args)
    ok <- res$table$status == "ok"
    expect_equal(dim(res$cov), c(cc$n_eval, cc$n_eval))
    expect_false(all(ok))
    expect_true(all(is.na(res$cov[!ok, ])) && all(is.na(res$cov[, !ok])))
    expect_equal(sqrt(diag(res$cov)[ok]), res$table$se_rb[ok], tolerance = 1e-10)
    expect_false("cb_lo" %in% names(res$table))
    expect_gt(res$crit, 1.96) # stored without bands, for LBATE
  }
})

test_that("WBATE: one-hot weights reproduce a point; equal weights use the covariance", {
  res <- rdnet(y_mock(), mock_coords(), h = 600, seed = 0)
  t <- res$table
  ok <- t$status == "ok"
  j <- which(ok)[3]
  w <- numeric(nrow(t))
  w[j] <- 5 # normalised to one
  one <- wbate(res, w)
  for (k in c("tau", "tau_bc", "se_rb", "ci_lo", "ci_hi", "pvalue")) {
    expect_equal(unname(one[k]), t[[k]][j], tolerance = 1e-10)
  }
  avg <- wbate(res, as.numeric(ok))
  expect_equal(unname(avg["tau_bc"]), mean(t$tau_bc[ok]))
  Sigma <- res$cov[ok, ok]
  expect_equal(unname(avg["se_rb"]), sqrt(sum(Sigma)) / sum(ok))
  # positively correlated neighbours: averaging helps less than independence would suggest
  expect_gt(avg[["se_rb"]], sqrt(sum(t$se_rb[ok]^2)) / sum(ok))
  expect_error(wbate(res, rep(1, nrow(t))), "not estimated")
  expect_error(wbate(res, rep(1, 3)), "one entry")
  expect_error(wbate(res, numeric(nrow(t))), "nonzero")
  expect_error(wbate(res, c(NaN, ok[-1])), "finite")
})

test_that("LBATE", {
  y <- y_mock()
  res <- rdnet(y, mock_coords(), h = 600, seed = 0)
  t <- res$table[res$table$status == "ok", ]
  lb <- lbate(res)
  expect_equal(lb[["tau_bc"]], max(t$tau_bc))
  expect_equal(lb[["tau"]], max(t$tau))
  expect_true(lb[["ci_lo"]] <= lb[["tau_bc"]] && lb[["tau_bc"]] <= lb[["ci_hi"]])
  expect_equal(unname(lb[c("ci_lo", "ci_hi")]), c(max(t$cb_lo), max(t$cb_hi)))
  # same interval without bands in the table (same seed, same critical value)
  nb <- lbate(rdnet(y, mock_coords(), h = 600, seed = 0, cbands = FALSE))
  expect_equal(nb[c("ci_lo", "ci_hi")], lb[c("ci_lo", "ci_hi")])
  s <- summary(res, wbate = as.numeric(res$table$status == "ok"), lbate = TRUE)
  expect_output(print(s), "WBATE")
  expect_output(print(s), "LBATE")
})

test_that("covariates with side bandwidths weight each side by 1 / h^2, as rd2d does", {
  cc <- mock_coords()
  set.seed(1)
  X <- rnorm(cc$n_obs) + 0.002 * cc$s
  y <- 0.0005 * cc$s + 2 * cc$treated + 3 * X + rnorm(cc$n_obs)
  h0 <- 400
  h1 <- 700
  j <- 4
  res <- rdnet(y, cc, h = c(h0, h1), covs = X, cbands = FALSE)
  at <- coords_at(cc, j)
  hs <- ifelse(cc$treated, h1, h0)
  u <- at$d / hs
  w <- pmax(1 - abs(u), 0) / hs^2
  use <- w > 0 & is.finite(at$d) & is.finite(at$o) & is.finite(cc$s)
  df <- data.frame(y = y, T = as.numeric(cc$treated), C = as.numeric(!cc$treated), s = cc$s / hs,
                   o = at$o / hs, X = X, w = w)[use, ]
  fit <- lm(y ~ 0 + T + T:s + T:o + C + C:s + C:o + X, data = df, weights = w)
  expect_equal(res$table$tau[j], unname(coef(fit)["T"] - coef(fit)["C"]), tolerance = 1e-8)
})

test_that("the bandwidth floor counts locations, not repeated measurements", {
  city <- mock_city()
  net <- city$network("length")
  xy <- as.matrix(city$obs[, c("x", "y")])
  c1 <- network_coordinates(net, city$boundary, xy, n_eval = 9)
  c3 <- network_coordinates(net, city$boundary, xy[rep(seq_len(nrow(xy)), each = 3), ], n_eval = 9)
  y <- city$obs$outcome
  b1 <- rdbwnet(y, c1)
  b3 <- rdbwnet(rep(y, each = 3), c3)
  ok <- b1$status == "ok"
  expect_equal(b3$hmin0[ok], b1$hmin0[ok], tolerance = 1e-12)
  expect_equal(b3$hmin1[ok], b1$hmin1[ok], tolerance = 1e-12)
})

test_that("a one-point band uses the normal quantile", {
  city <- mock_city()
  cc <- network_coordinates(city$network("length"), city$boundary, city$obs, eval_points = 520)
  res <- rdnet(city$obs$outcome, cc, h = 600)
  expect_equal(res$crit, qnorm(0.975), tolerance = 1e-12)
})

test_that("min_obs counts locations, not repeated measurements", {
  city <- mock_city()
  net <- city$network("length")
  xy <- as.matrix(city$obs[, c("x", "y")])
  y <- city$obs$outcome
  c1 <- network_coordinates(net, city$boundary, xy, n_eval = 9)
  c3 <- network_coordinates(net, city$boundary, xy[rep(seq_len(nrow(xy)), each = 3), ], n_eval = 9)
  s1 <- rdnet(y, c1, h = 180, cbands = FALSE)$table$status
  s3 <- rdnet(rep(y, each = 3), c3, h = 180, cbands = FALSE)$table$status
  expect_true(any(s1 == "insufficient"))
  expect_equal(s3, s1)
})

test_that("plot accepts a scalar or function truth", {
  res <- rdnet(mock_city()$obs$outcome, mock_coords(), h = 600, cbands = FALSE)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(res, true_tau = 2))
  expect_silent(plot(res, true_tau = function(arc) 1.5 + 0 * arc))
})

test_that("few clusters within the bandwidth warn", {
  city <- mock_city()
  expect_warning(rdnet(city$obs$outcome, mock_coords(), h = 600, cluster = city$obs$x %/% 700, seed = 0),
                 "few clusters")
})

test_that("the selector treats badly conditioned local designs as singular", {
  expect_error(rd2dnet:::.inv_gram(diag(c(1, 1e-10))))
  expect_equal(rd2dnet:::.inv_gram(diag(c(1, 1e-6))), diag(c(1, 1e6)))
})

test_that("the identifiable basis drops aliased and leverage-one columns", {
  set.seed(0)
  s <- runif(40)
  o <- c(rep(0, 20), rep(1, 20))
  X <- cbind(1, s, o, s * s, s * o, o * o)
  expect_equal(rd2dnet:::.identifiable(X, rep(1, 40)), c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
  o3 <- c(rep(0, 20), rep(1, 19), 2)
  X3 <- cbind(1, s, o3, s * s, s * o3, o3 * o3)
  expect_false(rd2dnet:::.identifiable(X3, rep(1, 40))[6])
})

test_that("a singular bias-correction fit does not change the main estimate", {
  city <- mock_city()
  cc <- network_coordinates(city$network("length"), city$boundary, city$obs, n_eval = 15)
  T <- cc$treated
  s <- cc$s
  arc <- cc$crossing_arc
  y <- 1 + 2.5 * T + 0.004 * s + 3e-6 * s^2 - 2e-6 * T * s^2 - 0.0015 * arc + 1e-6 * arc^2
  t <- rdnet(y, cc, h = 500, p = 2, q = 3, cbands = FALSE)$table
  ok <- t$status == "ok"
  expect_equal(t$tau[ok], rep(2.5, sum(ok)), tolerance = 1e-8)
  expect_equal(t$tau_bc[ok], rep(2.5, sum(ok)), tolerance = 1e-6)
})
