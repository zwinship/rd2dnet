# Live OpenStreetMap download, the counterpart of the Python package's
# `pytest -m online` test. Skipped unless RD2DNET_ONLINE=true (it needs
# internet access, dodgr, osmdata and sf); see .github/workflows/online.yml.

test_that("a real OpenStreetMap network downloads and works end to end", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RD2DNET_ONLINE"), "true"), "set RD2DNET_ONLINE=true to run online tests")
  skip_if_not_installed("dodgr")
  skip_if_not_installed("sf")
  net <- network_graph_from_place("Midtown Manhattan, New York", offline_fallback = FALSE)
  expect_gt(net$n_nodes, 100)
  expect_match(net$crs, "^EPSG:326")
  expect_true(all(is.finite(net$edges$length) & net$edges$length >= 0))
  # a boundary across the middle of the network, with points on both sides
  cx <- stats::median(net$coords[, 1])
  cy <- stats::median(net$coords[, 2])
  span <- diff(range(net$coords[, 1]))
  b <- boundary(rbind(c(cx - span, cy), c(cx + span, cy)), treated_side = "left", tol = 25)
  set.seed(1)
  pts <- cbind(cx + stats::runif(400, -300, 300), cy + stats::runif(400, -300, 300))
  cc <- suppressWarnings(network_coordinates(net, b, pts, n_eval = 5))
  expect_gt(mean(is.finite(cc$s)), 0.9)
  expect_true(any(cc$treated) && any(!cc$treated))
})
