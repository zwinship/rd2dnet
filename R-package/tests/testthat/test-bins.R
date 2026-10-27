# Nearest-crossing bins (mirrors tests/test_bins.py).

# Python reference on the mock city: bin_coordinates(coords), rdnet(h = 700, min_obs = 10)
# at the estimated bins, and rdnet_fe(h = 700): tau, se, tau_bc, se_rb.
.py_bins_ref <- list(
  tau_bc = c(3.9983786466343787, 2.110269906086497, 4.586325042139279, -0.15269496419695328, -0.5751017355618021,
             0.001954879349396532),
  se_rb = c(0.7942751421472519, 2.384314290443727, 2.5618744452279025, 1.1378721847563114, 1.6622027716007328,
            1.508649493107737),
  h_S = 3.953410654928521,
  fe = c(1.4329741053916103, 0.2693587968689204, 1.1824789662250355, 0.48376984245544763)
)

# Vertical streets crossing the boundary y = 0 (treated below), joined far away.
# Each street splits 200 m from the boundary on each side: it continues straight
# and a spur leaves at an angle.
comb_fixture <- function(n_streets = 3, spacing = 2000) {
  nodes <- list()
  edges <- list()
  tops <- bots <- character()
  for (k in seq_len(n_streets) - 1) {
    x <- k * spacing
    for (sg in c(1, -1)) {
      tag <- if (sg > 0) "n" else "s"
      j <- paste0(tag, "j", k)
      e <- paste0(tag, "e", k)
      p <- paste0(tag, "p", k)
      nodes[[length(nodes) + 1]] <- data.frame(name = c(j, e, p), x = c(x, x, x + 300),
                                               y = sg * c(200, 600, 500))
      edges[[length(edges) + 1]] <- data.frame(from = c(j, j), to = c(e, p), length = c(400, sqrt(2) * 300),
                                               lanes = c(2, 1))
      if (sg > 0) tops <- c(tops, e) else bots <- c(bots, e)
    }
    edges[[length(edges) + 1]] <- data.frame(from = paste0("nj", k), to = paste0("sj", k), length = 400, lanes = 4)
  }
  for (ends in list(tops, bots)) {
    edges[[length(edges) + 1]] <- data.frame(from = ends[-length(ends)], to = ends[-1], length = spacing, lanes = 1)
  }
  nodes <- do.call(rbind, nodes)
  edges <- do.call(rbind, edges)
  bnd <- boundary(rbind(c(-500, 0), c((n_streets - 1) * spacing + 500, 0)),
                  treated_region = rect_polygon(-1e5, -1e5, 1e5, 0))
  # locations along the streets near the crossings
  set.seed(1)
  xy <- NULL
  for (r in seq_len(nrow(edges))) {
    a <- nodes[nodes$name == edges$from[r], c("x", "y")]
    b <- nodes[nodes$name == edges$to[r], c("x", "y")]
    if (grepl("^.e", edges$from[r]) && grepl("^.e", edges$to[r])) next
    t <- stats::runif(60, 0.02, 0.98)
    xy <- rbind(xy, cbind(a$x + t * (b$x - a$x), a$y + t * (b$y - a$y)))
  }
  xy <- xy[abs(xy[, 2]) > 1, ]
  net <- network_graph(edges, nodes, metric = "length")
  coords <- network_coordinates(net, bnd, xy, n_eval = 3)
  list(nodes = nodes, edges = edges, coords = coords, xy = xy)
}

comb <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- comb_fixture()
    cache
  }
})

test_that("bins are the nearest crossings", {
  cb <- comb()
  b <- bin_coordinates(cb$coords, features = "lanes", topology = FALSE)
  expect_equal(b$n_eval, 3L)
  expect_equal(b$eval_node[b$bin], cb$coords$crossing_node)
  at <- coords_at(b, 2)
  m <- b$bin == 2
  expect_equal(at$d[m], abs(cb$coords$s[m]))
  expect_true(all(is.infinite(at$d[!m])))
  expect_true(all(at$o[m] == 0) && all(is.na(at$o[!m])))
  tab <- bins_table(b)
  expect_equal(sum(tab$n0 + tab$n1), cb$coords$n_obs)
})

test_that("branches split at the junction", {
  cb <- comb()
  b <- bin_coordinates(cb$coords, features = "lanes", topology = FALSE)
  y <- cb$xy[, 2]
  xoff <- cb$xy[, 1] - round(cb$xy[, 1] / 2000) * 2000
  expect_true(all(b$branch[abs(y) < 200 - 1e-6] == 0))
  straight <- abs(y) > 200 + 1e-6 & abs(xoff) < 1e-6
  spur <- abs(y) > 200 + 1e-6 & xoff > 1e-6
  expect_true(all(b$branch[straight] > 0) && all(b$branch[spur] > 0))
  for (m in 1:3) {
    for (side in c(TRUE, FALSE)) {
      sel <- b$bin == m & b$treated == side
      expect_length(unique(b$branch[sel & straight]), 1)
      expect_length(unique(b$branch[sel & spur]), 1)
      expect_false(unique(b$branch[sel & straight]) == unique(b$branch[sel & spur]))
    }
  }
  expect_equal(lengths(b$bin_branches), c(4L, 4L, 4L))
})

test_that("S is the Mahalanobis distance of the lanes", {
  cb <- comb()
  b <- bin_coordinates(cb$coords, features = "lanes", topology = FALSE)
  sd_l <- stats::sd(cb$edges$lanes)
  y <- cb$xy[, 2]
  xoff <- cb$xy[, 1] - round(cb$xy[, 1] / 2000) * 2000
  expect_equal(b$similarity, ifelse(abs(y) < 200, 0, ifelse(abs(xoff) < 1e-6, 2, 3)) / sd_l)
  expect_equal(b$feature_rank, 1L)
  expect_equal(b$h_S, sqrt(2 * stats::qchisq(0.95, 1)))
  expect_equal(b$kernel_weights, pmax(1 - b$similarity / b$h_S, 0))
  expect_null(with_similarity_bandwidth(b, Inf)$kernel_weights)
  expect_error(with_similarity_bandwidth(b, -1))
})

test_that("branch slopes recover the jump exactly", {
  cb <- comb()
  b <- bin_coordinates(cb$coords, features = "lanes", topology = FALSE)
  Tn <- as.numeric(cb$coords$treated)
  s <- cb$coords$s
  brs <- sort(unique(b$branch[b$branch > 0]))
  sl <- c(0.002, 0.001 * (seq_along(brs) + 2) * ifelse(seq_along(brs) %% 2 == 0, 1, -1))
  slope <- sl[match(b$branch, c(0, brs))]
  y <- c(5, -1, 0.5)[b$bin] + c(1, 2, 3)[b$bin] * Tn + slope * s
  res <- rdnet(y, with_similarity_bandwidth(b, Inf), h = 1500, min_obs = 5)
  expect_true(all(res$table$status == "ok"))
  expect_equal(res$table$tau, c(1, 2, 3), tolerance = 1e-8)
  expect_equal(res$table$tau_bc, c(1, 2, 3), tolerance = 1e-8)
  plain <- bin_coordinates(cb$coords, similarity = FALSE, branches = FALSE)
  res2 <- rdnet(y, plain, h = 1500, min_obs = 5)
  expect_gt(max(abs(res2$table$tau - c(1, 2, 3))), 0.05)
})

test_that("topology features on a triangle with a tail", {
  nodes <- data.frame(name = c("a", "b", "c", "d"), x = c(0, 1, 0.5, 0.5), y = c(0, 0, 1, 2))
  edges <- data.frame(from = c("a", "b", "c", "c"), to = c("b", "c", "a", "d"), length = 1)
  net <- network_graph(edges, nodes, metric = "length")
  E <- edge_frame(net)
  ef <- .edge_feature_matrix(net, E, character(), TRUE)
  expect_equal(ef$names, c("n_adjacent", "n_second", "degree_max", "clustering"))
  lab <- net$nodes$name
  key <- vapply(seq_len(nrow(E)), function(r) paste(sort(lab[c(E$u[r], E$v[r])]), collapse = ""), "")
  Z <- ef$Z[ef$orig, , drop = FALSE]
  expect_equal(unname(Z[key == "ab", ]), c(2, 1, 2, 1))
  expect_equal(unname(Z[key == "bc", ]), c(3, 0, 3, 2 / 3))
  expect_equal(unname(Z[key == "ac", ]), c(3, 0, 3, 2 / 3))
  expect_equal(unname(Z[key == "cd", ]), c(2, 1, 3, 1 / 6))
})

test_that("split edges keep the topology of the network as given", {
  cb <- comb()
  net <- cb$coords$network
  ef <- .edge_feature_matrix(net, edge_frame(net), character(), TRUE)
  expect_equal(nrow(ef$Z), nrow(cb$edges))
  ref <- network_graph(cb$edges, cb$nodes, metric = "length")
  ef0 <- .edge_feature_matrix(ref, edge_frame(ref), character(), TRUE)
  expect_equal(apply(ef$Z, 2, sort), apply(ef0$Z, 2, sort))
})

test_that("wbate density and gls rules", {
  city <- mock_city()
  coords <- mock_coords()
  res <- rdnet(city$obs$outcome, coords, h = 600)
  ok <- res$table$status == "ok"
  n <- res$table$n0 + res$table$n1
  expect_equal(wbate(res, "density")[["tau_bc"]], wbate(res, ifelse(ok, n, 0))[["tau_bc"]])
  S <- res$cov[ok, ok]
  g <- solve(S, rep(1, sum(ok)))
  g <- g / sum(g)
  gls <- expect_warning(wbate(res, "gls"), "negative")
  expect_equal(gls[["tau_bc"]], sum(g * res$table$tau_bc[ok]))
  expect_equal(gls[["se_rb"]], sqrt(drop(g %*% S %*% g)))
  expect_lte(gls[["se_rb"]], wbate(res, "density")[["se_rb"]] + 1e-12)
  expect_equal(wbate(res, "place"), wbate(res))
  expect_error(wbate(res, "median"))
})

test_that("rdnet_fe recovers the mass-weighted average effect", {
  cb <- comb()
  Tn <- as.numeric(cb$coords$treated)
  plain <- bin_coordinates(cb$coords, similarity = FALSE, branches = FALSE)
  tau_m <- c(1, 2, 4)
  y <- c(5, -1, 0.5)[plain$bin] + tau_m[plain$bin] * Tn
  fe <- rdnet_fe(y, plain, h = 500, p = 0, q = 0, bin_slopes = FALSE)
  W <- pmax(1 - abs(cb$coords$s) / 500, 0)
  mass <- as.numeric(tapply(W, plain$bin, sum))
  expect_equal(fe$tau, sum(mass * tau_m) / sum(mass), tolerance = 1e-10)
  w1 <- as.numeric(tapply(W * Tn, plain$bin, sum))
  ph <- w1 / mass
  fe0 <- rdnet_fe(y, plain, h = 500, p = 0, q = 0, bin_slopes = FALSE, ipw = FALSE)
  om <- mass * ph * (1 - ph)
  expect_equal(fe0$tau, sum(om * tau_m) / sum(om), tolerance = 1e-10)
  b <- bin_coordinates(cb$coords, features = "lanes")
  y2 <- c(5, -1, 0.5)[b$bin] + 2.5 * Tn + 0.003 * cb$coords$s - 0.001 * Tn * cb$coords$s
  fe2 <- rdnet_fe(y2, b, h = 700)
  expect_equal(fe2$tau, 2.5, tolerance = 1e-8)
  expect_equal(fe2$tau_bc, 2.5, tolerance = 1e-8)
  expect_true("R*S" %in% fe2$terms)
})

test_that("bins run with every variance option", {
  city <- mock_city()
  coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
  b <- bin_coordinates(coords)
  y <- city$obs$outcome
  opts <- list(list(), list(cluster = city$obs$block), list(conley_cutoff = 100),
               list(fuzzy = city$obs$treated), list(covs = city$obs$far))
  for (o in opts) {
    res <- suppressWarnings(do.call(rdnet, c(list(y, b, h = 700, min_obs = 10), o)))
    expect_true(any(res$table$status == "ok"))
    expect_true(is.finite(wbate(res, "density")[["se_rb"]]))
  }
  for (o in list(list(), list(cluster = city$obs$block), list(conley_cutoff = 100), list(vce = "hc1"))) {
    fe <- suppressWarnings(do.call(rdnet_fe, c(list(y, b, h = 700), o)))
    expect_true(is.finite(fe$se_rb) && fe$se_rb > 0)
  }
  sens <- similarity_sensitivity(y, b, h = 700, min_obs = 10)
  expect_equal(nrow(sens), 4L)
  expect_true(is.infinite(sens$h_S[4]))
  expect_output(print(b), "bins=")
})

test_that("bins match the Python reference on the mock city", {
  city <- mock_city()
  coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
  b <- bin_coordinates(coords)
  res <- rdnet(city$obs$outcome, b, h = 700, min_obs = 10)
  ok <- res$table$status == "ok"
  expect_equal(res$table$tau_bc[ok], .py_bins_ref$tau_bc, tolerance = 1e-8)
  expect_equal(res$table$se_rb[ok], .py_bins_ref$se_rb, tolerance = 1e-8)
  expect_equal(b$h_S, .py_bins_ref$h_S, tolerance = 1e-12)
  fe <- rdnet_fe(city$obs$outcome, b, h = 700)
  expect_equal(c(fe$tau, fe$se, fe$tau_bc, fe$se_rb), .py_bins_ref$fe, tolerance = 1e-8)
})
