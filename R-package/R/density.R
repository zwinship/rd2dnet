# Network-adapted density (sorting) test at the boundary.
#
# The estimator and jackknife variance are ported from rddensity (Cattaneo,
# Jansson and Ma; GPL-3) and adapted to networks: the empirical CDF of the
# signed network distance s is regressed on a polynomial in s separately on
# each side within a bandwidth, and the density at the boundary is the slope
# at s = 0. With normalize = "network", s is re-expressed as cumulative street
# length, so the density is observations per metre of street.

# Local polynomial CDF fit: raw coefficients, their jackknife covariance and
# window counts. Columns are ordered [1_L, 1_R, u_L, u_R, u_L^2, u_R^2, ...];
# the density on each side is the coefficient on u (positions 3 and 4).
.cjm <- function(x, om, hl, hr, p, kfun) {
  o <- order(x, method = "radix")
  x <- x[o]
  om <- om[o]
  denom <- sum(om) - mean(om)
  F <- (cumsum(om) - om) / denom
  win <- x >= -hl & x <= hr
  xw <- x[win]
  Fw <- F[win]
  omw <- om[win]
  left <- xw < 0
  u <- ifelse(left, xw / hl, xw / hr)
  K <- ifelse(left, kfun(abs(u)) / hl, kfun(abs(u)) / hr)
  keep <- K > 0
  xw <- xw[keep]
  Fw <- Fw[keep]
  omw <- omw[keep]
  left <- left[keep]
  u <- u[keep]
  K <- K[keep]
  nl <- sum(left)
  nr <- sum(!left)
  if (min(nl, nr) <= p + 1) stop("too few observations within the bandwidth on one side", call. = FALSE)
  Xp <- matrix(0, length(xw), 2 * p + 2)
  Hp <- numeric(2 * p + 2)
  for (k in 0:p) {
    Xp[, 2 * k + 1] <- ifelse(left, u^k, 0)
    Xp[, 2 * k + 2] <- ifelse(!left, u^k, 0)
    Hp[2 * k + 1] <- hl^k
    Hp[2 * k + 2] <- hr^k
  }
  XpW <- Xp * K
  Sinv <- solve(crossprod(XpW, Xp))
  b <- drop(Sinv %*% crossprod(XpW, Fw)) / Hp
  n <- nrow(XpW)
  tail <- apply(XpW[n:1, , drop = FALSE], 2, cumsum)[n:1, , drop = FALSE] # sum over k with x_k >= x_i
  L <- tail * (omw / denom)
  Hinv <- diag(1 / Hp, length(Hp))
  V <- Hinv %*% Sinv %*% crossprod(L) %*% Sinv %*% Hinv
  list(b = b, V = V, nl = nl, nr = nr)
}

# MSE-optimal bandwidth for one side's density (plug-in, rd2d-style constants).
.side_bandwidth <- function(x, om, right, p, kfun, scaleregul, n_min) {
  side <- if (right) x >= 0 else x < 0
  a <- abs(x[side])
  if (length(a) < n_min) stop("too few observations on one side for bandwidth selection", call. = FALSE)
  sdx <- sqrt(mean((x - mean(x))^2))
  lo <- .kth(a, n_min)
  hi <- max(a)
  h0 <- min(max(2 * sdx * length(x)^(-1 / (2 * p + 1)), lo), hi)
  bn <- min(max(stats::median(a), lo), hi)
  col <- if (right) 4 else 3 # density coefficient
  # Variance constant at the pilot.
  Cv <- .cjm(x, om, h0, h0, p, kfun)$V[col, col] * h0
  # Bias constant: degree p+1 coefficient from an order p+1 fit at bn, projected on the slope.
  f1 <- .cjm(x, om, bn, bn, p + 1, kfun)
  gi <- 2 * (p + 1) + (col - 2)
  gamma <- f1$b[gi]
  var_gamma <- f1$V[gi, gi]
  # Projection of u^(p+1) onto the slope at h0, in scaled units.
  xs <- x[x >= -h0 & x <= h0]
  u <- xs / h0
  K <- kfun(abs(u)) / h0
  on <- if (right) xs >= 0 else xs < 0
  us <- u[on]
  Ks <- K[on]
  X <- outer(us, 0:p, `^`)
  proj <- solve(crossprod(X * Ks, X), crossprod(X * Ks, us^(p + 1)))[2]
  Cb2 <- (proj * gamma)^2 + scaleregul * proj^2 * var_gamma
  h <- (Cv / (2 * p * max(Cb2, 1e-300)))^(1 / (2 * p + 1))
  min(max(h, lo), hi)
}

# Cumulative street length within distance a of the boundary, per side:
# list(left, right) of functions.
# Distance along the boundary between arc positions (the shorter way round a closed ring).
.arc_gap <- function(b, arc, center) {
  gap <- abs(arc - center)
  if (isTRUE(b$closed)) gap <- pmin(gap, b$length - gap)
  gap
}

.street_length_transform <- function(coords, center, window) {
  net <- coords$network
  E <- edge_frame(net)
  u <- E$u
  v <- E$v
  on <- coords$node_on_boundary
  nd <- ifelse(on, 0, coords$node_dist)
  side <- ifelse(on[u], coords$node_treated[v], coords$node_treated[u])
  cost <- E$cost
  length <- E$length
  keep <- !(on[u] & on[v]) & is.finite(nd[u]) & is.finite(nd[v])
  src <- coords$node_source
  su <- ifelse(is.na(src[u]), 1L, src[u])
  sv <- ifelse(is.na(src[v]), 1L, src[v])
  arc <- (boundary_project(coords$boundary, net$coords[su, , drop = FALSE]) +
            boundary_project(coords$boundary, net$coords[sv, , drop = FALSE])) / 2
  keep <- keep & !boundary_is_masked(coords$boundary, arc)
  if (!is.null(center)) keep <- keep & .arc_gap(coords$boundary, arc, center) <= window
  exact <- is.finite(cost)
  med <- if (any(is.finite(cost[keep]))) {
    stats::median(cost[keep], na.rm = TRUE)
  } else {
    stats::median(abs(nd[u] - nd[v])[keep], na.rm = TRUE)
  }
  lapply(c(FALSE, TRUE), function(flag) {
    ks <- which(keep & side == flag)
    if (!length(ks)) stop("no streets on one side of the boundary in this window", call. = FALSE)
    a <- nd[u[ks]]
    b <- nd[v[ks]]
    ex <- exact[ks]
    cc <- ifelse(ex, cost[ks], abs(b - a))
    # Sample each street at m points (more on long ones) and spread its length over them.
    m <- pmax(4, ceiling(16 * cc / max(med, 1e-12)))
    e <- rep(seq_along(ks), m)
    tt <- (seq_len(sum(m)) - 1 - rep(cumsum(m) - m, m) + 0.5) / m[e]
    xs <- ifelse(ex[e], pmin(a[e] + tt * cc[e], b[e] + (1 - tt) * cc[e]), a[e] + tt * (b[e] - a[e]))
    ws <- (length[ks] / m)[e]
    o <- order(xs, method = "radix")
    grid <- c(0, xs[o])
    cum <- c(0, cumsum(ws[o]))
    function(z) .interp(z, grid, cum)
  })
}

#' Density (sorting) test at a network boundary
#'
#' Tests for a discontinuity in the density of observations at the boundary
#' with the local polynomial density estimator of Cattaneo, Jansson and Ma
#' (as in \pkg{rddensity}), with jackknife standard errors and bias-corrected
#' inference, applied to the signed network distance `s`.
#'
#' With `normalize = "network"`, each observation's distance is
#' re-expressed as cumulative street length: `G(a)` is the total length of
#' street on that side lying within network distance `a` of the boundary, and
#' the test runs on `sign * G(|s|)`. The density in that coordinate is
#' observations per metre of street, so jumps in street supply at the
#' boundary (a highway on one side, a dense grid on the other, a park
#' frontage) do not masquerade as sorting. Bandwidths are then in
#' street-length units. Masked stretches of the boundary are always excluded.
#'
#' @param x Output of [network_coordinates()].
#' @param eval_index,window Restrict to observations (and streets) whose
#'   boundary crossing lies within `window` (arc length) of evaluation point
#'   `eval_index` (1-based). By default the whole boundary is pooled.
#' @param h Bandwidth, scalar or `c(h_control, h_treated)`, in metric units
#'   (street-length units with `normalize = "network"`). Selected by an MSE
#'   plug-in rule per side when omitted.
#' @param p,q Local polynomial order of the CDF fit (the density is its slope)
#'   and of the bias-corrected fit used for inference (default `p + 1`).
#' @param normalize The running variable. `"euclidean"` (default): the signed
#'   straight-line distance from each location to the boundary line; sorting is
#'   about where units locate, and this is the test that held its size (0.03-0.06
#'   at 5%) on straight, round, wavy and irregular boundaries in our Monte Carlo,
#'   with the most power against one-sided sorting. `"none"`: the signed network
#'   distance `s`; its density jumps at curved or irregular boundaries without any
#'   sorting. `"network"`: observations per metre of street; valid only when
#'   observations lie on the streets of a straight boundary. No density test holds
#'   its size everywhere when the observations lie on or next to the network.
#' @param kernel Kernel name or function (see [get_kernel()]).
#' @param weights Optional non-negative observation weights.
#' @param bwselect `"common"` (one bandwidth, the smaller of the two side-wise
#'   MSE-optimal choices, for both sides) or `"each"` (separate bandwidths).
#'   On street grids the density of `s` is jagged, and unequal bandwidths let
#'   unequal smoothing bias pass for a discontinuity.
#' @param scaleregul Weight on the variance of the bias estimate in the
#'   bandwidth rule.
#' @param unique_locations Count each location once (the first observation
#'   there). Use this with repeated measurements (panels, repeat sales, sensor
#'   readings): the sorting question is about where units are, and repeated
#'   rows at the same place otherwise look like a density spike.
#' @param ... Ignored.
#'
#' @return An object of class `rdnet_density`: the density jump `diff`
#'   (`f+ - f-`, order `p`) and `diff_bc` (bias-corrected), its jackknife
#'   standard error `se`, `z`, `pvalue`, the descriptive log ratio `theta`
#'   (`log(f+/f-)`), side-wise densities `f_left`/`f_right`, bandwidths
#'   `h_left`/`h_right`, sample sizes `n_left`/`n_right` and within the
#'   bandwidth `n_left_h`/`n_right_h`, `p`, `q` and `normalize`.
#' @references
#' Cattaneo, M. D., Jansson, M. and Ma, X. (2020). Simple Local Polynomial
#' Density Estimators. Journal of the American Statistical Association
#' 115(531): 1449-1455.
#'
#' Cattaneo, M. D., Jansson, M. and Ma, X. (2018). Manipulation Testing Based
#' on Density Discontinuity. Stata Journal 18(1): 234-261.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' density(coords)
#' density(coords, normalize = "none")
#' @importFrom stats density
#' @export
density.rd2dnet_coords <- function(x, eval_index = NULL, window = NULL, h = NULL, p = 2, q = NULL,
                                   normalize = "euclidean", kernel = "triangular", weights = NULL,
                                   bwselect = "common", scaleregul = 1, unique_locations = FALSE, ...) {
  coords <- x
  if (!normalize %in% c("euclidean", "none", "network")) {
    stop("normalize must be 'euclidean', 'none' or 'network'", call. = FALSE)
  }
  .warn_if_on_network(coords)
  q <- if (is.null(q)) p + 1 else as.integer(q)
  if (p < 1 || q < p) stop("need 1 <= p <= q", call. = FALSE)
  if (!bwselect %in% c("common", "each")) stop("bwselect must be 'common' or 'each'", call. = FALSE)
  kfun <- get_kernel(kernel)
  T <- coords$treated
  om <- if (is.null(weights)) rep(1, coords$n_obs) else .check_weights(weights, coords$n_obs)
  b <- coords$boundary
  if (normalize == "euclidean") {
    # straight-line distance to the boundary line, where it is reached in the plane
    s <- boundary_distance(b, coords$xy)
    arc <- boundary_project(b, coords$xy)
    reach <- !boundary_is_masked(b, arc)
    if (!isTRUE(b$closed)) reach <- reach & arc > 1e-9 * b$length & arc < b$length * (1 - 1e-9)
  } else {
    s <- abs(coords$s)
    arc <- coords$crossing_arc
    a0 <- arc
    a0[is.na(a0)] <- -1
    reach <- !boundary_is_masked(b, a0)
  }
  # Masked stretches (parks, rivers) are void by construction: leave them out.
  sel <- is.finite(s) & om > 0 & reach
  first <- !duplicated(coords$loc)
  if (isTRUE(unique_locations)) {
    sel <- sel & first
  } else if (any(sel & !first)) {
    warning("several observations share a location (repeated measurements?); the density test counts ",
            "each one. Use unique_locations = TRUE to test the density of locations", call. = FALSE)
  }
  center <- NULL
  if (!is.null(eval_index)) {
    if (is.null(window)) stop("give a window (arc length) together with eval_index", call. = FALSE)
    center <- coords$eval_node_arc[eval_index]
    sel <- sel & !is.na(arc) & .arc_gap(coords$boundary, arc, center) <= window
  }
  if (normalize == "network") {
    G <- .street_length_transform(coords, center, window)
    sv <- s
    sv[!is.finite(sv)] <- 0
    s <- ifelse(T, G[[2]](sv), G[[1]](sv))
  }
  # controls exactly on the boundary (s = 0) stay on the control side
  xs <- ifelse(T, s, ifelse(s == 0, -.Machine$double.xmin, -s))[sel]
  om <- om[sel]
  if (length(xs) < 20 || sum(xs >= 0) < 10 || sum(xs < 0) < 10) {
    stop("too few observations for a density test", call. = FALSE)
  }
  if (is.null(h)) {
    n_min <- max(10, 2 * (q + 2))
    hl <- .side_bandwidth(xs, om, FALSE, p, kfun, scaleregul, n_min)
    hr <- .side_bandwidth(xs, om, TRUE, p, kfun, scaleregul, n_min)
    if (bwselect == "common") hl <- hr <- min(hl, hr)
  } else {
    hh <- rep_len(as.numeric(h), 2)
    hl <- hh[1]
    hr <- hh[2]
  }
  fp <- .cjm(xs, om, hl, hr, p, kfun)
  fq <- .cjm(xs, om, hl, hr, q, kfun)
  fl <- fp$b[3]
  fr <- fp$b[4]
  diff <- fr - fl
  diff_bc <- fq$b[4] - fq$b[3]
  Vlr <- fq$V[3:4, 3:4]
  se <- sqrt(max(Vlr[1, 1] + Vlr[2, 2] - 2 * Vlr[1, 2], 0))
  zstat <- if (se > 0) diff_bc / se else NA_real_
  structure(list(
    diff = diff, diff_bc = diff_bc, se = se, z = zstat, pvalue = 2 * stats::pnorm(-abs(zstat)),
    theta = if (fl > 0 && fr > 0) log(fr / fl) else NA_real_, f_left = fl, f_right = fr,
    h_left = hl, h_right = hr, n_left = sum(xs < 0), n_right = sum(xs >= 0),
    n_left_h = fp$nl, n_right_h = fp$nr, p = p, q = q, normalize = normalize
  ), class = "rdnet_density")
}

#' @export
print.rdnet_density <- function(x, ...) {
  cat(sprintf(paste0("rd2dnet density test (%s): jump=%.3g (bc %.3g) se=%.3g z=%.3f p=%.4f ",
                     "log-ratio=%.3f\n  [n=%d/%d, in window %d/%d, h=%.1f/%.1f]\n"),
              x$normalize, x$diff, x$diff_bc, x$se, x$z, x$pvalue, x$theta, x$n_left, x$n_right,
              x$n_left_h, x$n_right_h, x$h_left, x$h_right))
  invisible(x)
}
