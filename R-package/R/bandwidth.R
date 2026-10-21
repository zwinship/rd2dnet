# Data-driven bandwidth selection for network boundary discontinuity estimation.
#
# Ported from rd2d (Cattaneo, Titiunik and Yu; GPL-3; rdbw2d's direct plug-in,
# method = "dpi") and adapted to network coordinates. Kernel weights are
# K(d_i(b) / h) / h^2 with d the network distance (rd2d's radial kernel) and
# the local basis is in (s / h, o / h). At each evaluation point and side:
#   1. a rule-of-thumb pilot dn from a bivariate normal reference for (s, o(b)),
#      clipped to [bwcheck-th nearest, farthest] network distance;
#   2. the variance constant V from the order-p fit at dn;
#   3. a DPI pilot bn for the leading bias, chosen by the same MSE rule with its
#      bias taken from an order p + 2 fit at the median distance;
#   4. the bias constant B, its variance Reg1 (weighted by scaleregul) and the
#      next-order bias Reg2 (weighted by scalebiascrct) from an order p + 2 fit
#      at the median distance.
# Then h = (2 V / ((2p + 2) ((B + Reg2)^2 + scaleregul Reg1)))^(1 / (2p + 4)).

.BWSELECT <- c("mserd", "msetwo", "imserd", "imsetwo", "cerrd", "certwo", "icerrd", "icertwo")

# (mu_2(K)^2, R(K)) style constants used by rd2d's rule of thumb.
.ROT_CONST <- list(triangular = c(3 / 20, 3 / (2 * pi)), epanechnikov = c(1 / 6, 4 / (3 * pi)),
                   uniform = c(1 / 4, 1 / pi), gaussian = c(1, 1 / (4 * pi)))

.default_min_obs <- function(q) 20L + as.integer(.n_terms(q))

# Inverse of a local Gram matrix; errors when the local design is degenerate.
# The bases are scaled by the bandwidth, so well-posed designs have condition
# numbers far below 1e8; above it the plug-in constants are unreliable (they can
# change sign with rounding) and the design counts as singular.
.inv_gram <- function(A) {
  if (!all(is.finite(A)) || .cond(A) > 1e8) stop(.singular())
  chol2inv(chol(A))
}

# rd2d's bivariate normal-reference rule of thumb (a univariate one with dim = 1).
.rot <- function(s, o, kname, n_eff, dim = 2) {
  if (dim == 1) return(1.06 * stats::sd(s) * n_eff^(-1 / 5))
  k <- .ROT_CONST[[kname]] %||% .ROT_CONST$triangular
  mu2 <- k[1]
  l2 <- k[2]
  S <- stats::cov(cbind(s, o))
  inv <- .pinv(S, rcond = 1e-20)
  det <- det(S)
  sqrt_det <- if (is.finite(det) && det > 0) {
    sqrt(det)
  } else {
    sqrt(max(prod(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 1e-300))
  }
  trace_const <- 1 / (2^4 * pi * sqrt_det) * (2 * sum(diag(inv %*% inv)) + sum(diag(inv))^2)
  ((2 * l2) / (n_eff * mu2 * trace_const))^(1 / 6)
}

# Local data on one side of the boundary at one evaluation point. Kernel
# weights and the basis depend only on h, and the order-`order` basis is a
# column prefix of any higher-order one, so they are cached per h.
.side <- function(s, o, d, y, om, kfun, cluster = NULL, dim = 2) {
  self <- new.env(parent = emptyenv())
  self$dim <- dim
  self$s <- s
  self$o <- o
  self$d <- d
  self$y <- y
  self$om <- om
  self$kfun <- kfun
  self$cluster <- cluster
  self$basis <- list()
  self
}

.side_design <- function(sd, h, order, y = NULL) {
  key <- sprintf("%.17g", h)
  hit <- sd$basis[[key]]
  if (is.null(hit) || hit$order < order) {
    w <- sd$kfun(sd$d / h) / h^2 * sd$om
    keep <- which(w > 0)
    hit <- list(keep = keep, w = w[keep], R = .poly_basis(sd$s[keep] / h, sd$o[keep] / h, order, sd$dim),
                order = order)
    sd$basis[[key]] <- hit
  }
  yy <- (if (is.null(y)) sd$y else y)[hit$keep]
  cl <- if (is.null(sd$cluster)) NULL else sd$cluster[hit$keep]
  list(R = hit$R[, seq_len(.n_terms(order, sd$dim)), drop = FALSE], w = hit$w, y = yy, cl = cl)
}

.bw_sigma <- function(WR, resd, h, cl) {
  scores <- WR * resd
  if (is.null(cl)) {
    return(crossprod(scores) * h^2)
  }
  summed <- rowsum(scores, cl)
  n <- length(cl)
  k <- ncol(scores)
  g <- nrow(summed)
  adj <- if (g > 1 && n > k) ((n - 1) / (n - k)) * (g / (g - 1)) else 1
  crossprod(summed) * h^2 * adj
}

.bw_adjust <- function(resd, R, w, invG, vce) {
  n <- nrow(R)
  k <- ncol(R)
  if (vce == "hc1" && n > k) {
    return(resd * sqrt(n / (n - k)))
  }
  if (vce %in% c("hc2", "hc3")) {
    hii <- w * rowSums((R %*% invG) * R)
    return(if (vce == "hc2") resd / sqrt(pmax(1 - hii, 1e-8)) else resd / pmax(1 - hii, 1e-8))
  }
  resd
}

# Raw-scale coefficients (and their scaled covariance constant).
.side_lm <- function(sd, h, order, vce = "hc1", varr = FALSE, y = NULL) {
  ds <- .side_design(sd, h, order, y)
  R <- ds$R
  w <- ds$w
  RtW <- t(R * w)
  invG <- .inv_gram(RtW %*% R)
  beta_s <- drop(invG %*% (RtW %*% ds$y))
  e <- .poly_exponents(order, sd$dim)
  beta <- beta_s / h^(e[, "a"] + e[, "c"])
  cov <- NULL
  if (varr) {
    resd <- .bw_adjust(ds$y - drop(R %*% beta_s), R, w, invG, vce)
    cov <- invG %*% .bw_sigma(R * w, resd, h, ds$cl) %*% invG
  }
  list(beta = beta, cov = cov)
}

# Projection of the degree p + 1 terms onto the target `vec`: rd2d's vec_q.
.side_coeff <- function(sd, dn, p, vec) {
  ds <- .side_design(sd, dn, p + 1)
  k <- .n_terms(p, sd$dim)
  Rp <- ds$R[, seq_len(k), drop = FALSE]
  S <- ds$R[, -seq_len(k), drop = FALSE]
  RtW <- t(Rp * ds$w)
  invG <- .inv_gram(RtW %*% Rp)
  c(rep(0, k), drop(vec %*% invG %*% (RtW %*% S)))
}

.side_consts <- function(sd, p, vec, dn, bn1, bn2, vce, y = NULL) {
  ds <- .side_design(sd, dn, p + (if (is.null(bn2)) 1 else 2), y)
  kp <- .n_terms(p, sd$dim)
  kp1 <- .n_terms(p + 1, sd$dim)
  R <- ds$R
  w <- ds$w
  Rp <- R[, seq_len(kp), drop = FALSE]
  S <- R[, (kp + 1):kp1, drop = FALSE]
  WR <- Rp * w
  invG <- .inv_gram(crossprod(WR, Rp))
  vec_q <- c(rep(0, kp), drop(vec %*% invG %*% crossprod(WR, S)))
  beta_s <- drop(invG %*% crossprod(WR, ds$y))
  resd <- .bw_adjust(ds$y - drop(Rp %*% beta_s), Rp, w, invG, vce)
  V <- drop(vec %*% invG %*% .bw_sigma(WR, resd, dn, ds$cl) %*% invG %*% vec)
  fit1 <- .side_lm(sd, bn1, p + 1, vce, varr = TRUE, y = y)
  B <- sum(vec_q * fit1$beta)
  reg1 <- drop(vec_q %*% fit1$cov %*% vec_q) / bn1^(2 + 2 * (p + 1))
  # a quadratic form in a covariance matrix: negative values are rounding
  # (a degenerate bias direction), and their sign must not decide the fit
  if (is.finite(reg1) && reg1 < 0) reg1 <- 0
  reg2 <- NA_real_
  if (!is.null(bn2)) {
    Tm <- R[, -seq_len(kp1), drop = FALSE]
    vec_t <- c(rep(0, kp1), drop(vec %*% invG %*% crossprod(WR, Tm)))
    beta2 <- .side_lm(sd, bn2, p + 2, vce, y = y)$beta
    reg2 <- dn * sum(vec_t * beta2)
  }
  c(V = V, B = B, Reg1 = reg1, Reg2 = reg2)
}

# rd2d's shrinkage from the MSE-optimal to the coverage-error-optimal bandwidth.
.cer_factor <- function(n, p) {
  n <- max(as.integer(n), 1L)
  n^(1 / (2 * p + 4) - 1 / (p + 4))
}

# Covariate coefficients common to both sides, partialling out side-specific polynomials.
.covariate_gamma <- function(sides, Zs, dns, p, outcomes) {
  A <- 0
  b <- 0
  for (t in seq_along(sides)) {
    sd <- sides[[t]]
    dn <- dns[t]
    w <- sd$kfun(sd$d / dn) / dn^2 * sd$om
    k <- w > 0
    R <- .poly_basis(sd$s[k] / dn, sd$o[k] / dn, p, sd$dim)
    W <- w[k]
    Zk <- Zs[[t]][k, , drop = FALSE]
    yk <- outcomes[[t]][k]
    RtW <- t(R * W)
    invG <- .pinv(RtW %*% R, rcond = 1e-12)
    uz <- RtW %*% Zk
    A <- A + crossprod(Zk * W, Zk) - t(uz) %*% invG %*% uz
    b <- b + crossprod(Zk * W, yk) - t(uz) %*% invG %*% (RtW %*% yk)
  }
  drop(.lstsq(A, b))
}

#' Bandwidth selection for network boundary discontinuity designs
#'
#' Selects MSE- or coverage-error-optimal bandwidths (in metric units) at each
#' evaluation point. A port of rd2d's direct plug-in selector (`rdbw2d`) to
#' network coordinates: a rule-of-thumb pilot, a separate pilot bandwidth for
#' the bias, regularisation (`scaleregul`) and a higher-order bias term
#' (`scalebiascrct`).
#'
#' `deriv`/`tangvec` target a derivative jump and `covs` are partialled out
#' (with coefficients common to both sides, estimated at the pilot
#' bandwidth), as in rd2d. The `cer*` variants shrink the MSE-optimal
#' bandwidth by rd2d's coverage-error rate factor `n^(1/(2p+4) - 1/(p+4))`.
#' Fuzzy designs use the linearised outcome `y / fs - itt / fs^2 * D`, as rd2d does.
#'
#' @param y Outcome, one value per observation in `coords`.
#' @param coords Output of [network_coordinates()].
#' @param p Polynomial order of the local fit.
#' @param deriv,tangvec,covs,kernel,cluster,fuzzy,weights As in [rdnet()].
#' @param bwselect `"mserd"` (one MSE-optimal bandwidth for both sides),
#'   `"msetwo"` (one per side), `"imserd"`/`"imsetwo"` (one bandwidth for all
#'   evaluation points), or the coverage-error-optimal `"cerrd"`, `"certwo"`,
#'   `"icerrd"`, `"icertwo"`.
#' @param vce `"hc0"`-`"hc3"` variance used inside the selector
#'   (`"cluster"` and `"conley"` select with `"hc0"` residuals, cluster-robust
#'   when `cluster` is given, as [rdnet()] does with those standard errors).
#' @param bwcheck Every bandwidth is at least the distance to the
#'   `bwcheck`-th nearest observation on each side (default `20 +` the number
#'   of order-`p + 1` terms).
#' @param scaleregul,scalebiascrct Weights on the bias-variance
#'   regularisation and the higher-order bias term.
#'
#' @return A data frame with one row per evaluation point: the
#'   control/treated bandwidths `h0`/`h1` (and `h` when they are equal), the
#'   pilots `dn`, `bn0`, `bn1`, the MSE constants (`V`, `B`, `Reg1`, `Reg2`
#'   per side, suffixed 0 for control and 1 for treated), the bounds
#'   `hmin0`/`hmin1`/`hmax0`/`hmax1`, the side sample sizes `n0`/`n1`, the
#'   local `basis` (`"s,o"`, or `"s"` where few boundary crossings leave no
#'   variation in the offset `o`), and a `status`.
#'
#'   [rdnet()] selects with `scaleregul = 3` (as rd2d does for estimation), so
#'   its bandwidths differ from this function's defaults; they are stored in
#'   the `bandwidths` element of its result.
#' @references Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). rd2d: Causal
#'   Inference in Boundary Discontinuity Designs. arXiv:2505.07989.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 7)
#' bw <- rdbwnet(city$obs$outcome, coords)
#' bw[, c("arc", "h0", "h1", "n0", "n1", "status")]
#' @export
rdbwnet <- function(y, coords, p = 1, deriv = c(0, 0), tangvec = NULL, covs = NULL,
                    kernel = "triangular", bwselect = "mserd", vce = "hc3", cluster = NULL,
                    fuzzy = NULL, weights = NULL, bwcheck = NULL, scaleregul = 1,
                    scalebiascrct = 1) {
  .check_coords(coords)
  bwselect <- tolower(bwselect)
  if (!bwselect %in% .BWSELECT) {
    stop(sprintf("bwselect must be one of %s", paste(.BWSELECT, collapse = ", ")), call. = FALSE)
  }
  if (!vce %in% .VCE_TYPES) {
    stop("vce must be one of 'hc0', 'hc1', 'hc2', 'hc3' (or 'cluster'/'conley', which select with 'hc0')",
         call. = FALSE)
  }
  if (!vce %in% c("hc0", "hc1", "hc2", "hc3")) vce <- "hc0"
  y <- .check_obs(y, coords$n_obs, "y")
  kfun <- get_kernel(kernel)
  kname <- if (is.character(kernel)) .kernel_name(kernel) else "triangular"
  om <- if (is.null(weights)) rep(1, length(y)) else .check_weights(weights, coords$n_obs)
  D <- if (is.null(fuzzy)) NULL else .check_obs(fuzzy, coords$n_obs, "fuzzy")
  cl <- if (is.null(cluster)) NULL else .check_cluster(cluster, coords$n_obs)
  Z <- NULL
  if (!is.null(covs)) {
    Z <- as.matrix(covs)
    storage.mode(Z) <- "double"
    .check_collinear(Z)
  }
  bwcheck <- if (is.null(bwcheck)) .default_min_obs(p + 1) else as.integer(bwcheck)
  T <- coords$treated
  s <- coords$s
  fin <- is.finite(s)
  # Effective sample sizes count distinct locations (mass points count once), as in rd2d.
  n_eff <- length(unique(coords$loc[fin]))
  n_eff_side <- c(length(unique(coords$loc[fin & !T])), length(unique(coords$loc[fin & T])))
  deriv <- as.integer(deriv)
  if (length(deriv) != 2 || min(deriv) < 0 || sum(deriv) > p) {
    stop("deriv must be c(a, c) with a, c >= 0 and a + c <= p", call. = FALSE)
  }
  J <- coords$n_eval
  tang <- if (is.null(tangvec)) NULL else .tangent_matrix(tangvec, J)
  v_ord <- if (!is.null(tang)) 1 else sum(deriv)
  expo_b <- 1 / (2 * p + 6)

  point <- function(j) {
    rec <- point_dim(j, 2)
    if (rec$status == "singular" && deriv[2] == 0 && is.null(tang)) {
      # Few boundary crossings near this point (a river, rail line or highway
      # crossing the border): the offset o barely varies, so fit in s alone.
      rec1 <- point_dim(j, 1)
      if (rec1$status == "ok") return(rec1)
    }
    rec
  }

  point_dim <- function(j, dim) {
    rec <- list(status = "ok", basis = if (dim == 2) "s,o" else "s")
    if (coords$eval_masked[j]) {
      rec$status <- "masked"
      return(rec)
    }
    vec <- rep(0, .n_terms(p, dim))
    vec[.exp_index(p, deriv[1], deriv[2], dim)] <- factorial(deriv[1]) * factorial(deriv[2])
    if (!is.null(tang)) {
      vec[] <- 0
      vec[.exp_index(p, 1, 0, dim)] <- tang[j, 1]
      vec[.exp_index(p, 0, 1, dim)] <- tang[j, 2]
    }
    at <- coords_at(coords, j)
    d <- at$d
    o <- at$o
    ok <- is.finite(d) & is.finite(o) & is.finite(s) & is.finite(y) & om > 0
    if (!is.null(D)) ok <- ok & is.finite(D)
    if (!is.null(Z)) ok <- ok & rowSums(!is.finite(Z)) == 0
    masks <- list(ok & !T, ok & T)
    # bwcheck counts distinct locations (repeated measurements count once), as
    # rdrobust and rd2d count unique observations
    dloc <- lapply(masks, function(m) d[m][!duplicated(coords$loc[m])])
    if (min(lengths(dloc)) < bwcheck) {
      rec$status <- "insufficient"
      return(rec)
    }
    sides <- lapply(masks, function(m) {
      .side(s[m], o[m], d[m], y[m], om[m], kfun, if (is.null(cl)) NULL else cl[m], dim)
    })
    dn <- .rot(s[ok], o[ok], kname, n_eff, dim)
    # nudged up so the bwcheck-th observation gets positive kernel weight
    bounds <- lapply(1:2, function(t) c(.kth(dloc[[t]], bwcheck) * (1 + 1e-9), max(sides[[t]]$d)))
    dns <- vapply(bounds, function(b) min(max(dn, b[1]), b[2]), numeric(1))
    thr <- vapply(sides, function(sd) stats::median(sd$d), numeric(1))
    res <- tryCatch({
      Dside <- lapply(masks, function(m) if (is.null(D)) NULL else D[m])
      if (!is.null(Z)) {
        zr <- apply(Z[ok, , drop = FALSE], 2, function(z) diff(range(z)))
        zc <- zr > 1e-12
        if (any(zc)) {
          Zs <- lapply(masks, function(m) Z[m, zc, drop = FALSE])
          g_y <- .covariate_gamma(sides, Zs, dns, p, lapply(sides, function(sd) sd$y))
          for (t in 1:2) sides[[t]]$y <- sides[[t]]$y - drop(Zs[[t]] %*% g_y)
          if (!is.null(D)) {
            g_d <- .covariate_gamma(sides, Zs, dns, p, Dside)
            Dside <- lapply(1:2, function(t) Dside[[t]] - drop(Zs[[t]] %*% g_d))
          }
        }
      }
      ybw <- list(NULL, NULL)
      if (!is.null(D)) {
        itt <- 0
        fs <- 0
        for (t in 1:2) {
          sgn <- if (t == 2) 1 else -1
          itt <- itt + sgn * sum(vec * .side_lm(sides[[t]], dns[t], p)$beta)
          fs <- fs + sgn * sum(vec * .side_lm(sides[[t]], dns[t], p, y = Dside[[t]])$beta)
        }
        if (is.finite(fs) && abs(fs) > 1e-8) {
          ybw <- lapply(1:2, function(t) sides[[t]]$y / fs - itt / fs^2 * Dside[[t]])
        }
      }
      hc <- list()
      bn <- numeric(2)
      for (t in 1:2) {
        sd <- sides[[t]]
        vq <- .side_coeff(sd, dns[t], p, vec)
        bc <- .side_consts(sd, p + 1, vq, dns[t], thr[t], NULL, vce, y = ybw[[t]])
        if (!(bc[["V"]] > 0 && bc[["Reg1"]] >= 0 && is.finite(bc[["B"]]))) {
          stop(.singular("invalid bias-pilot constants"))
        }
        b_t <- ((2 + 2 * (p + 1)) * bc[["V"]] / (2 * (bc[["B"]]^2 + scaleregul * bc[["Reg1"]])))^expo_b
        b_t <- min(max(b_t, bounds[[t]][1]), bounds[[t]][2])
        bn[t] <- b_t
        hc[[t]] <- .side_consts(sd, p, vec, dns[t], b_t, thr[t], vce, y = ybw[[t]])
        if (!(hc[[t]][["V"]] > 0 && hc[[t]][["Reg1"]] >= 0 && is.finite(hc[[t]][["B"]]))) {
          stop(.singular("invalid MSE constants"))
        }
      }
      list(hc = hc, bn = bn)
    }, error = function(e) NULL)
    if (is.null(res)) {
      rec$status <- "singular"
      return(rec)
    }
    rec <- c(rec, list(dn = dn, bn0 = res$bn[1], bn1 = res$bn[2], n0 = sum(masks[[1]]),
                       n1 = sum(masks[[2]]), hmin0 = bounds[[1]][1], hmin1 = bounds[[2]][1],
                       hmax0 = bounds[[1]][2], hmax1 = bounds[[2]][2]))
    for (t in 1:2) {
      h <- res$hc[[t]]
      rec[paste0(names(h), t - 1)] <- as.list(h)
    }
    rec
  }

  rows <- lapply(seq_len(J), point)
  cols <- c("dn", "bn0", "bn1", "n0", "n1", "hmin0", "hmin1", "hmax0", "hmax1",
            "V0", "B0", "Reg10", "Reg20", "V1", "B1", "Reg11", "Reg21")
  out <- data.frame(eval = seq_len(J), arc = coords$eval_arc,
                    status = vapply(rows, `[[`, "", "status"),
                    basis = vapply(rows, `[[`, "", "basis"), stringsAsFactors = FALSE)
  for (cn in cols) {
    out[[cn]] <- vapply(rows, function(r) if (is.null(r[[cn]])) NA_real_ else as.numeric(r[[cn]]), numeric(1))
  }
  out$h0 <- NA_real_
  out$h1 <- NA_real_
  good <- out$status == "ok"
  if (!any(good)) {
    out$h <- NA_real_
    return(out)
  }
  g <- out[good, , drop = FALSE]
  reg2 <- function(t) {
    r <- g[[paste0("Reg2", t)]]
    r[is.na(r)] <- 0
    r
  }
  Bc <- list(g$B0 + scalebiascrct * reg2(0), g$B1 + scalebiascrct * reg2(1))
  expo <- 1 / (2 * p + 4)
  base <- switch(bwselect, cerrd = "mserd", certwo = "msetwo", icerrd = "imserd",
                 icertwo = "imsetwo", bwselect)
  cer <- bwselect != base
  num_c <- 2 + 2 * v_ord
  den_c <- 2 * p + 2 - 2 * v_ord
  if (base %in% c("mserd", "imserd")) {
    num <- num_c * (g$V0 + g$V1)
    den <- den_c * ((Bc[[1]] - Bc[[2]])^2 + scaleregul * (g$Reg10 + g$Reg11))
    h <- if (base == "imserd") (mean(num) / mean(den))^expo else (num / den)^expo
    if (cer) h <- h * .cer_factor(n_eff, p)
    h <- pmin(pmax(h, pmax(g$hmin0, g$hmin1)), pmax(g$hmax0, g$hmax1))
    out$h0[good] <- h
    out$h1[good] <- h
  } else {
    for (t in 0:1) {
      num <- num_c * g[[paste0("V", t)]]
      den <- den_c * (Bc[[t + 1]]^2 + scaleregul * g[[paste0("Reg1", t)]])
      h <- if (base == "imsetwo") (mean(num) / mean(den))^expo else (num / den)^expo
      if (cer) h <- h * .cer_factor(n_eff_side[t + 1], p)
      out[[paste0("h", t)]][good] <- pmin(pmax(h, g[[paste0("hmin", t)]]), g[[paste0("hmax", t)]])
    }
  }
  out$h <- ifelse(.isclose(out$h0, out$h1), out$h0, NA_real_)
  out
}

# Directions (t_s, t_o) as a J x 2 matrix (one shared direction or one per point).
.tangent_matrix <- function(tangvec, J) {
  tv <- as.matrix(tangvec)
  if (length(tangvec) == 2 && is.null(dim(tangvec))) tv <- matrix(as.numeric(tangvec), 1, 2)
  if (ncol(tv) != 2 || !(nrow(tv) %in% c(1, J))) {
    stop("tangvec must be c(t_s, t_o) or a J x 2 matrix", call. = FALSE)
  }
  tv[rep_len(seq_len(nrow(tv)), J), , drop = FALSE]
}
