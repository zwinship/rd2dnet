# Network-constrained boundary average treatment effect curves (BATEC).
#
# At each boundary evaluation point b the estimator fits, separately on the
# treated and control sides, a weighted local polynomial of order p in the
# network coordinates (s_i, o_i(b)) with kernel weights K(d_i(b) / h) (times
# optional intensity weights). tau(b) is the difference of the two intercepts;
# in fuzzy designs it is the ratio of the outcome and take-up jumps. Robust
# bias-corrected inference uses the order q = p + 1 fit, as in rd2d.

#' Boundary discontinuity estimation on a network
#'
#' Estimates the boundary average treatment effect curve (BATEC) at each
#' evaluation point along the boundary, with robust bias-corrected inference
#' and a uniform confidence band. This carries rd2d's location-based
#' estimator over to network coordinates: `(s, o(b))` replace
#' `(X1 - b1, X2 - b2)` and the network distance `d(b)` replaces the Euclidean
#' distance in the kernel.
#'
#' @param y Outcome, one value per observation in `coords`. Missing values
#'   are skipped.
#' @param coords Output of [network_coordinates()].
#' @param h Bandwidth in metric units: a scalar, `c(h_control, h_treated)`
#'   (a length-2 vector always means this), one value per evaluation point, a
#'   `J x 2` matrix (use this for per-point bandwidths when `J = 2`), or
#'   `NULL` to select with [rdbwnet()] using `bwselect`.
#' @param p,q Polynomial order for estimation and for bias-corrected
#'   inference (default `p + 1`).
#' @param deriv `c(a, c)`: estimate the jump in `d^(a+c) mu / ds^a do^c`
#'   instead of the level. `c(1, 0)` is the jump in the slope across the
#'   boundary (a kink), `c(0, 1)` the jump in the slope along it. Requires
#'   `a + c <= p`.
#' @param tangvec Direction `c(t_s, t_o)` (or a `J x 2` matrix, one per
#'   evaluation point): estimate the jump in the directional derivative
#'   `t_s d/ds + t_o d/do`.
#' @param covs Covariates (vector or `n x m` matrix) entering the local fit
#'   with coefficients common to both sides (rd2d's `covs_eff`). Covariates
#'   that are constant in a local sample are dropped there.
#' @param kernel `"triangular"`, `"epanechnikov"`, `"uniform"`, `"gaussian"`
#'   or a function (see [get_kernel()]).
#' @param bwselect Bandwidth selector used when `h` is `NULL`; see
#'   [rdbwnet()]. Selection uses `scaleregul = 3`, as rd2d does for
#'   estimation, so it differs from a call to [rdbwnet()] with its defaults;
#'   the selected bandwidths are stored in the `bandwidths` element.
#' @param vce `"hc0"`-`"hc3"` (default `"hc3"`); `"cluster"` (with `cluster`) and
#'   `"conley"` (with `conley_cutoff`) are implied when those are given. rd2d uses
#'   HC1; HC3 is the default here because local samples on networks often contain
#'   high-leverage points (few distinct offsets), where HC1 understates the variance.
#' @param cluster Cluster labels, one per observation; implies `vce = "cluster"`
#'   (CR1 with rd2d's `(n - 1)/(n - k) * g/(g - 1)` factor). Points with fewer
#'   than two clusters on a side get status `"few clusters"`. With repeated
#'   measurements (panels), clusters much smaller than the bandwidth, such as
#'   the unit itself, can still understate the variance of spatially
#'   correlated shocks; `conley_cutoff` is the safer choice.
#' @param conley_cutoff,conley_distance Spatial HAC with Bartlett weights over
#'   `"network"` (default) or `"euclidean"` distance, used for pointwise
#'   errors and for the uniform band.
#' @param conley_resid_scale,conley_resid_order Conley residuals come from an
#'   auxiliary fit of order `conley_resid_order` (default `q`) at
#'   `conley_resid_scale * h` (default 3). Residuals from the estimation fit
#'   itself absorb spatially correlated shocks at scales below `h`, which
#'   makes the HAC variance far too small.
#' @param fuzzy Treatment take-up `D` (0/1 or a share). The estimand becomes
#'   the ratio of the outcome jump to the take-up jump at each boundary point.
#' @param weights Non-negative intensity weights (e.g. floor area)
#'   multiplying the kernel weights.
#' @param pooled Fit one regression `y ~ 1 + T + poly(s, o)` with common
#'   slopes instead of side-specific polynomials (levels only).
#' @param level Confidence level in percent (e.g. 95).
#' @param cbands Add the uniform confidence band `cb_lo`/`cb_hi` to the table.
#'   The joint covariance and its critical value are stored either way.
#' @param repp Gaussian draws used to simulate the uniform critical value.
#' @param min_obs Minimum local observations per side (default `20 +` the
#'   number of order-`q` terms); points with fewer are reported as
#'   `"insufficient"`.
#' @param scaleregul,scalebiascrct Bandwidth-selection regularisation and
#'   higher-order bias weights, passed to [rdbwnet()] when `h` is selected.
#' @param seed Seed for the critical-value simulation (default 0, so bands and
#'   LBATE are reproducible; `NULL` draws from the current stream). The global
#'   random number stream is left untouched.
#'
#' @return An object of class `rdnet`, a list with
#'   \describe{
#'     \item{`table`}{one row per evaluation point: `eval`, `arc`, `node_arc`,
#'       bandwidths `h0`, `h1` (and `h` when equal), local sample sizes `n0`,
#'       `n1`, the order-`p` estimate `tau` and its standard error `se`, the
#'       bias-corrected `tau_bc` with robust `se_rb`, the robust interval
#'       `ci_lo`/`ci_hi`, `pvalue`, a `status` (`"ok"`, `"masked"`,
#'       `"insufficient"` (fewer than `min_obs` observations on a side),
#'       `"few clusters"`, `"singular"`, `"weak first stage"`, `"no
#'       bandwidth"`), the local `basis` (`"s,o"` or `"s"`, see Details), and
#'       optionally the uniform band `cb_lo`/`cb_hi`. Fuzzy designs add `itt`,
#'       `fs` and `fs_se`.}
#'     \item{`cov`}{the `J x J` covariance of `tau_bc` across evaluation
#'       points (rd2d's `params_cov`), `NA` for points not estimated; its
#'       diagonal is `se_rb^2`.}
#'     \item{`crit`}{the uniform critical value implied by `cov`.}
#'     \item{`bandwidths`}{the [rdbwnet()] output when `h` was selected.}
#'   }
#'   plus `p`, `q`, `kernel`, `vce`, `level`, `metric`, `fuzzy` and `coords`.
#'   Use [wbate()] and [lbate()] for aggregate effects, `summary()`, `plot()`,
#'   `coef()`, `vcov()` and `confint()` for the usual summaries.
#' @details The local fit is a polynomial in `(s, o)`. Where the network
#'   meets the boundary at only a few points (rivers, rail lines, highways
#'   crossing a border), the offset `o` hardly varies near an evaluation point
#'   and that fit is singular; the point is then fitted in `s` alone and its
#'   `basis` column reads `"s"` instead of `"s,o"`.
#' @references
#' Cattaneo, M. D., Titiunik, R. and Yu, R. R. (2026). rd2d: Causal Inference
#' in Boundary Discontinuity Designs. arXiv:2505.07989.
#'
#' Calonico, S., Cattaneo, M. D. and Titiunik, R. (2014). Robust Nonparametric
#' Confidence Intervals for Regression-Discontinuity Designs. Econometrica
#' 82(6): 2295-2326.
#'
#' Conley, T. G. (1999). GMM Estimation with Cross Sectional Dependence.
#' Journal of Econometrics 92(1): 1-45.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 9)
#' res <- rdnet(city$obs$outcome, coords, weights = city$obs$far, seed = 1)
#' res
#' ok <- res$table$status == "ok"
#' wbate(res, as.numeric(ok))
#' plot(res)
#' @export
rdnet <- function(y, coords, h = NULL, p = 1, q = NULL, deriv = c(0, 0), tangvec = NULL,
                  covs = NULL, kernel = "triangular", bwselect = "mserd", vce = "hc3",
                  cluster = NULL, conley_cutoff = NULL, conley_distance = "network",
                  conley_resid_scale = 3, conley_resid_order = NULL, fuzzy = NULL,
                  weights = NULL, pooled = FALSE, level = 95, cbands = TRUE, repp = 2000,
                  min_obs = NULL, scaleregul = 3, scalebiascrct = 1, seed = 0) {
  .check_coords(coords)
  y <- .check_obs(y, coords$n_obs, "y")
  if (!is.numeric(level) || length(level) != 1 || !(level > 1 && level < 100)) {
    stop(sprintf("level is a percentage, e.g. 95 (got %s)", format(level)), call. = FALSE)
  }
  q <- if (is.null(q)) p + 1 else as.integer(q)
  if (p < 0 || q < p) stop("need 0 <= p <= q", call. = FALSE)
  J <- coords$n_eval
  deriv <- as.integer(deriv)
  if (length(deriv) != 2 || min(deriv) < 0 || sum(deriv) > p) {
    stop("deriv must be c(a, c) with a, c >= 0 and a + c <= p", call. = FALSE)
  }
  tang <- NULL
  if (!is.null(tangvec)) {
    if (sum(deriv) != 0) stop("give either deriv or tangvec, not both", call. = FALSE)
    if (p < 1) stop("tangvec needs p >= 1", call. = FALSE)
    tang <- .tangent_matrix(tangvec, J)
  }
  if (pooled && (sum(deriv) > 0 || !is.null(tang))) {
    stop("pooled = TRUE estimates level jumps only", call. = FALSE)
  }
  kfun <- get_kernel(kernel)
  kname <- if (is.character(kernel)) .kernel_name(kernel) else "custom"
  intensity <- if (is.null(weights)) rep(1, length(y)) else .check_weights(weights, length(y))
  D <- if (is.null(fuzzy)) NULL else .check_obs(fuzzy, length(y), "fuzzy")
  Z <- NULL
  if (!is.null(covs)) {
    Z <- as.matrix(covs)
    storage.mode(Z) <- "double"
    if (nrow(Z) != length(y)) stop("covs must have one row per observation", call. = FALSE)
    .check_collinear(Z)
  }
  if (!is.null(cluster)) {
    vce <- "cluster"
    cluster <- .check_cluster(cluster, length(y))
  }
  if (!is.null(conley_cutoff)) {
    if (!is.null(cluster)) {
      warning("both cluster and conley_cutoff given: standard errors use Conley; cluster only ",
              "affects bandwidth selection", call. = FALSE)
    }
    vce <- "conley"
  }
  if (!vce %in% .VCE_TYPES) {
    stop(sprintf("vce must be one of %s", paste(.VCE_TYPES, collapse = ", ")), call. = FALSE)
  }
  if (vce == "cluster" && is.null(cluster)) stop("vce = 'cluster' needs a cluster variable", call. = FALSE)
  if (vce == "conley" && (is.null(conley_cutoff) || conley_cutoff <= 0)) {
    stop("vce = 'conley' needs a positive conley_cutoff", call. = FALSE)
  }
  if (!conley_distance %in% c("network", "euclidean")) {
    stop("conley_distance must be 'network' or 'euclidean'", call. = FALSE)
  }
  base_vce <- if (vce %in% c("hc1", "hc2", "hc3")) vce else "hc0"
  min_obs <- if (is.null(min_obs)) .default_min_obs(q) else as.integer(min_obs)
  aux_order <- if (is.null(conley_resid_order)) q else as.integer(conley_resid_order)

  bw <- NULL
  if (is.null(h)) {
    bw <- rdbwnet(y, coords, p = p, deriv = deriv, tangvec = tangvec, covs = Z, kernel = kernel,
                  bwselect = bwselect, vce = base_vce, cluster = cluster, fuzzy = D,
                  weights = weights, bwcheck = min_obs, scaleregul = scaleregul,
                  scalebiascrct = scalebiascrct)
    H <- cbind(bw$h0, bw$h1)
    bases <- bw$basis
  } else {
    H <- .bandwidth_matrix(h, J)
    bases <- rep(NA_character_, J)
  }

  z <- stats::qnorm(0.5 + level / 200)
  T <- coords$treated
  Tf <- as.numeric(T)
  s <- coords$s
  finite <- is.finite(s) & is.finite(y)
  if (!is.null(D)) finite <- finite & is.finite(D)
  if (!is.null(Z)) finite <- finite & rowSums(!is.finite(Z)) == 0
  targets <- if (is.null(D)) list(y) else list(y, D)
  conley <- vce == "conley"

  tab <- data.frame(eval = seq_len(J), arc = coords$eval_arc, node_arc = coords$eval_node_arc,
                    h0 = H[, 1], h1 = H[, 2], n0 = 0L, n1 = 0L, tau = NA_real_, se = NA_real_,
                    tau_bc = NA_real_, se_rb = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                    pvalue = NA_real_, status = "ok", basis = "", stringsAsFactors = FALSE)
  if (!is.null(D)) {
    tab$itt <- NA_real_
    tab$fs <- NA_real_
    tab$fs_se <- NA_real_
  }
  psi_store <- list()
  # Conley: all variances come from one HAC pass over every local sample at the end.
  hac_cols <- list()

  local_var <- function(idx, psi, k) {
    if (conley) {
      return(NA_real_)
    }
    v <- .variance(psi, vce, cluster = if (is.null(cluster)) NULL else cluster[idx])
    # .variance() applies g / (g - 1); rd2d's CR1 also has (n - 1) / (n - k)
    if (vce == "cluster") v <- v * .cluster_dof(length(idx), k, cluster[idx])
    v
  }

  g_min <- integer() # clusters per side at each point (fewest side)
  for (j in seq_len(J)) {
    if (coords$eval_masked[j]) {
      tab$status[j] <- "masked"
      next
    }
    if (!(all(is.finite(H[j, ])) && all(H[j, ] > 0))) {
      # the selector's reason (insufficient, singular) when it has one
      tab$status[j] <- if (!is.null(bw) && bw$status[j] != "ok") bw$status[j] else "no bandwidth"
      next
    }
    at <- coords_at(coords, j)
    d <- at$d
    o <- at$o
    base <- is.finite(d) & is.finite(o) & finite & intensity > 0
    hside <- ifelse(T, H[j, 2], H[j, 1])
    # each side's kernel weights scaled by 1 / h_side^2, as in rd2d: this matters
    # for the coefficients shared across sides (covariates, pooled slopes)
    w <- kfun(d / hside) / hside^2 * intensity
    use <- w > 0 & base
    i1 <- which(use & T)
    i0 <- which(use & !T)
    tab$n1[j] <- length(i1)
    tab$n0[j] <- length(i0)
    # count distinct locations: repeated measurements at one site add no design points
    if (min(length(unique(coords$loc[i1])), length(unique(coords$loc[i0]))) < min_obs) {
      tab$status[j] <- "insufficient"
      next
    }
    idx <- c(i1, i0)
    zcols <- NULL
    if (!is.null(Z)) {
      zr <- apply(Z[idx, , drop = FALSE], 2, function(v) diff(range(v)))
      zcols <- zr > 1e-12 # drop locally constant covariates
      if (!any(zcols)) zcols <- NULL
    }
    Zi <- if (is.null(zcols)) NULL else Z[idx, zcols, drop = FALSE]

    if (vce == "cluster") {
      g_side <- min(length(unique(cluster[i1])), length(unique(cluster[i0])))
      g_min <- c(g_min, g_side)
      if (g_side < 2) {
        tab$status[j] <- "few clusters"
        next
      }
    }

    fit_all <- function() {
      dim <- 2
      aux <- NULL
      if (conley) {
        hs <- hside * conley_resid_scale
        wk <- kfun(d / hs) / hs^2 * intensity
        m <- which(wk > 0 & base)
        Xa <- .design(s[m] / hs[m], o[m] / hs[m], Tf[m], aux_order,
                      if (is.null(zcols)) NULL else Z[m, zcols, drop = FALSE], FALSE, dim)
        Xi <- .design(s[idx] / hs[idx], o[idx] / hs[idx], Tf[idx], aux_order, Zi, FALSE, dim)
        ka <- .identifiable(Xa, wk[m])
        aux <- list(Xa = Xa[, ka, drop = FALSE], wa = wk[m], m = m, Xi = Xi[, ka, drop = FALSE])
      }
      out <- list()
      ks <- c()
      for (order in unique(c(p, q))) {
        X <- .design(s[idx] / hside[idx], o[idx] / hside[idx], Tf[idx], order, Zi, pooled, dim)
        cvec <- rep(0, ncol(X))
        if (pooled) {
          cvec[2] <- 1
        } else {
          dv <- .deriv_vector(order, deriv, if (is.null(tang)) NULL else tang[j, ], dim)
          k <- length(dv$vec)
          cvec[seq_len(k)] <- dv$vec / H[j, 2]^dv$order
          cvec[k + seq_len(k)] <- -dv$vec / H[j, 1]^dv$order
        }
        # keep the identifiable columns (few boundary crossings leave too little
        # variation in the offset o for its higher powers, or any)
        keep <- .identifiable(X, w[idx])
        if (any(cvec[!keep] != 0)) stop(.singular())
        if (order == p) {
          ex <- .poly_exponents(order)
          oexp <- if (pooled) ex[-1, "c"] else rep(ex[, "c"], 2)
          first <- if (pooled) 2 else 0
          has_o <- any(keep[first + seq_along(oexp)][oexp > 0])
        }
        X <- X[, keep, drop = FALSE]
        cvec <- cvec[keep]
        jumps <- numeric(0)
        psis <- list()
        for (target in targets) {
          f <- .wls(X, target[idx], w[idx])
          jumps <- c(jumps, sum(cvec * f$beta))
          if (is.null(aux)) {
            psis[[length(psis) + 1]] <- .influence(f, cvec, base_vce)
          } else {
            beta_a <- .wls(aux$Xa, target[aux$m], aux$wa)$beta
            r_aux <- target[idx] - drop(aux$Xi %*% beta_a)
            psis[[length(psis) + 1]] <- .loadings(f, cvec) * r_aux
          }
        }
        if (is.null(D)) {
          tau <- jumps[1]
          psi <- psis[[1]]
          extra <- list()
        } else {
          itt <- jumps[1]
          fs <- jumps[2]
          if (abs(fs) < 1e-8) stop(structure(class = c("rd2dnet_weak", "error", "condition"),
                                             list(message = "weak first stage", call = NULL)))
          tau <- itt / fs
          psi <- (psis[[1]] - tau * psis[[2]]) / fs
          extra <- list(itt = itt, fs = fs, psi_fs = psis[[2]])
        }
        ks[as.character(order)] <- ncol(X)
        out[[as.character(order)]] <- list(tau = tau, se = sqrt(max(local_var(idx, psi, ncol(X)), 0)),
                                           psi = psi, extra = extra)
      }
      # "s" when no term in the offset o is identified (a river, rail line or
      # highway crossing the border at one point): the fit is in s alone
      list(out = out, ks = ks, basis = if (has_o) "s,o" else "s")
    }

    fit <- tryCatch(fit_all(), rd2dnet_singular = function(e) "singular",
                    rd2dnet_weak = function(e) "weak first stage")
    if (is.character(fit)) {
      tab$status[j] <- fit
      next
    }
    fp <- fit$out[[as.character(p)]]
    fq <- fit$out[[as.character(q)]]
    kp <- fit$ks[[as.character(p)]]
    kq <- fit$ks[[as.character(q)]]
    if (!is.null(D)) {
      tab$itt[j] <- fp$extra$itt
      tab$fs[j] <- fp$extra$fs
      tab$fs_se[j] <- sqrt(max(local_var(idx, fp$extra$psi_fs, kp), 0))
      if (conley) hac_cols[[length(hac_cols) + 1]] <- list(j = j, name = "fs_se", idx = idx, psi = fp$extra$psi_fs)
    }
    if (conley) {
      hac_cols[[length(hac_cols) + 1]] <- list(j = j, name = "se", idx = idx, psi = fp$psi)
      hac_cols[[length(hac_cols) + 1]] <- list(j = j, name = "se_rb", idx = idx, psi = fq$psi)
    }
    tab$basis[j] <- fit$basis
    tab$tau[j] <- fp$tau
    tab$se[j] <- fp$se
    tab$tau_bc[j] <- fq$tau
    tab <- .set_rb(tab, j, fq$se, z)
    psi_store[[as.character(j)]] <- list(j = j, idx = idx, psi = fq$psi, k = kq)
  }

  crit <- NULL
  cov <- NULL
  if (conley && length(hac_cols)) {
    union <- sort(unique(unlist(lapply(hac_cols, `[[`, "idx"))))
    pos <- integer(coords$n_obs)
    pos[union] <- seq_along(union)
    Psi <- matrix(0, length(union), length(hac_cols))
    for (cc in seq_along(hac_cols)) Psi[pos[hac_cols[[cc]]$idx], cc] <- hac_cols[[cc]]$psi
    V <- .hac_quadratic(coords, union, Psi, conley_cutoff, conley_distance)
    se_all <- sqrt(pmax(diag(V), 0))
    for (cc in seq_along(hac_cols)) {
      hcj <- hac_cols[[cc]]
      if (hcj$name == "se_rb") {
        tab <- .set_rb(tab, hcj$j, se_all[cc], z)
      } else {
        tab[[hcj$name]][hcj$j] <- se_all[cc]
      }
    }
    qc <- which(vapply(hac_cols, `[[`, "", "name") == "se_rb")
    qj <- vapply(hac_cols[qc], `[[`, 0L, "j")
    cov <- matrix(NA_real_, J, J)
    cov[qj, qj] <- V[qc, qc]
  } else if (length(psi_store)) {
    cov <- .joint_covariance(psi_store, coords, vce, cluster, J)
  }
  g_ok <- g_min[g_min >= 2]
  if (length(g_ok) && min(g_ok) < 10) {
    warning(sprintf(paste0("cluster-robust standard errors rest on few clusters: as few as %d per side within ",
                           "the bandwidth at some evaluation points; with fewer than about 10 per side they are ",
                           "too small (coverage near 0.7 in Monte Carlo runs with 5-6). Consider conley_cutoff."),
                    min(g_ok)), call. = FALSE)
  }
  tab$h <- ifelse(.isclose(tab$h0, tab$h1), tab$h0, NA_real_)
  if (!is.null(cov)) {
    ok <- is.finite(diag(cov))
    crit <- .uniform_critical_value(cov[ok, ok, drop = FALSE], level, repp, seed)
    if (cbands) {
      tab$cb_lo <- tab$tau_bc - crit * tab$se_rb
      tab$cb_hi <- tab$tau_bc + crit * tab$se_rb
    }
  }
  structure(list(table = tab, p = p, q = q, kernel = kname, vce = vce, level = level,
                 metric = format(coords$network$metric), fuzzy = !is.null(D), crit = crit,
                 cov = cov, bandwidths = bw, coords = coords),
            class = "rdnet")
}

# User bandwidths as a J x 2 matrix of (control, treated).
.bandwidth_matrix <- function(h, J) {
  if (!is.numeric(h) || anyNA(h) || any(h <= 0)) stop("bandwidths must be positive", call. = FALSE)
  if (is.matrix(h)) {
    if (all(dim(h) == c(J, 2))) {
      return(unname(h) + 0)
    }
    stop("h must be a scalar, c(h0, h1), one value per evaluation point, or a J x 2 matrix", call. = FALSE)
  }
  h <- as.numeric(h)
  if (length(h) == 1) {
    return(matrix(h, J, 2))
  }
  if (length(h) == 2) {
    return(matrix(h, J, 2, byrow = TRUE))
  }
  if (length(h) == J) {
    return(cbind(h, h, deparse.level = 0))
  }
  stop("h must be a scalar, c(h0, h1), one value per evaluation point, or a J x 2 matrix", call. = FALSE)
}

# Target vector over the order-`order` basis, and the derivative order it picks.
.deriv_vector <- function(order, deriv, tangvec_j, dim = 2) {
  vec <- rep(0, .n_terms(order, dim))
  if (!is.null(tangvec_j)) {
    vec[.exp_index(order, 1, 0, dim)] <- tangvec_j[1]
    vec[.exp_index(order, 0, 1, dim)] <- tangvec_j[2]
    return(list(vec = vec, order = 1))
  }
  vec[.exp_index(order, deriv[1], deriv[2], dim)] <- factorial(deriv[1]) * factorial(deriv[2])
  list(vec = vec, order = sum(deriv))
}

# Joint local design: side-specific polynomials (or pooled slopes) plus common covariates.
# Columns of X kept for the weighted fit: the largest identifiable sub-basis.
# Columns are taken in order (lower degree first; s before o within a degree);
# a column is dropped when its weighted residual after projecting on the columns
# already kept is below tol times its norm (as lm() drops aliased terms), or when
# it would give some observation leverage 1 (a coefficient determined by that
# observation alone).
.identifiable <- function(X, w, tol = 1e-4) {
  A <- X * sqrt(w)
  Q <- matrix(0, nrow(A), 0)
  keep <- logical(ncol(A))
  lev <- numeric(nrow(A))
  for (j in seq_len(ncol(A))) {
    v <- A[, j]
    n0 <- sqrt(sum(v * v))
    if (n0 == 0) next
    for (rep in 1:2) { # Gram-Schmidt, repeated for accuracy
      for (k in seq_len(ncol(Q))) v <- v - sum(Q[, k] * v) * Q[, k]
    }
    n1 <- sqrt(sum(v * v))
    if (n1 <= tol * n0) next
    u <- v / n1
    if (max(lev + u * u) > 1 - 1e-8) next
    keep[j] <- TRUE
    Q <- cbind(Q, u)
    lev <- lev + u * u
  }
  keep
}

.design <- function(us, uo, Ti, order, Zi, pooled, dim = 2) {
  R <- .poly_basis(us, uo, order, dim)
  X <- if (pooled) cbind(1, Ti, R[, -1, drop = FALSE]) else cbind(R * Ti, R * (1 - Ti))
  X <- unname(X)
  if (is.null(Zi)) X else cbind(X, Zi, deparse.level = 0)
}

# rd2d's (n - 1) / (n - k) cluster degrees-of-freedom factor (1 with one cluster or n <= k).
.cluster_dof <- function(n, k, cl) {
  if (n > k && length(unique(cl)) > 1) (n - 1) / (n - k) else 1
}

# Robust standard error, confidence interval and p-value from the bias-corrected fit.
.set_rb <- function(tab, j, se_q, z) {
  tau_q <- tab$tau_bc[j]
  tab$se_rb[j] <- se_q
  tab$ci_lo[j] <- tau_q - z * se_q
  tab$ci_hi[j] <- tau_q + z * se_q
  tab$pvalue[j] <- if (isTRUE(se_q > 0)) 2 * stats::pnorm(-abs(tau_q / se_q)) else NA_real_
  tab
}

# J x J covariance of the bias-corrected estimates, NA where a point was not
# estimated. Its diagonal reproduces the pointwise variances exactly; for
# clustering the small-sample factor uses each point's own local sample.
.joint_covariance <- function(psi_store, coords, vce, cluster, J) {
  keys <- vapply(psi_store, `[[`, 0L, "j")
  union <- sort(unique(unlist(lapply(psi_store, `[[`, "idx"))))
  pos <- integer(coords$n_obs)
  pos[union] <- seq_along(union)
  Psi <- matrix(0, length(union), length(keys))
  for (cc in seq_along(psi_store)) Psi[pos[psi_store[[cc]]$idx], cc] <- psi_store[[cc]]$psi
  if (vce == "cluster") {
    G <- rowsum(Psi, cluster[union])
    f <- vapply(psi_store, function(ps) {
      g <- length(unique(cluster[ps$idx]))
      sqrt((if (g > 1) g / (g - 1) else 1) * .cluster_dof(length(ps$idx), ps$k, cluster[ps$idx]))
    }, numeric(1))
    V <- crossprod(G) * outer(f, f)
  } else {
    V <- crossprod(Psi)
  }
  out <- matrix(NA_real_, J, J)
  out[keys, keys] <- V
  out
}

# Critical value for a uniform band: quantile of max |Z| with Z ~ N(0, corr).
.uniform_critical_value <- function(Sigma, level, repp, seed) {
  # one point: the normal quantile itself, no simulation
  if (nrow(Sigma) == 1) return(stats::qnorm(0.5 + level / 200))
  sd <- sqrt(pmax(diag(Sigma), 1e-300))
  corr <- Sigma / outer(sd, sd)
  ev <- eigen(corr, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(sd))
  draws <- .with_seed(seed, matrix(stats::rnorm(repp * length(sd)), repp, length(sd))) %*% t(L)
  unname(stats::quantile(apply(abs(draws), 1, max), level / 100))
}

#' Aggregate boundary treatment effects
#'
#' `wbate()` is the weighted boundary average treatment effect `w' tau`
#' (rd2d's WBATE): `weights` has one entry per evaluation point and is
#' normalised to sum to one. By default every distinct place where an estimated
#' point lies gets equal weight: evaluation points that snapped to the same place
#' (`node_arc`) repeat one estimate and share its weight. Points that were not
#' estimated must get weight 0. Its standard error `sqrt(w' Sigma w)` uses the
#' joint covariance `x$cov`, so it accounts for the correlation between nearby points.
#'
#' `lbate()` is the largest boundary average treatment effect
#' `max_j tau(b_j)` (rd2d's LBATE), with interval
#' `[max_j (tau_bc - crit se_rb), max_j (tau_bc + crit se_rb)]` over the
#' estimated points, `crit` the uniform critical value, so it inherits the
#' uniform band's joint coverage. `kind = "min"` gives the smallest effect
#' (`[min_j lo_j, min_j hi_j]`) and `kind = "absmax"` the largest absolute effect
#' (`[max_j dist(0, [lo_j, hi_j]), max_j max(|lo_j|, |hi_j|)]`).
#'
#' @param x An [rdnet()] result.
#' @param weights One weight per evaluation point (default: equal weights on
#'   the distinct places of the estimated points).
#' @param kind For `lbate()`: `"max"` (default), `"min"` or `"absmax"`.
#' @return A named numeric vector: `tau`, `tau_bc`, `se_rb`, `ci_lo`,
#'   `ci_hi`, `pvalue` (the bias-corrected estimate drives the interval and
#'   p-value; `se_rb` and `pvalue` are `NA` for LBATE).
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 9)
#' res <- rdnet(city$obs$outcome, coords, h = 600, seed = 1)
#' ok <- res$table$status == "ok"
#' wbate(res, as.numeric(ok))
#' lbate(res)
#' @export
wbate <- function(x, weights = NULL) {
  if (!inherits(x, "rdnet")) stop("expected an rdnet() result", call. = FALSE)
  J <- nrow(x$table)
  w <- if (is.null(weights)) .place_weights(x$table) else as.numeric(weights)
  if (length(w) != J) {
    stop(sprintf("WBATE weights must have one entry per evaluation point (%d)", J), call. = FALSE)
  }
  if (!all(is.finite(w))) stop("WBATE weights must be finite", call. = FALSE)
  if (sum(w) == 0) stop("WBATE weights must have a nonzero sum", call. = FALSE)
  ok <- x$table$status == "ok"
  bad <- which(!ok & w != 0)
  if (length(bad)) {
    stop(sprintf("WBATE puts weight on evaluation points that were not estimated: %s; give them weight 0",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  w <- (w / sum(w))[ok]
  t <- x$table[ok, , drop = FALSE]
  tau_q <- sum(w * t$tau_bc)
  se <- sqrt(max(drop(w %*% x$cov[ok, ok, drop = FALSE] %*% w), 0))
  z <- stats::qnorm(0.5 + x$level / 200)
  c(tau = sum(w * t$tau), tau_bc = tau_q, se_rb = se, ci_lo = tau_q - z * se,
    ci_hi = tau_q + z * se, pvalue = if (se > 0) 2 * stats::pnorm(-abs(tau_q / se)) else NA_real_)
}

#' @rdname wbate
#' @export
lbate <- function(x, kind = c("max", "min", "absmax")) {
  if (!inherits(x, "rdnet")) stop("expected an rdnet() result", call. = FALSE)
  kind <- match.arg(kind)
  ok <- x$table$status == "ok"
  if (!any(ok) || is.null(x$crit)) stop("LBATE needs at least one estimated evaluation point", call. = FALSE)
  t <- x$table[ok, , drop = FALSE]
  lo <- t$tau_bc - x$crit * t$se_rb
  hi <- t$tau_bc + x$crit * t$se_rb
  v <- switch(kind,
    max = c(max(t$tau), max(t$tau_bc), max(lo), max(hi)),
    min = c(min(t$tau), min(t$tau_bc), min(lo), min(hi)),
    absmax = c(max(abs(t$tau)), max(abs(t$tau_bc)), max(pmax(0, lo, -hi)), max(pmax(abs(lo), abs(hi)))))
  c(tau = v[1], tau_bc = v[2], se_rb = NA_real_, ci_lo = v[3], ci_hi = v[4], pvalue = NA_real_)
}

# Equal weight per distinct place among the estimated points, shared by repeats.
.place_weights <- function(tab) {
  ok <- tab$status == "ok"
  w <- numeric(nrow(tab))
  if (!any(ok)) return(w)
  place <- round(tab$node_arc[ok], 9)
  w[ok] <- 1 / as.numeric(table(place)[as.character(place)])
  w
}

#' Summaries of rdnet results
#'
#' @param object,x An [rdnet()] result (or its summary).
#' @param wbate Optional WBATE weights (see [wbate()]) to add an aggregate row.
#' @param lbate Add the LBATE row (see [lbate()]).
#' @param digits Digits to print.
#' @param parm Evaluation points (indices) for `confint()`; default all.
#' @param level Confidence level for `confint()` as a proportion.
#' @param type For `coef()`: `"tau_bc"` (bias-corrected, default) or `"tau"`.
#' @param band Draw the uniform band.
#' @param true_tau Optional true effects at the evaluation points, drawn dashed: one
#'   value per point, a single value, or a function of the arc position.
#' @param ... Further arguments (passed to `plot()` for the plot method).
#' @return `summary()` returns an object of class `summary.rdnet`; `coef()` the
#'   estimates by evaluation point; `vcov()` the joint covariance `x$cov`;
#'   `confint()` robust intervals `tau_bc -/+ z se_rb`; `as.data.frame()` the
#'   result table; `plot()` draws the BATEC curve with pointwise robust
#'   intervals and the uniform band and returns `x` invisibly.
#' @name rdnet-methods
NULL

#' @rdname rdnet-methods
#' @export
summary.rdnet <- function(object, wbate = NULL, lbate = FALSE, ...) {
  agg <- NULL
  if (!is.null(wbate)) agg <- rbind(agg, WBATE = rd2dnet::wbate(object, wbate))
  if (lbate) agg <- rbind(agg, LBATE = rd2dnet::lbate(object))
  structure(list(object = object, aggregate = agg), class = "summary.rdnet")
}

#' @rdname rdnet-methods
#' @export
print.summary.rdnet <- function(x, digits = 3, ...) {
  r <- x$object
  t <- r$table
  cat(sprintf("rd2dnet: network boundary discontinuity estimates (%s)\n", if (r$fuzzy) "fuzzy" else "sharp"))
  cat(sprintf("  metric=%s  kernel=%s  p=%d  q=%d  vce=%s\n", r$metric, r$kernel, r$p, r$q, r$vce))
  cat(sprintf("  evaluation points: %d  estimated: %d\n", nrow(t), sum(t$status == "ok")))
  if ("cb_lo" %in% names(t)) cat(sprintf("  uniform band critical value: %.3f\n", r$crit))
  cols <- c("arc", "h0", "h1", "n0", "n1", "tau", "se", "tau_bc", "ci_lo", "ci_hi", "pvalue", "status")
  if (r$fuzzy) cols <- append(cols, "fs", after = 5)
  print(format(t[cols], digits = digits, nsmall = digits), row.names = FALSE)
  if (!is.null(x$aggregate)) {
    cat("\naggregate effects (tau_bc and CI are bias-corrected):\n")
    print(round(x$aggregate, digits), na.print = "")
  }
  invisible(x)
}

#' @rdname rdnet-methods
#' @export
print.rdnet <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @rdname rdnet-methods
#' @export
coef.rdnet <- function(object, type = c("tau_bc", "tau"), ...) {
  type <- match.arg(type)
  stats::setNames(object$table[[type]], object$table$eval)
}

#' @rdname rdnet-methods
#' @export
vcov.rdnet <- function(object, ...) {
  v <- object$cov
  if (!is.null(v)) dimnames(v) <- list(object$table$eval, object$table$eval)
  v
}

#' @rdname rdnet-methods
#' @export
confint.rdnet <- function(object, parm, level = object$level / 100, ...) {
  t <- object$table
  if (missing(parm)) parm <- seq_len(nrow(t))
  z <- stats::qnorm(0.5 + level / 2)
  out <- cbind(t$tau_bc - z * t$se_rb, t$tau_bc + z * t$se_rb)[parm, , drop = FALSE]
  dimnames(out) <- list(t$eval[parm], sprintf("%g %%", 100 * c((1 - level) / 2, 1 - (1 - level) / 2)))
  out
}

#' @rdname rdnet-methods
#' @param row.names,optional Ignored.
#' @export
as.data.frame.rdnet <- function(x, row.names = NULL, optional = FALSE, ...) x$table

#' @rdname rdnet-methods
#' @export
plot.rdnet <- function(x, band = TRUE, true_tau = NULL, ...) {
  t <- x$table
  if (is.function(true_tau)) true_tau <- true_tau(t$arc)
  if (!is.null(true_tau)) true_tau <- rep_len(as.numeric(true_tau), nrow(t))
  rng <- range(c(t$ci_lo, t$ci_hi, if (band && "cb_lo" %in% names(t)) c(t$cb_lo, t$cb_hi), t$tau,
                 true_tau, 0), na.rm = TRUE, finite = TRUE)
  args <- utils::modifyList(list(x = t$arc, y = t$tau, type = "n", ylim = rng,
                                 xlab = "position along boundary", ylab = "treatment effect"),
                            list(...))
  do.call(graphics::plot, args)
  graphics::abline(h = 0, col = "grey60", lwd = 0.8)
  if (band && "cb_lo" %in% names(t)) {
    ok <- is.finite(t$cb_lo) & is.finite(t$cb_hi)
    runs <- split(which(ok), cumsum(!ok)[ok])
    for (r in runs) {
      graphics::polygon(c(t$arc[r], rev(t$arc[r])), c(t$cb_lo[r], rev(t$cb_hi[r])),
                        col = grDevices::adjustcolor("steelblue", 0.15), border = NA)
    }
  }
  ok <- is.finite(t$ci_lo)
  graphics::arrows(t$arc[ok], t$ci_lo[ok], t$arc[ok], t$ci_hi[ok], angle = 90, code = 3,
                   length = 0.03, col = "steelblue4", lwd = 0.8)
  graphics::lines(t$arc, t$tau, type = "o", pch = 19, cex = 0.6, col = "darkorange3")
  leg <- c(if (band && "cb_lo" %in% names(t)) "uniform band", sprintf("%g%% robust CI", x$level), "tau(b)")
  col <- c(if (band && "cb_lo" %in% names(t)) grDevices::adjustcolor("steelblue", 0.3), "steelblue4", "darkorange3")
  lty <- c(if (band && "cb_lo" %in% names(t)) NA, 1, 1)
  pch <- c(if (band && "cb_lo" %in% names(t)) 15, NA, 19)
  if (!is.null(true_tau)) {
    graphics::lines(t$arc, true_tau, lty = 2, col = "black")
    leg <- c(leg, "truth")
    col <- c(col, "black")
    lty <- c(lty, 2)
    pch <- c(pch, NA)
  }
  graphics::legend("topright", legend = leg, col = col, lty = lty, pch = pch, bty = "n", cex = 0.8)
  invisible(x)
}
