# Manipulation (sorting) test at the boundary, after Bugni, Crippa and Restrepo
# (2026). Mirrors the Python package's manipulation.py step for step, so the
# clusters, the relabelling draws and the results are identical.

.MT_M <- 2147483647 # 2^31 - 1
.MT_A <- 16807
.MT_HALF <- 1073741823 # draws 1..2^30-1 are "treated": exactly half of 1..2^31-2

# reps x q Bernoulli(1/2) labels from the portable Park-Miller stream.
.lcg_labels <- function(seed, reps, q) {
  x <- (as.numeric(seed) %% (.MT_M - 1)) + 1
  out <- matrix(0L, reps, q)
  for (r in seq_len(reps)) {
    for (j in seq_len(q)) {
      x <- (x * .MT_A) %% .MT_M
      out[r, j] <- as.integer(x <= .MT_HALF)
    }
  }
  out
}

# Optimal 1-D k-means partitions of sorted positions for k = 1..k_max; for each
# k, the (1-based) index into pos of the first element of every cluster.
.kmeans_1d_partitions <- function(pos, k_max) {
  vals <- unique(pos)
  start <- match(vals, pos)
  w <- as.numeric(tabulate(match(pos, vals), length(vals)))
  m <- length(vals)
  k_max <- min(k_max, m)
  s1 <- c(0, cumsum(w * vals))
  s2 <- c(0, cumsum(w * vals * vals))
  sw <- c(0, cumsum(w))
  cost <- function(i, j) { # distinct values i..j (1-based, inclusive); i may be a vector
    n <- sw[j + 1] - sw[i]
    a <- s1[j + 1] - s1[i]
    pmax((s2[j + 1] - s2[i]) - a * a / n, 0)
  }
  D <- matrix(Inf, k_max, m)
  arg <- matrix(1L, k_max, m)
  D[1, ] <- cost(rep(1, m), seq_len(m))
  if (k_max >= 2) {
    for (k in 2:k_max) {
      for (j in k:m) {
        i <- k:j # the last cluster is values i..j
        cc <- D[k - 1, i - 1] + cost(i, j)
        best <- min(cc)
        bi <- i[which(cc <= best + 1e-12 * max(1, abs(best)))[1]]
        D[k, j] <- cc[bi - k + 1]
        arg[k, j] <- bi
      }
    }
  }
  lapply(seq_len(k_max), function(k) {
    firsts <- integer(0)
    j <- m
    kk <- k
    while (kk >= 1) {
      i <- if (kk > 1) arg[kk, j] else 1L
      firsts <- c(i, firsts)
      j <- i - 1
      kk <- kk - 1
    }
    start[firsts]
  })
}

# S for each row of labels (draws x q, in position order).
.mt_statistic <- function(labels, firsts, q) {
  P <- cbind(0, t(apply(labels, 1, cumsum)))
  S <- numeric(nrow(labels))
  for (f in firsts) {
    ends <- c(f[-1], q + 1)
    n <- ends - f
    tt <- P[, ends, drop = FALSE] - P[, f, drop = FALSE]
    dev <- sweep(abs(sweep(tt, 2, n, "/") - 0.5), 2, sqrt(n), "*")
    S <- pmax(S, apply(dev, 1, max))
  }
  S
}

# Warn when most locations lie on the network itself (no sorting test is valid then).
.warn_if_on_network <- function(coords) {
  depth <- coords$snap$depth[!duplicated(coords$loc)]
  span <- apply(coords$network$coords, 2, function(z) diff(range(z)))
  extent <- sqrt(sum(span^2))
  if (length(depth) && mean(depth <= 1e-6 * max(extent, 1)) > 0.5) {
    warning(paste0("most locations lie on the network itself (stations on a river, points on streets); ",
                   "sorting tests assume units spread over the plane, and in Monte Carlo runs none of them ",
                   "held its size for points on streets near a curved boundary"), call. = FALSE)
  }
}

#' Manipulation (sorting) test at the boundary
#'
#' The test of Bugni, Crippa and Restrepo (2026): without sorting across the
#' boundary, the locations closest to it are about equally likely to lie on
#' either side, anywhere along the boundary. The `q` locations nearest the
#' boundary (straight-line distance in the plane: sorting is about where units
#' locate, not about the network) are projected onto it, grouped into `k`
#' contiguous clusters for every `k` from 1 to `k_max`, and the statistic is the
#' largest standardised imbalance `max sqrt(n) |share treated - 1/2|`. Its null
#' distribution comes from relabelling the `q` locations as independent
#' Bernoulli(1/2) draws with the clusters held fixed.
#'
#' Clusters are the exact optimal 1-D k-means partitions of the positions along
#' the boundary (a closed ring is cut at its largest gap), and the draws come
#' from the Park-Miller generator, so the Python, R and Stata versions give the
#' same answer. Locations count once; locations on the boundary line, in masked
#' stretches, or projecting onto an end of an open boundary are left out.
#'
#' The test has power against sorting that offsets along the boundary, which a
#' density test of the distance to the boundary cannot detect; against sorting in
#' one direction everywhere, [density()] is more powerful. It is not valid for
#' observations that lie on the network itself (stations along a river, points
#' on streets).
#'
#' @param coords Output of [network_coordinates()].
#' @param q Number of locations nearest the boundary, or a vector of values (the
#'   result's `path` then has one row per value). Defaults to 200, or a tenth of
#'   the usable locations when there are fewer than 2,000.
#' @param k_max Largest number of clusters along the boundary (default `q`).
#' @param reps Relabelling draws for the null distribution.
#' @param level Confidence level (percent) for the reported critical value.
#' @param seed Seed of the relabelling draws.
#' @return A list of class `rd2dnet_manipulation`: `statistic`, `pvalue`, `crit`,
#'   `q`, `n_treated`, `distance` (of the farthest selected location), `k_max`,
#'   `reps`, `level` and `path` (one row per `q`).
#' @references Bugni, F. A., Crippa, F. and Restrepo, D. (2026). Manipulation
#'   Testing in Boundary Discontinuity Designs. arXiv:2609.00350.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' manipulation_test(coords, q = c(50, 200), reps = 500)
#' @export
manipulation_test <- function(coords, q = NULL, k_max = NULL, reps = 10000, level = 95, seed = 0) {
  .check_coords(coords)
  .warn_if_on_network(coords)
  if (!(is.numeric(level) && length(level) == 1 && level > 1 && level < 100)) {
    stop(sprintf("level is a percentage, e.g. 95 (got %s)", format(level)), call. = FALSE)
  }
  if (reps < 100) stop("reps must be at least 100", call. = FALSE)
  b <- coords$boundary
  first <- which(!duplicated(coords$loc))
  xy <- coords$xy[first, , drop = FALSE]
  T <- as.logical(coords$treated[first])
  dist <- boundary_distance(b, xy)
  arc <- boundary_project(b, xy)
  L <- b$length
  ok <- is.finite(dist) & dist > b$tol & !boundary_is_masked(b, arc)
  if (!isTRUE(b$closed)) ok <- ok & arc > 1e-9 * L & arc < L * (1 - 1e-9)
  idx <- which(ok)
  n_use <- length(idx)
  qs <- if (is.null(q)) (if (n_use >= 2000) 200 else max(n_use %/% 10, 1)) else q
  rows <- list()
  for (qq in qs) {
    qq <- as.integer(qq)
    if (qq < 2 || qq > n_use) {
      stop(sprintf("q must be between 2 and the number of usable locations (%d), got %d", n_use, qq), call. = FALSE)
    }
    if (qq < 20) warning(sprintf("q = %d locations give a test with little power", qq), call. = FALSE)
    if (qq > 0.4 * n_use) {
      warning(sprintf(paste0("q = %d is a large share of the %d usable locations; the test's approximation ",
                             "assumes the selected locations are close to the boundary"), qq, n_use), call. = FALSE)
    }
    sel <- idx[order(dist[idx], idx)][seq_len(qq)] # nearest first; ties by location order
    pos <- arc[sel]
    if (isTRUE(b$closed)) { # cut the ring at the largest gap between neighbouring positions
      srt <- sort(unique(pos))
      gaps <- diff(c(srt, srt[1] + L))
      g <- which.max(gaps)
      cut <- srt[(g %% length(srt)) + 1]
      pos <- (pos - cut) %% L
    }
    po <- order(pos, sel) # position order; ties by location order
    pos <- pos[po]
    lab <- as.integer(T[sel][po])
    km <- if (is.null(k_max)) qq else min(as.integer(k_max), qq)
    firsts <- .kmeans_1d_partitions(pos, km)
    S <- .mt_statistic(matrix(lab, 1), firsts, qq)
    null <- .mt_statistic(.lcg_labels(seed, reps, qq), firsts, qq)
    rows[[length(rows) + 1]] <- data.frame(
      q = qq, statistic = S, pvalue = mean(null >= S - 1e-12 * max(1, S)),
      crit = sort(null)[ceiling(level / 100 * length(null))], n_treated = sum(lab),
      distance = max(dist[sel]), k_max = length(firsts))
  }
  path <- do.call(rbind, rows)
  out <- c(as.list(path[1, ]), list(reps = as.integer(reps), level = level, path = path))
  class(out) <- "rd2dnet_manipulation"
  out
}

#' @export
print.rd2dnet_manipulation <- function(x, ...) {
  cat(sprintf(paste0("rd2dnet manipulation test (Bugni, Crippa and Restrepo): q=%d locations within %.4g of ",
                     "the boundary, %d treated; S=%.3f, p=%.4f (critical value %.3f at %g%%, %d draws)\n"),
              x$q, x$distance, x$n_treated, x$statistic, x$pvalue, x$crit, x$level, x$reps))
  if (nrow(x$path) > 1) print(x$path, row.names = FALSE)
  invisible(x)
}
