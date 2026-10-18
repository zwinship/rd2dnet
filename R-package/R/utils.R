# Small internal helpers shared across the package.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Evaluate `code` with a temporary RNG seed, restoring the caller's RNG state.
.with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  set.seed(seed)
  code
}

# numpy.isclose(a, b) with its default tolerances.
.isclose <- function(a, b, rtol = 1e-5, atol = 1e-8) {
  out <- abs(a - b) <= atol + rtol * abs(b)
  out[is.na(out)] <- FALSE
  out
}

# numpy.interp: piecewise-linear interpolation on a non-decreasing grid that
# may contain ties (the value jumps across a tie), constant beyond the ends.
.interp <- function(x, grid, values) {
  n <- length(grid)
  j <- findInterval(x, grid)
  out <- numeric(length(x))
  lo <- j == 0
  hi <- j >= n
  out[lo] <- values[1]
  out[hi] <- values[n]
  mid <- !lo & !hi
  jm <- j[mid]
  t <- (x[mid] - grid[jm]) / (grid[jm + 1] - grid[jm])
  out[mid] <- values[jm] + t * (values[jm + 1] - values[jm])
  out
}

# Condition number in the 2-norm (numpy.linalg.cond).
.cond <- function(A) {
  if (!all(is.finite(A))) {
    return(Inf)
  }
  d <- svd(A, nu = 0, nv = 0)$d
  if (d[length(d)] == 0) Inf else d[1] / d[length(d)]
}

# Moore-Penrose pseudo-inverse with numpy's relative cutoff `rcond`.
.pinv <- function(A, rcond = 1e-15) {
  A <- as.matrix(A)
  s <- svd(A)
  keep <- s$d > rcond * max(s$d)
  if (!any(keep)) {
    return(matrix(0, ncol(A), nrow(A)))
  }
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# Minimum-norm least-squares solution (numpy.linalg.lstsq with rcond=None).
.lstsq <- function(A, b) {
  A <- as.matrix(A)
  .pinv(A, rcond = .Machine$double.eps * max(dim(A))) %*% b
}

# Inverse of a symmetric matrix after a conditioning check (singular local designs).
.inv_checked <- function(G, cond = 1e12) {
  if (.cond(G) > cond) {
    stop(.singular())
  }
  Gi <- solve(G)
  (Gi + t(Gi)) / 2
}

.singular <- function(msg = "local design matrix is singular") {
  structure(class = c("rd2dnet_singular", "error", "condition"),
            list(message = msg, call = NULL))
}

# k-th smallest value (1-based k clamped to the data).
.kth <- function(d, k) {
  k <- as.integer(min(max(k, 1), length(d)))
  sort(d, partial = k)[k]
}

.check_obs <- function(v, n, name) {
  v <- as.numeric(v)
  if (length(v) != n) {
    stop(sprintf("%s must have one value per observation (%d), got %d", name, n, length(v)), call. = FALSE)
  }
  v
}

.check_weights <- function(weights, n) {
  w <- .check_obs(weights, n, "weights")
  if (!all(is.finite(w)) || any(w < 0)) stop("weights must be finite and non-negative", call. = FALSE)
  w
}

.check_cluster <- function(cluster, n) {
  cl <- as.vector(cluster)
  if (is.factor(cluster)) cl <- as.character(cluster)
  if (length(cl) != n) {
    stop(sprintf("cluster must have one value per observation (%d), got %d", n, length(cl)), call. = FALSE)
  }
  if (anyNA(cl)) stop("cluster has missing values", call. = FALSE)
  cl
}

# Column sums by group: rowsum() that keeps a matrix.
.group_sum <- function(X, g) {
  rowsum(as.matrix(X), g, reorder = TRUE)
}

.check_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be a single number", name), call. = FALSE)
  }
  x
}
