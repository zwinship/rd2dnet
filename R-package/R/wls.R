# Local weighted least squares in network coordinates and its variance.

.VCE_TYPES <- c("hc0", "hc1", "hc2", "hc3", "cluster", "conley")

# Exponents (a, c) of s^a * o^c for all terms of total degree <= p, as a
# two-column matrix; lower orders are row prefixes.
.poly_exponents <- function(p) {
  a <- unlist(lapply(0:p, function(k) k:0))
  k <- unlist(lapply(0:p, function(k) rep(k, k + 1)))
  cbind(a = a, c = k - a)
}

.n_terms <- function(p) (p + 1) * (p + 2) / 2

# Columns s^a * o^c in .poly_exponents() order.
.poly_basis <- function(s, o, p) {
  e <- .poly_exponents(p)
  s <- as.numeric(s)
  o <- as.numeric(o)
  sp <- vapply(0:p, function(k) s^k, numeric(length(s)))
  op <- vapply(0:p, function(k) o^k, numeric(length(o)))
  sp <- matrix(sp, nrow = length(s))
  op <- matrix(op, nrow = length(o))
  sp[, e[, "a"] + 1, drop = FALSE] * op[, e[, "c"] + 1, drop = FALSE]
}

.exp_index <- function(p, a, c) {
  e <- .poly_exponents(p)
  which(e[, "a"] == a & e[, "c"] == c)
}

.wls <- function(X, y, w) {
  XtW <- t(X * w)
  Ginv <- .inv_checked(XtW %*% X)
  beta <- drop(Ginv %*% (XtW %*% y))
  list(beta = beta, Ginv = Ginv, resid = y - drop(X %*% beta), X = X, w = w)
}

# Weights a_i with c' beta = sum_i a_i y_i.
.loadings <- function(f, c) drop(f$X %*% (f$Ginv %*% c)) * f$w

.leverage <- function(f) f$w * rowSums((f$X %*% f$Ginv) * f$X)

# Per-observation influence psi_i with c' beta - truth ~ sum_i psi_i.
.influence <- function(f, c, vce = "hc1") {
  e <- f$resid
  n <- length(e)
  k <- length(f$beta)
  if (vce == "hc1") {
    e <- e * sqrt(n / max(n - k, 1))
  } else if (vce == "hc2") {
    e <- e / sqrt(pmax(1 - .leverage(f), 1e-12))
  } else if (vce == "hc3") {
    e <- e / pmax(1 - .leverage(f), 1e-12)
  }
  .loadings(f, c) * e
}

# Variance of sum(psi) under independence or clustering.
.variance <- function(psi, vce, cluster = NULL) {
  if (vce == "cluster") {
    sums <- rowsum(psi, cluster)
    g <- length(sums)
    return(sum(sums^2) * (if (g > 1) g / (g - 1) else 1))
  }
  sum(psi^2)
}
