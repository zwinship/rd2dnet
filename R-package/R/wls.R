# Local weighted least squares in network coordinates and its variance.

.VCE_TYPES <- c("hc0", "hc1", "hc2", "hc3", "cluster", "conley")

# Exponents (a, c) of s^a * o^c for all terms of total degree <= p, as a
# two-column matrix; lower orders are row prefixes. With dim = 1 the basis is a
# polynomial in s alone (no offset terms).
.poly_exponents <- function(p, dim = 2) {
  if (dim == 1) return(cbind(a = 0:p, c = rep(0L, p + 1)))
  a <- unlist(lapply(0:p, function(k) k:0))
  k <- unlist(lapply(0:p, function(k) rep(k, k + 1)))
  cbind(a = a, c = k - a)
}

.n_terms <- function(p, dim = 2) if (dim == 1) p + 1 else (p + 1) * (p + 2) / 2

# Columns s^a * o^c in .poly_exponents() order.
.poly_basis <- function(s, o, p, dim = 2) {
  e <- .poly_exponents(p, dim)
  s <- as.numeric(s)
  o <- as.numeric(o)
  sp <- vapply(0:p, function(k) s^k, numeric(length(s)))
  op <- vapply(0:p, function(k) o^k, numeric(length(o)))
  sp <- matrix(sp, nrow = length(s))
  op <- matrix(op, nrow = length(o))
  sp[, e[, "a"] + 1, drop = FALSE] * op[, e[, "c"] + 1, drop = FALSE]
}

.exp_index <- function(p, a, c, dim = 2) {
  e <- .poly_exponents(p, dim)
  which(e[, "a"] == a & e[, "c"] == c)
}

.wls <- function(X, y, w) {
  XtW <- t(X * w)
  G <- XtW %*% X
  # conditioning of the column-equilibrated Gram matrix: covariates in their own
  # units (income in dollars) must not make a well-posed fit look singular
  dg <- sqrt(pmax(diag(G), 1e-300))
  if (!all(is.finite(G)) || .cond(G / outer(dg, dg)) > 1e12) stop(.singular())
  Ginv <- solve(G / outer(dg, dg)) / outer(dg, dg)
  Ginv <- (Ginv + t(Ginv)) / 2
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
    e <- e / sqrt(pmax(1 - .leverage(f), 1e-8))
  } else if (vce == "hc3") {
    e <- e / pmax(1 - .leverage(f), 1e-8)
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

# Bias-reduced cluster adjustment (CR2) and Bell-McCaffrey degrees of freedom.
# Returns R with psi_i = R_i * resid_i (the variance is the sum over clusters
# of the squared cluster sums, with no small-sample factor) and the
# Satterthwaite degrees of freedom of that variance. R_g = D_g^(-1/2) a_g with
# a the loadings and D_g = ((I - H)(I - H)')_gg, H = X B, B = (X'WX)^(-1) X'W
# (Bell and McCaffrey 2002; Pustejovsky and Tipton 2018, working model the
# identity). D_g - I = U M U' with U = [X_g, W_g X_g], so D_g^(-1/2) comes from
# a small eigenproblem; directions with eigenvalue 0 (leverage one) get weight 0.
.cr2 <- function(f, cvec, cluster) {
  X <- f$X
  w <- f$w
  Gi <- f$Ginv
  k <- ncol(X)
  a <- .loadings(f, cvec)
  codes <- match(cluster, unique(cluster))
  G <- max(codes)
  Q <- Gi %*% crossprod(X * w) %*% Gi # B B'
  M <- rbind(cbind(Q, -Gi), cbind(-Gi, matrix(0, k, k)))
  R <- numeric(length(a))
  for (g in seq_len(G)) {
    m <- which(codes == g)
    U <- cbind(X[m, , drop = FALSE], X[m, , drop = FALSE] * w[m])
    qrU <- qr(U)
    Uo <- qr.Q(qrU)
    Rm <- qr.R(qrU)
    piv <- order(qrU$pivot) # undo column pivoting: U = Uo Rm[, piv]
    Rm <- Rm[, piv, drop = FALSE]
    e <- eigen(Rm %*% M %*% t(Rm), symmetric = TRUE)
    ev <- 1 + e$values # eigenvalues of D_g on span(U); 1 elsewhere
    fct <- ifelse(ev > 1e-10 * max(max(ev), 1), 1 / sqrt(pmax(ev, 1e-300)), 0) - 1
    am <- a[m]
    R[m] <- am + drop(Uo %*% (e$vectors %*% (fct * crossprod(e$vectors, crossprod(Uo, am)))))
  }
  # Bell-McCaffrey: P = Rb' (I - H)(I - H)' Rb with Rb the n x G block matrix of R
  Xr <- rowsum(X * R, codes, reorder = TRUE)
  Br <- rowsum(X * (R * w), codes, reorder = TRUE) %*% Gi
  P <- diag(as.numeric(rowsum(R * R, codes, reorder = TRUE)), G) - Xr %*% t(Br) - Br %*% t(Xr) +
    Xr %*% Q %*% t(Xr)
  den <- sum(P * P)
  list(R = R, df = if (den > 0) sum(diag(P))^2 / den else Inf)
}
