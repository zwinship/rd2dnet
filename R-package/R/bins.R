# Nearest-crossing bins: an alternative design to the default location-based one.
# Mirrors src/rd2dnet/bins.py (same steps, same tie rules), so the two agree to
# rounding.

#' Nearest-crossing bins
#'
#' An alternative to the default location-based design (which pools every
#' observation near an evaluation point and fits a polynomial in `(s, o)`).
#' Each observation belongs to the bin of the boundary crossing nearest to it
#' along the network, and [rdnet()] on the result fits each crossing
#' separately:
#'
#' 1. **Bins.** `R = s` is an observation's signed distance to its crossing;
#'    the local sample of bin `m` is its members, at kernel distance `|s|`.
#'    Where a street runs along the boundary, each of its nodes is a crossing.
#' 2. **Edge similarity.** Each street segment has a feature vector: the edge
#'    attributes in `features` plus its topology on the network as given
#'    (`n_adjacent`: edges sharing an end; `n_second`: edges two steps away;
#'    `degree_max`: the larger end degree; `clustering`: the mean clustering
#'    coefficient of its ends). `S` is the Mahalanobis distance between an
#'    observation's segment and the crossing's segment(s), with the covariance
#'    of the features over all segments.
#' 3. **Bin fits.** Weights are the product kernel `K1(R / h_R) * K2(S / h_S)`.
#'    The jump is common to the bin; slopes in `R` differ by side and by branch
#'    (the edge where the path from the crossing first passes a junction), and
#'    `R * S` enters with a slope per side.
#' 4. **Aggregation.** `wbate(res, "density")` or `wbate(res, "gls")`.
#'
#' Branches follow shortest paths over the metric's edge costs ignoring edge
#' directions (exact for undirected cost metrics; an approximation with directed
#' networks, transition penalties or precomputed metrics, where lengths are used
#' when edge costs are not available). [rdbwnet()] selects `h_R` per bin from the
#' members' fit in `R` with the `K2` weights, without the branch and `S` slopes.
#'
#' `h_S = "rot"` is [similarity_bandwidth()]. There is no data-driven selector
#' for `h_S`: `S` is not a running variable extrapolated to a cutoff, so the
#' MSE theory behind [rdbwnet()] does not cover it. Report
#' [similarity_sensitivity()].
#'
#' @param coords Output of [network_coordinates()] (its evaluation points are not used).
#' @param features Edge attributes describing street segments (capacity, speed,
#'   lanes, ...). Use attributes of the street, not of its piece: lengths are
#'   divided when edges are split, other attributes are copied. A segment's
#'   value is the mean over the stored edges between its ends; missing values
#'   are filled with the attribute's mean.
#' @param topology Add the topology features.
#' @param similarity Use the dissimilarity `S` (the kernel `K2` and the `R * S`
#'   slopes); `FALSE` gives plain bins.
#' @param branches Give each branch its own slopes on each side.
#' @param min_branch A branch gets its own slope only with at least this many
#'   distinct locations on its side of the local sample; smaller branches share
#'   the reference slope (the trunk's, or else the largest branch's).
#' @param h_S `"rot"`, a positive number, or `Inf` (no kernel in `S`).
#' @param kernel_S Kernel `K2` (default triangular).
#' @return Coordinates of class `rd2dnet_bins` (also `rd2dnet_coords`) with one
#'   evaluation point per bin, ordered along the boundary, and elements `bin`,
#'   `branch` (0 on the trunk), `similarity`, `kernel_weights`, `h_S`,
#'   `features` and `reference` (each bin's crossing features).
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' bins <- bin_coordinates(coords)
#' bins_table(bins)
#' res <- rdnet(city$obs$outcome, bins, h = 700, min_obs = 10)
#' wbate(res, "density")
#' @export
bin_coordinates <- function(coords, features = character(), topology = TRUE, similarity = TRUE,
                            branches = TRUE, min_branch = 10, h_S = "rot", kernel_S = "triangular") {
  .check_coords(coords)
  if (inherits(coords, "rd2dnet_bins")) {
    stop("coords are already bins; pass the output of network_coordinates()", call. = FALSE)
  }
  features <- as.character(features)
  net <- coords$network
  n <- nrow(net$coords)
  E <- edge_frame(net)
  node_on <- coords$node_on_boundary

  # one bin per crossing (boundary node) that some location reaches first
  first <- coords$first
  cross <- coords$lcross_node
  ok_loc <- !is.na(cross) & is.finite(coords$s[first])
  used <- sort(unique(cross[ok_loc]))
  arc <- if (length(used)) boundary_project(coords$boundary, net$coords[used, , drop = FALSE]) else numeric(0)
  o <- order(arc, method = "radix")
  rep_nodes <- used[o]
  arc <- arc[o]
  M <- length(rep_nodes)
  bin_node <- rep(NA_integer_, n)
  bin_node[rep_nodes] <- seq_len(M)
  bin_loc <- rep(NA_integer_, length(cross))
  bin_loc[ok_loc] <- bin_node[cross[ok_loc]]

  tree <- .crossing_tree(net, node_on, coords$node_treated, E)
  sn <- .at_nodes(coords$lsnap)
  Tl <- coords$treated[first]
  branch_loc <- if (branches) .location_branches(tree, sn, Tl) else integer(length(first))
  # a location whose tree path starts in another bin (ties, or a directed metric): trunk
  via_u <- sn$cost_u + tree$dist[sn$u]
  end <- ifelse(via_u <= sn$cost_v + tree$dist[sn$v], sn$u, sn$v)
  src <- tree$source[end]
  tb <- ifelse(src > 0, bin_node[pmax(src, 1L)], NA_integer_)
  other <- !is.na(bin_loc) & (is.na(tb) | tb != bin_loc)
  branch_loc[other] <- 0L

  out <- coords
  out$cache <- new.env(parent = emptyenv())
  out$eval_arc <- arc
  out$eval_node_arc <- arc
  out$eval_node <- rep_nodes
  out$eval_xy <- net$coords[rep_nodes, , drop = FALSE]
  out$eval_snap <- .snap(rep_nodes, rep_nodes, numeric(M), numeric(M), numeric(M))
  out$eval_masked <- boundary_is_masked(coords$boundary, arc)
  out$n_eval <- M
  out$bin <- bin_loc[coords$loc]
  out$branch <- branch_loc[coords$loc]
  out$bin_branches <- lapply(seq_len(M), function(m) {
    # in order of first appearance among the locations, as in Python and Stata
    unique(branch_loc[!is.na(bin_loc) & bin_loc == m & branch_loc > 0])
  })

  names_f <- character()
  if (similarity && (length(features) || topology)) {
    ef <- .edge_feature_matrix(net, E, features, topology)
    Z <- ef$Z
    names_f <- ef$names
    Omega <- if (nrow(Z) > 1) stats::cov(Z) else matrix(0, ncol(Z), ncol(Z))
    pinv <- .pinv_sym(Omega)
    # each bin's reference: mean features of the segments at its crossing
    ref <- matrix(0, M, ncol(Z))
    inc_bin <- c(bin_node[E$u], bin_node[E$v])
    inc_orig <- c(ef$orig, ef$orig)
    keep <- !is.na(inc_bin)
    pairs <- unique(cbind(inc_bin[keep], inc_orig[keep]))
    for (m in seq_len(M)) {
      ids <- pairs[pairs[, 1] == m, 2]
      if (length(ids)) ref[m, ] <- colMeans(Z[ids, , drop = FALSE])
    }
    row_loc <- .location_edge_rows(sn, tree)
    has_row <- !is.na(row_loc)
    dz <- matrix(0, length(first), ncol(Z))
    in_bin <- has_row & !is.na(bin_loc)
    dz[in_bin, ] <- Z[ef$orig[row_loc[in_bin]], , drop = FALSE] - ref[bin_loc[in_bin], , drop = FALSE]
    S_loc <- sqrt(pmax(rowSums((dz %*% pinv$inv) * dz), 0))
    S_loc[is.na(bin_loc)] <- NA_real_
    out$similarity <- S_loc[coords$loc]
    out$reference <- ref
    out$omega <- Omega
    out$feature_rank <- pinv$rank
  } else {
    out$similarity <- rep(NA_real_, coords$n_obs)
    out$reference <- matrix(0, M, 0)
    out$omega <- matrix(0, 0, 0)
    out$feature_rank <- 0L
  }
  colnames(out$reference) <- names_f
  out$features <- names_f
  out$use_similarity <- isTRUE(similarity)
  if (!is.numeric(min_branch) || length(min_branch) != 1 || min_branch < 1 || min_branch != round(min_branch)) {
    stop("min_branch must be a positive integer", call. = FALSE)
  }
  out$min_branch <- as.integer(min_branch)
  get_kernel(kernel_S)
  out$kernel_S <- kernel_S
  if (is.character(h_S)) {
    if (!identical(h_S, "rot")) stop("h_S must be 'rot', a positive number or Inf", call. = FALSE)
    h_S <- similarity_bandwidth(out$feature_rank)
  }
  class(out) <- c("rd2dnet_bins", class(coords))
  if (M == 0) warning("no observation reaches the boundary: there are no bins", call. = FALSE)
  with_similarity_bandwidth(out, h_S)
}

#' Rule-of-thumb similarity bandwidth
#'
#' `sqrt(2 * qchisq(level, k))`: if two segments' features were independent
#' draws with covariance `Omega`, `S^2 / 2` would follow a chi-squared
#' distribution with `k` (the rank of `Omega`) degrees of freedom. The default
#' keeps all but the 5% most dissimilar pairs under that reference. It is a
#' scale-free starting point, not an optimal choice.
#'
#' @param k Number of independent features.
#' @param level Quantile.
#' @return A bandwidth (`Inf` when `k = 0`).
#' @examples
#' similarity_bandwidth(4)
#' @export
similarity_bandwidth <- function(k, level = 0.95) {
  if (k <= 0) return(Inf)
  sqrt(2 * stats::qchisq(level, k))
}

#' Change the similarity bandwidth of bins
#'
#' @param bins Output of [bin_coordinates()].
#' @param h_S A positive number, or `Inf` for no similarity kernel.
#' @return The same bins with `h_S` and `kernel_weights` updated.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' bins <- with_similarity_bandwidth(bin_coordinates(coords), Inf)
#' is.null(bins$kernel_weights)
#' @export
with_similarity_bandwidth <- function(bins, h_S) {
  if (!inherits(bins, "rd2dnet_bins")) stop("expected bin_coordinates()", call. = FALSE)
  if (!is.numeric(h_S) || length(h_S) != 1 || is.na(h_S) || h_S <= 0) {
    stop("h_S must be positive (Inf for no similarity kernel)", call. = FALSE)
  }
  bins$h_S <- as.numeric(h_S)
  bins$cache <- new.env(parent = emptyenv())
  if (.has_similarity(bins) && is.finite(bins$h_S)) {
    k <- get_kernel(bins$kernel_S)
    S <- ifelse(is.finite(bins$similarity), bins$similarity, 0)
    bins$kernel_weights <- as.numeric(k(S / bins$h_S)) / as.numeric(k(0))
  } else {
    bins$kernel_weights <- NULL
  }
  bins
}

.has_similarity <- function(bins) isTRUE(bins$use_similarity) && length(bins$features) > 0

#' Table of nearest-crossing bins
#'
#' @param bins Output of [bin_coordinates()].
#' @return One row per bin: its position `arc`, crossing `node`, members per
#'   side `n0`, `n1`, number of branches, whether it is masked, and the
#'   crossing's features.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' bins_table(bin_coordinates(coords))
#' @export
bins_table <- function(bins) {
  if (!inherits(bins, "rd2dnet_bins")) stop("expected bin_coordinates()", call. = FALSE)
  M <- bins$n_eval
  b <- bins$bin
  T <- bins$treated
  out <- data.frame(bin = seq_len(M), arc = bins$eval_arc, node = bins$network$nodes$name[bins$eval_node],
                    n0 = tabulate(b[!is.na(b) & !T], M), n1 = tabulate(b[!is.na(b) & T], M),
                    n_branches = lengths(bins$bin_branches), masked = bins$eval_masked)
  for (k in seq_along(bins$features)) out[[bins$features[k]]] <- bins$reference[, k]
  out
}

#' @export
print.rd2dnet_bins <- function(x, ...) {
  cat(sprintf("<rd2dnet bins> n=%d, bins=%d, features=%d, h_S=%.4g\n", x$n_obs, x$n_eval,
              length(x$features), x$h_S))
  invisible(x)
}

# Extra slope regressors of bin j at `rows` (multiplied by s / h in the fit):
# one indicator per branch with its own slope (decided on `basis`, default the
# rows themselves), then S. NULL for other coordinates.
.local_slopes <- function(coords, j, rows, basis = rows) {
  if (!inherits(coords, "rd2dnet_bins")) return(NULL)
  cols <- lapply(.slope_branches(coords, j, basis), function(b) as.numeric(coords$branch[rows] == b))
  if (.has_similarity(coords)) cols[[length(cols) + 1]] <- coords$similarity[rows]
  if (!length(cols)) return(NULL)
  do.call(cbind, cols)
}

# Branches of bin j with their own slope in the local sample `rows`: on each side,
# those with at least min_branch distinct locations; the reference (no indicator)
# is the trunk when it qualifies, else the largest qualifying branch (the first
# to appear on ties). Smaller branches share the reference slope. Independent of
# column order, as in Python and Stata.
.slope_branches <- function(bins, j, rows) {
  order <- bins$bin_branches[[j]]
  br <- bins$branch[rows]
  Tr <- bins$treated[rows]
  loc <- bins$loc[rows]
  out <- integer()
  for (side in c(TRUE, FALSE)) {
    on <- Tr == side
    labs <- c(0L, order)
    cnt <- vapply(labs, function(x) length(unique(loc[on & br == x])), integer(1))
    ok <- order[cnt[-1] >= bins$min_branch]
    if (cnt[1] >= bins$min_branch || !length(ok)) {
      ref <- 0L
    } else {
      cok <- cnt[-1][cnt[-1] >= bins$min_branch]
      ref <- ok[which.max(cok)]
    }
    out <- c(out, setdiff(ok, ref))
  }
  order[order %in% out]
}

# Pseudo-inverse and rank of a symmetric PSD matrix (eigenvalues above rcond x the largest).
.pinv_sym <- function(A, rcond = 1e-10) {
  if (!length(A)) return(list(inv = A, rank = 0L))
  e <- eigen(A, symmetric = TRUE)
  top <- max(e$values, 0)
  keep <- if (top > 0) e$values > rcond * top else rep(FALSE, length(e$values))
  V <- e$vectors[, keep, drop = FALSE]
  list(inv = V %*% (t(V) / e$values[keep]), rank = sum(keep))
}

# Shortest-path tree from the boundary nodes and the branch of every node. A
# node's parent is, among its neighbours on a shortest path from the crossings,
# the one with the smallest index; candidates come before the node in
# (distance, index) order, so there are no cycles. A branch starts at the first
# junction on the path from the crossing (a node with three or more neighbours,
# or a crossing with two or more neighbours off the boundary on that side) and is
# labelled by the edge-frame row of its first edge; 0 is the trunk.
.crossing_tree <- function(net, node_on, node_T, E) {
  n <- nrow(net$coords)
  eu <- E$u
  ev <- E$v
  cost <- ifelse(is.finite(E$cost), E$cost, E$length)
  real <- eu != ev
  a <- c(eu[real], ev[real])
  b <- c(ev[real], eu[real])
  cc <- c(cost[real], cost[real])
  sources <- which(node_on)
  eng <- .new_engine(.csr(a, b, cc, n), n, seq_len(n) - 1L, seq_len(n) - 1L, TRUE)
  dist <- if (length(sources)) .engine_nearest(eng, sources)$dist else rep(Inf, n)
  tol <- 1e-9 * pmax(1, dist[b])
  da <- dist[a]
  db <- dist[b]
  cand <- is.finite(db) & !node_on[b] & (da + cc <= db + tol) &
    ((da < db - tol) | (abs(da - db) <= tol & a < b))
  pred <- integer(n)
  if (any(cand)) {
    m <- tapply(a[cand], b[cand], min)
    pred[as.integer(names(m))] <- as.integer(m)
  }
  lo <- pmin(eu, ev)[real]
  hi <- pmax(eu, ev)[real]
  key <- unique((lo - 1) * n + (hi - 1))
  lo <- key %/% n + 1
  hi <- key %% n + 1
  deg <- tabulate(c(lo, hi), n)
  aa <- c(lo, hi)
  bb <- c(hi, lo)
  off <- node_on[aa] & !node_on[bb]
  side_T <- tabulate(aa[off & node_T[bb]], n)
  side_C <- tabulate(aa[off & !node_T[bb]], n)
  full_key <- (pmin(eu, ev) - 1) * n + (pmax(eu, ev) - 1)
  row_of <- function(p, x) match((pmin(p, x) - 1) * n + (pmax(p, x) - 1), full_key)
  splits <- function(p, ts) ifelse(node_on[p], ifelse(ts, side_T[p], side_C[p]) >= 2, deg[p] >= 3)
  label <- integer(n)
  child <- which(pred > 0)
  st <- splits(pred[child], node_T[child])
  label[child[st]] <- row_of(pred[child[st]], child[st])
  root <- ifelse(pred > 0, pred, seq_len(n))
  anc <- pred
  repeat {
    m <- which(anc > 0)
    if (!length(m)) break
    up <- anc[m]
    label[m] <- ifelse(label[up] > 0, label[up], label[m])
    root[m] <- root[up]
    anc[m] <- anc[up]
  }
  root <- ifelse(is.finite(dist) & node_on[root], root, 0L)
  list(dist = dist, pred = pred, label = label, splits = splits, row_of = row_of, source = root)
}

# Locations that sit exactly on an end node of their segment, placed on that
# node: such a location lies on every segment that meets the node, and which one
# the snapper picked is a tie (broken differently in Python, R and Stata).
.at_nodes <- function(sn) {
  tiny <- 1e-9 * pmax(1, sn$cost_u + sn$cost_v)
  at_u <- sn$u != sn$v & sn$cost_u <= tiny
  at_v <- sn$u != sn$v & !at_u & sn$cost_v <= tiny
  at <- at_u | at_v
  node <- ifelse(at_u, sn$u, sn$v)
  .snap(ifelse(at, node, sn$u), ifelse(at, node, sn$v), ifelse(at, 0, sn$cost_u), ifelse(at, 0, sn$cost_v), sn$depth)
}

# Branch of each location (one per element of the location snap `sn`).
.location_branches <- function(tree, sn, Tl) {
  swap <- sn$cost_u + tree$dist[sn$u] > sn$cost_v + tree$dist[sn$v]
  a <- ifelse(swap, sn$v, sn$u) # the end the location is reached through
  b <- ifelse(swap, sn$u, sn$v)
  out <- tree$label[a]
  new <- out == 0 & a != b & tree$splits(a, Tl)
  out[new] <- tree$row_of(a[new], b[new])
  out
}

# Edge-frame row of the segment each location sits on (NA at a crossing node).
.location_edge_rows <- function(sn, tree) {
  out <- rep(NA_integer_, length(sn$u))
  on_edge <- sn$u != sn$v
  out[on_edge] <- tree$row_of(sn$u[on_edge], sn$v[on_edge])
  at_node <- !on_edge & tree$pred[sn$u] > 0 # on a node: the segment it was reached along
  out[at_node] <- tree$row_of(tree$pred[sn$u[at_node]], sn$u[at_node])
  out
}

# Original-segment id of each edge-frame row (nodes the package inserted, at
# boundary crossings and at locations, with two neighbours are skipped) and the
# features of each original segment.
.edge_feature_matrix <- function(net, E, attrs, topology) {
  n <- nrow(net$coords)
  eu <- E$u
  ev <- E$v
  real <- eu != ev
  deg_split <- tabulate(c(eu[real], ev[real]), n)
  labels <- as.character(net$nodes$name)
  inserted <- grepl("^rdnet_(crossing|point)_", labels)
  transparent <- inserted & deg_split == 2
  x <- c(eu, ev)
  y <- c(ev, eu)
  tr <- transparent[x]
  nb <- matrix(NA_integer_, n, 2)
  if (any(tr)) {
    xs <- x[tr]
    ys <- y[tr]
    pos <- stats::ave(seq_along(xs), xs, FUN = seq_along)
    keep <- pos <= 2
    nb[cbind(xs[keep], pos[keep])] <- ys[keep]
  }
  walk <- function(a, prev) {
    for (i in seq_len(n)) {
      if (!transparent[a]) return(a)
      nxt <- if (!is.na(nb[a, 1]) && nb[a, 1] == prev) nb[a, 2] else nb[a, 1]
      prev <- a
      a <- nxt
    }
    a
  }
  oa <- eu
  ob <- ev
  for (r in which(transparent[eu] | transparent[ev])) {
    oa[r] <- walk(eu[r], ev[r])
    ob[r] <- walk(ev[r], eu[r])
  }
  okey <- (pmin(oa, ob) - 1) * n + (pmax(oa, ob) - 1)
  ukeys <- sort(unique(okey))
  orig <- match(okey, ukeys)
  m <- length(ukeys)
  A <- ukeys %/% n + 1
  B <- ukeys %% n + 1
  cols <- list()
  St <- net$edges
  skey <- (pmin(St$from, St$to) - 1) * n + (pmax(St$from, St$to) - 1)
  rkey <- (pmin(eu, ev) - 1) * n + (pmax(eu, ev) - 1)
  for (attr in attrs) {
    if (!attr %in% names(St)) stop(sprintf("edge attribute '%s' is missing on every segment", attr), call. = FALSE)
    sv <- suppressWarnings(.first_value(St[[attr]]))
    ok <- is.finite(sv)
    # a segment's value: the mean over the stored edges between its ends
    rs <- tapply(sv[ok], skey[ok], mean)
    vals <- unname(rs[match(rkey, as.numeric(names(rs)))])
    vals <- as.numeric(vals)
    fin <- is.finite(vals)
    sums <- as.numeric(tapply(vals[fin], factor(orig[fin], levels = seq_len(m)), sum))
    sums[is.na(sums)] <- 0
    cnt <- tabulate(orig[fin], m)
    col <- ifelse(cnt > 0, sums / cnt, NA_real_)
    miss <- !is.finite(col)
    if (all(miss)) stop(sprintf("edge attribute '%s' is missing on every segment", attr), call. = FALSE)
    if (any(miss)) {
      warning(sprintf("edge attribute '%s' is missing on %d of %d segments; filled with its mean", attr,
                      sum(miss), m), call. = FALSE)
      col[miss] <- mean(col[!miss])
    }
    cols[[attr]] <- col
  }
  if (topology) {
    rl <- A != B
    deg <- tabulate(c(A[rl], B[rl]), n)
    cols$n_adjacent <- ifelse(rl, deg[A] + deg[B] - 2, deg[A])
    cols$n_second <- .second_edges(A, B, rl, n)
    cols$degree_max <- pmax(deg[A], deg[B])
    cl <- .clustering(A[rl], B[rl], n)
    cols$clustering <- (cl[A] + cl[B]) / 2
  }
  Z <- if (length(cols)) do.call(cbind, unname(cols)) else matrix(0, m, 0)
  list(orig = orig, Z = Z, names = names(cols))
}

# Number of segments two steps away from each segment: those sharing an end with
# an adjacent segment, other than the segment itself and its adjacent segments.
.second_edges <- function(A, B, rl, n) {
  m <- length(A)
  e <- which(rl)
  inc_node <- c(A[e], B[e])
  inc_edge <- c(e, e)
  o <- order(inc_node, inc_edge)
  inc_node <- inc_node[o]
  inc_edge <- inc_edge[o]
  # adjacent pairs (e, f) sharing a node
  pe <- integer(0)
  pf <- integer(0)
  k <- 1L
  repeat {
    i <- seq_len(length(inc_node) - k)
    if (!length(i)) break
    same <- inc_node[i] == inc_node[i + k]
    if (!any(same)) break
    pe <- c(pe, inc_edge[i][same], inc_edge[i + k][same])
    pf <- c(pf, inc_edge[i + k][same], inc_edge[i][same])
    k <- k + 1L
  }
  ak <- unique((pe - 1) * m + (pf - 1))
  pe <- ak %/% m + 1
  pf <- ak %% m + 1
  o <- order(pe, pf)
  pe <- pe[o]
  pf <- pf[o]
  start <- match(seq_len(m), pe)
  cnt <- tabulate(pe, m)
  # two steps: (e, f) then (f, g)
  rep_e <- rep(pe, cnt[pf])
  g <- pf[rep(start[pf], cnt[pf]) + sequence(cnt[pf]) - 1L]
  key <- (rep_e - 1) * m + (g - 1)
  key <- unique(key[g != rep_e])
  key <- key[!key %in% ak]
  tabulate(key %/% m + 1, m)
}

# Local clustering coefficient of every node of a simple undirected graph (0 below degree 2).
.clustering <- function(a, b, n) {
  lo <- pmin(a, b)
  hi <- pmax(a, b)
  ekey <- unique((lo - 1) * n + (hi - 1))
  lo <- ekey %/% n + 1
  hi <- ekey %% n + 1
  x <- c(lo, hi)
  y <- c(hi, lo)
  o <- order(x, y)
  x <- x[o]
  y <- y[o]
  deg <- tabulate(x, n)
  tri <- numeric(n)
  k <- 1L
  repeat {
    i <- seq_len(length(x) - k)
    if (!length(i)) break
    same <- x[i] == x[i + k]
    if (!any(same)) break
    y1 <- y[i][same]
    y2 <- y[i + k][same]
    linked <- ((pmin(y1, y2) - 1) * n + (pmax(y1, y2) - 1)) %in% ekey
    tri <- tri + tabulate(x[i][same][linked], n)
    k <- k + 1L
  }
  den <- deg * (deg - 1) / 2
  ifelse(den > 0, tri / den, 0)
}

#' Pooled nearest-crossing BATE with crossing fixed effects
#'
#' For boundaries with many sparsely populated crossings: one weighted
#' regression over the members of all bins,
#' `y = a_m + T_B T + poly_p(R) + T poly_p(R) + R Z_m' th1 + T R Z_m' th2 + d1 R S + d2 T R S + covs`,
#' with `T` the 0/1 treatment indicator, `R = s / h` (each side's bandwidth),
#' `a_m` one intercept per bin (the jump is identified within bins only),
#' `Z_m` the bin's crossing features (standardised across bins) and `S` the edge
#' dissimilarity. Weights are the product kernel times, with `ipw = TRUE`,
#' `1 / p_m` for treated and `1 / (1 - p_m)` for control members, where `p_m` is
#' the kernel-weighted treated share of bin `m`. Plain fixed-effects OLS weights
#' a bin by `p_m (1 - p_m)` times its mass; after reweighting both sides of every
#' bin carry the bin's kernel mass, so `T_B` weights bins by their local sample
#' (the density-weighted BATE). `tau` comes from the order-`p` fit and the
#' bias-corrected `tau_bc` (with `se_rb`, the interval and the p-value) from
#' the order-`q` fit, as in [rdnet()].
#'
#' @inheritParams rdnet
#' @param bins Output of [bin_coordinates()].
#' @param h Bandwidth: a scalar, `c(h_control, h_treated)`, one value per bin, a
#'   `M x 2` matrix, or `NULL` to select per bin with [rdbwnet()] (bins whose
#'   selection fails use the median of the selected ones).
#' @param ipw Reweight bins to the density-weighted BATE (default).
#' @param bin_slopes Interact the slopes with the bins' crossing features.
#' @param min_obs Minimum distinct locations per side for a bin to enter.
#' @param bwselect Bandwidth selector used when `h` is `NULL` (see [rdbwnet()]).
#' @param scaleregul Bandwidth-selection regularisation passed to [rdbwnet()] when `h` is selected.
#' @return An object of class `rdnet_fe`: `tau`, `se`, `tau_bc`, `se_rb`,
#'   `ci_lo`, `ci_hi`, `pvalue`, `df`, the per-bin table `bins` (bandwidths,
#'   locations per side, `phat`, `used`) and the slope `terms`.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' rdnet_fe(city$obs$outcome, bin_coordinates(coords), h = 700)
#' @export
rdnet_fe <- function(y, bins, h = NULL, p = 1, q = NULL, kernel = "triangular", vce = "hc3", cluster = NULL,
                     cluster_vce = "cr2", conley_cutoff = NULL, conley_distance = "network", weights = NULL,
                     covs = NULL, ipw = TRUE, bin_slopes = TRUE, level = 95, min_obs = 1, bwselect = "mserd",
                     scaleregul = 3) {
  if (!inherits(bins, "rd2dnet_bins")) stop("rdnet_fe() needs bin_coordinates()", call. = FALSE)
  y <- .check_obs(y, bins$n_obs, "y")
  p <- .as_int(p, "p", 0L)
  q <- if (is.null(q)) p + 1L else .as_int(q, "q")
  if (q < p) stop("need 0 <= p <= q", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || !(level > 1 && level < 100)) {
    stop(sprintf("level is a percentage, e.g. 95 (got %s)", format(level)), call. = FALSE)
  }
  kfun <- get_kernel(kernel)
  intensity <- if (is.null(weights)) rep(1, length(y)) else .check_weights(weights, length(y))
  Zc <- NULL
  if (!is.null(covs)) {
    Zc <- as.matrix(covs)
    storage.mode(Zc) <- "double"
    if (nrow(Zc) != length(y)) stop("covs must have one row per observation", call. = FALSE)
  }
  if (!is.null(cluster)) {
    vce <- "cluster"
    cluster <- .check_cluster(cluster, length(y))
  }
  if (!is.null(conley_cutoff)) {
    vce <- "conley"
    if (conley_cutoff <= 0) stop("conley_cutoff must be positive", call. = FALSE)
  }
  if (!vce %in% .VCE_TYPES) stop("vce must be 'hc0'-'hc3', 'cluster' or 'conley'", call. = FALSE)
  if (vce == "cluster" && is.null(cluster)) stop("vce = 'cluster' needs a cluster variable", call. = FALSE)
  if (!cluster_vce %in% c("cr2", "cr1")) stop("cluster_vce must be 'cr2' or 'cr1'", call. = FALSE)
  M <- bins$n_eval
  if (is.null(h)) {
    bw <- rdbwnet(y, bins, p = p, kernel = kernel, bwselect = bwselect, weights = weights,
                  scaleregul = scaleregul, vce = if (vce %in% c("hc0", "hc1", "hc2", "hc3")) vce else "hc0")
    H <- as.matrix(bw[, c("h0", "h1")])
    good <- apply(is.finite(H) & H > 0, 1, all)
    if (!any(good)) stop("no bin's bandwidth could be selected; pass h", call. = FALSE)
    med <- apply(H[good, , drop = FALSE], 2, stats::median)
    H[!good, ] <- matrix(med, sum(!good), 2, byrow = TRUE)
  } else {
    H <- .bandwidth_matrix(h, M)
  }
  T <- bins$treated
  s <- bins$s
  b <- bins$bin
  bb <- ifelse(is.na(b), 1L, b)
  kw <- if (is.null(bins$kernel_weights)) rep(1, length(y)) else bins$kernel_weights
  ok <- !is.na(b) & is.finite(s) & is.finite(y) & intensity > 0 & !bins$eval_masked[bb]
  if (!is.null(Zc)) ok <- ok & stats::complete.cases(Zc)
  hside <- ifelse(T, H[bb, 2], H[bb, 1])
  W <- ifelse(ok, kfun(abs(s) / hside) * kw * intensity, 0)
  ok <- ok & W > 0
  nl <- function(sel) vapply(seq_len(M), function(m) length(unique(bins$loc[sel & b == m])), integer(1))
  n1 <- nl(ok & T)
  n0 <- nl(ok & !T)
  w1 <- as.numeric(tapply(c(W[ok & T], 0), factor(c(b[ok & T], NA), levels = seq_len(M)), sum))
  w0 <- as.numeric(tapply(c(W[ok & !T], 0), factor(c(b[ok & !T], NA), levels = seq_len(M)), sum))
  w1[is.na(w1)] <- 0
  w0[is.na(w0)] <- 0
  phat <- ifelse(w1 + w0 > 0, w1 / (w1 + w0), NA_real_)
  use_bin <- n0 >= min_obs & n1 >= min_obs
  tab <- data.frame(bin = seq_len(M), arc = bins$eval_arc, h0 = H[, 1], h1 = H[, 2], n0 = n0, n1 = n1,
                    phat = phat, used = use_bin)
  idx <- which(ok & use_bin[bb])
  if (!any(use_bin) || !length(idx)) stop("no bin has observations on both sides within the bandwidth", call. = FALSE)
  bi <- b[idx]
  Ti <- as.numeric(T[idx])
  w <- W[idx]
  if (ipw) w <- w * ifelse(T[idx], 1 / phat[bi], 1 / (1 - phat[bi]))
  R <- s[idx] / hside[idx]
  ub <- sort(unique(bi))
  binv <- match(bi, ub)
  FE <- matrix(0, length(idx), length(ub))
  FE[cbind(seq_along(idx), binv)] <- 1
  extra <- list()
  terms <- character()
  if (bin_slopes && ncol(bins$reference)) {
    Zm <- bins$reference[ub, , drop = FALSE]
    sdv <- apply(Zm, 2, function(v) sqrt(mean((v - mean(v))^2)))
    vary <- sdv > 1e-12 * pmax(apply(abs(Zm), 2, max), 1)
    if (any(vary)) {
      Zs <- sweep(sweep(Zm[, vary, drop = FALSE], 2, colMeans(Zm[, vary, drop = FALSE])), 2, sdv[vary], "/")
      Zi <- Zs[binv, , drop = FALSE]
      extra <- c(extra, list(R * Zi, (Ti * R) * Zi))
      nm <- bins$features[vary]
      terms <- c(terms, paste0("R*", nm), paste0("T*R*", nm))
    }
  }
  if (.has_similarity(bins)) {
    Si <- bins$similarity[idx]
    extra <- c(extra, list(R * Si, Ti * R * Si))
    terms <- c(terms, "R*S", "T*R*S")
  }
  z <- stats::qnorm(0.5 + level / 200)
  cr2 <- vce == "cluster" && cluster_vce == "cr2"
  base_vce <- if (vce %in% c("hc1", "hc2", "hc3")) vce else "hc0"
  fit <- function(order) {
    P <- if (order > 0) .poly_basis(R, numeric(length(R)), order, dim = 1)[, -1, drop = FALSE] else NULL
    X <- do.call(cbind, c(list(FE, Ti), if (!is.null(P)) list(P, Ti * P), extra, if (!is.null(Zc)) list(Zc[idx, , drop = FALSE])))
    X <- unname(X)
    keep <- .identifiable(X, w)
    tcol <- length(ub) + 1L
    if (!keep[tcol]) stop("the treatment effect is not identified (no bin has both sides)", call. = FALSE)
    cvec <- rep(0, sum(keep))
    cvec[sum(keep[seq_len(tcol)])] <- 1
    X <- X[, keep, drop = FALSE]
    f <- .wls(X, y[idx], w)
    tau <- sum(cvec * f$beta)
    df <- Inf
    if (vce == "conley") {
      psi <- .loadings(f, cvec) * f$resid
      v <- drop(.hac_quadratic(bins, idx, psi, conley_cutoff, conley_distance))
    } else if (cr2) {
      adj <- .cr2(f, cvec, cluster[idx])
      df <- adj$df
      v <- sum(tapply(adj$R * f$resid, cluster[idx], sum)^2)
    } else if (vce == "cluster") {
      v <- .variance(.influence(f, cvec, "hc0"), "cluster", cluster[idx]) *
        .cluster_dof(length(idx), ncol(X), cluster[idx])
    } else {
      v <- sum(.influence(f, cvec, base_vce)^2)
    }
    list(tau = tau, se = sqrt(max(v, 0)), df = df)
  }
  fp <- fit(p)
  fq <- fit(q)
  df <- fq$df
  fin <- is.finite(df) && df > 0
  zq <- if (fin) stats::qt(stats::pnorm(z), df) else z
  pv <- if (fq$se > 0) {
    if (fin) 2 * stats::pt(-abs(fq$tau / fq$se), df) else 2 * stats::pnorm(-abs(fq$tau / fq$se))
  } else {
    NA_real_
  }
  structure(list(tau = fp$tau, se = fp$se, tau_bc = fq$tau, se_rb = fq$se, ci_lo = fq$tau - zq * fq$se,
                 ci_hi = fq$tau + zq * fq$se, pvalue = pv, level = level, p = p, q = q, vce = vce, df = df,
                 bins = tab, terms = terms),
            class = "rdnet_fe")
}

#' @export
print.rdnet_fe <- function(x, digits = 3, ...) {
  cat("rd2dnet: pooled nearest-crossing BATE with bin fixed effects\n")
  u <- x$bins$used
  cat(sprintf("  p=%d  q=%d  vce=%s  bins used: %d of %d  n0=%d  n1=%d\n", x$p, x$q, x$vce, sum(u),
              nrow(x$bins), sum(x$bins$n0[u]), sum(x$bins$n1[u])))
  cat(sprintf("  slope terms: %s\n\n", if (length(x$terms)) paste(x$terms, collapse = ", ") else "none"))
  print(round(unlist(x[c("tau", "se", "tau_bc", "se_rb", "ci_lo", "ci_hi", "pvalue")]), digits))
  invisible(x)
}

#' BATE over a range of similarity bandwidths
#'
#' Refits [rdnet()] on the bins for each `h_S` (with `h` not given, `h_R` is
#' selected again each time) and aggregates with [wbate()].
#'
#' @param y Outcome.
#' @param bins Output of [bin_coordinates()].
#' @param h_S Similarity bandwidths; default half, one and two times the rule
#'   of thumb, and `Inf` (no similarity kernel).
#' @param wbate Aggregation rule or weights passed to [wbate()].
#' @param ... Passed to [rdnet()].
#' @return A data frame with one row per `h_S`: the number of bins estimated and
#'   the aggregate effect.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' similarity_sensitivity(city$obs$outcome, bin_coordinates(coords), h = 700, min_obs = 10)
#' @export
similarity_sensitivity <- function(y, bins, ..., h_S = NULL, wbate = "density") {
  if (!inherits(bins, "rd2dnet_bins")) stop("similarity_sensitivity() needs bin_coordinates()", call. = FALSE)
  if (is.null(h_S)) {
    base <- similarity_bandwidth(bins$feature_rank)
    h_S <- if (is.finite(base)) c(0.5, 1, 2) * base else numeric(0)
    h_S <- c(h_S, Inf)
  }
  agg <- wbate
  aggregate <- get("wbate", envir = asNamespace("rd2dnet"), mode = "function")
  rows <- lapply(h_S, function(hs) {
    res <- rdnet(y, with_similarity_bandwidth(bins, hs), ...)
    a <- tryCatch(aggregate(res, agg), error = function(e) {
      c(tau = NA_real_, tau_bc = NA_real_, se_rb = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, pvalue = NA_real_)
    })
    data.frame(h_S = hs, bins_ok = sum(res$table$status == "ok"), t(a))
  })
  do.call(rbind, rows)
}
