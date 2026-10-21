# Network-native coordinates relative to boundary evaluation points.
#
# For observation i and evaluation point b on the boundary:
#   s_i     signed network distance from i to the boundary (the metric distance
#           to the nearest boundary node c_i, plus snapping depth), negative on
#           the control side; it does not depend on b.
#   o_i(b)  signed lateral offset along the boundary: the metric distance from
#           b to c_i, signed by whether c_i lies before or after b on the line.
#   d_i(b)  network distance from b to i, used for kernel weights.
# (s_i, o_i(b)) play the role of (X1 - b1, X2 - b2) in rd2d's location-based
# local polynomial.

#' Network coordinates around a boundary
#'
#' Translates point locations into network coordinates around a boundary:
#' the signed network distance `s` to the boundary (positive on the treated
#' side) and, at every evaluation point `b` along the boundary, the network
#' distance `d(b)` used for kernel weights and the signed offset `o(b)` along
#' the boundary. These replace rd2d's Euclidean `(X1 - b1, X2 - b2)`.
#'
#' Each observation snaps to the street network **on its own side** of the
#' boundary; its off-network depth is charged at the network's `snap_scale`.
#' Repeated locations (panels, sensors) are processed once.
#'
#' @param network A [network_graph()], carrying its metric.
#' @param boundary A [boundary()].
#' @param points Observation locations: an `n x 2` matrix, a data frame with
#'   `x`/`y` columns (or exactly two columns), or sf points, in the network's
#'   projected coordinates (or longitude/latitude with `lonlat = TRUE`). sf
#'   points in another CRS are transformed to the network's `crs`.
#' @param treated Optional treatment indicator (logical or 0/1); defaults to
#'   the side of the boundary each point is on (points on the line count as
#'   treated with a `treated_region`).
#' @param eval_points Arc-length positions along the boundary, or a two-column
#'   matrix of coordinates to project onto it. Defaults to `n_eval` equally
#'   spaced points spanning the boundary nodes (around the whole ring for a
#'   closed boundary). Each evaluation point is placed
#'   on the street segment running along the boundary, or at the nearest
#'   boundary node when there is none.
#' @param n_eval Number of default evaluation points.
#' @param snap `"edge"` (attach to the nearest street segment) or `"node"`.
#' @param lonlat Points are longitude/latitude; project them into the
#'   network's `crs` (requires \pkg{sf}).
#' @param split_edges Insert nodes where edges cross the boundary (recommended).
#' @param extensive Edge columns divided in proportion when an edge is split
#'   at the boundary (the metric's own cost column always is); others are copied.
#'
#' @details On a directed network (`direction = "out"`), distances run from
#'   the boundary to each location, following edge directions: on a river
#'   drawn downstream, locations upstream of the boundary cannot be reached.
#'   Use `direction = "in"` for distances from each location to the boundary,
#'   or `"both"` to ignore direction.
#'
#' @return An object of class `rd2dnet_coords`, a list with (among others)
#'   `s`, `treated`, `xy`, `crossing_arc` (where each observation reaches the
#'   boundary), `eval_arc`, `eval_node_arc`, `eval_xy`, `eval_masked`, `loc`
#'   (location id of each observation), `n_obs`, `n_eval`, `n_locations`, and
#'   the (boundary-split) `network` and `boundary`. Use [coords_at()] for `d`
#'   and `o` at an evaluation point, [eval_table()] for the evaluation points,
#'   and `as.data.frame()` for a per-observation table.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 11)
#' coords
#' head(as.data.frame(coords, j = 5))
#' @export
network_coordinates <- function(network, boundary, points, treated = NULL, eval_points = NULL,
                                n_eval = 25, snap = "edge", lonlat = FALSE, split_edges = TRUE,
                                extensive = c("length", "travel_time")) {
  .check_network(network)
  .check_boundary(boundary)
  if (is.null(eval_points) && (length(n_eval) != 1 || !is.finite(n_eval) || n_eval < 1 || n_eval != round(n_eval))) {
    stop(sprintf("n_eval must be a positive integer, got %s", format(n_eval)), call. = FALSE)
  }
  xy <- .as_points(points, network$crs)
  if (lonlat) {
    if (is.null(network$crs)) stop("network has no crs to project lon/lat into", call. = FALSE)
    xy <- project_lonlat(xy[, 1], xy[, 2], network$crs)
  }
  if (!all(is.finite(xy))) stop("point coordinates must be finite", call. = FALSE)
  if (split_edges) network <- with_boundary_crossings(network, boundary, extensive)

  bdist <- boundary_distance(boundary, network$coords)
  node_on <- bdist <= boundary$tol
  near_tol <- 1e-4 * boundary$length
  near <- bdist > boundary$tol & bdist <= near_tol
  if (sum(near) >= 3) {
    warning(sprintf(paste0("%d network nodes lie within %.3g of the boundary but outside boundary(tol = %g); ",
                           "if the boundary follows a street, increase tol (in the network's coordinate ",
                           "units) so they count as on it"), sum(near), near_tol, boundary$tol), call. = FALSE)
  }
  node_T <- boundary_is_treated(boundary, network$coords) & !node_on
  bnodes <- which(node_on)
  if (!length(bnodes)) {
    stop("no network nodes lie on the boundary; check the coordinates or increase boundary(tol = ...)",
         call. = FALSE)
  }
  Tin <- if (is.null(treated)) NULL else .as_treated(treated, nrow(xy))
  # Repeated locations (panels, mass points) are processed once.
  ul <- .unique_locations(xy, Tin)
  loc <- ul$loc
  xy_l <- xy[ul$first, , drop = FALSE]
  T_l <- if (is.null(Tin)) boundary_is_treated(boundary, xy_l) else Tin[ul$first]
  if (!is.null(Tin)) {
    off <- T_l != boundary_is_treated(boundary, xy_l) & boundary_distance(boundary, xy_l) > boundary$tol
    if (any(off)) {
      n_rows <- sum(off[loc])
      warning(sprintf(paste0("%d observations (%.1f%%) have a treatment status that does not match their side ",
                             "of the boundary; they are attached to the network on the side `treated` gives ",
                             "(check the boundary's treated side, or data errors)"),
                      n_rows, 100 * n_rows / length(loc)), call. = FALSE)
    }
  }

  # Snap each location to the network on its own side of the boundary.
  n <- nrow(xy_l)
  u <- v <- integer(n)
  cu <- cv <- depth <- numeric(n)
  sides <- list(list(T_l, node_T | node_on), list(!T_l, !node_T | node_on))
  for (sd in sides) {
    m <- sd[[1]]
    if (any(m)) {
      sn <- snap_points(network, xy_l[m, , drop = FALSE], allowed = sd[[2]], mode = snap)
      u[m] <- sn$u
      v[m] <- sn$v
      cu[m] <- sn$cost_u
      cv[m] <- sn$cost_v
      depth[m] <- sn$depth
    }
  }
  snapped <- .snap(u, v, cu, cv, depth)
  depth_cost <- depth * .snap_rates(network, snapped)
  extent <- sqrt(sum(apply(network$coords, 2, function(z) diff(range(z)))^2))
  if (length(depth) && extent > 0 && stats::median(depth) > extent) {
    warning(sprintf(paste0("the points lie far from the network (median distance %.4g, network extent %.4g); ",
                           "check that they are in the network's CRS (lon/lat points need lonlat = TRUE)"),
                    stats::median(depth), extent), call. = FALSE)
  }

  bnode_arc <- boundary_project(boundary, network$coords[bnodes, , drop = FALSE])
  closed <- isTRUE(boundary$closed)
  if (is.null(eval_points)) {
    if (closed) { # around the ring, without repeating the start
      eval_arc <- (min(bnode_arc) + (seq_len(n_eval) - 1) * boundary$length / n_eval) %% boundary$length
    } else {
      eval_arc <- seq(min(bnode_arc), max(bnode_arc), length.out = n_eval)
    }
  } else {
    if (is.matrix(eval_points) || is.data.frame(eval_points) || .is_sf(eval_points)) {
      eval_arc <- boundary_project(boundary, eval_points)
    } else {
      eval_arc <- as.numeric(eval_points)
    }
    if (!length(eval_arc)) stop("eval_points is empty", call. = FALSE)
    if (!all(is.finite(eval_arc))) stop("eval_points must be finite", call. = FALSE)
    outside <- eval_arc < min(bnode_arc) - 1e-9 | eval_arc > max(bnode_arc) + 1e-9
    if (any(outside) && !closed) {
      warning(sprintf(paste0("%d evaluation points lie beyond the first or last place the network meets ",
                             "the boundary; they are moved to the nearest boundary node"), sum(outside)),
              call. = FALSE)
    }
  }
  eval_snap <- .snap_eval_points(network, boundary, eval_arc, node_on, bnodes, bnode_arc)

  if (.needs_point_nodes(network)) {
    # One-way edges and class-change penalties: put every location and
    # evaluation point on a node, so that paths respect both.
    k_obs <- length(u)
    pu <- c(u, eval_snap$u)
    pv <- c(v, eval_snap$v)
    pc <- c(cu, eval_snap$cost_u)
    tot <- pc + c(cv, eval_snap$cost_v)
    frac <- ifelse(tot > 0, pc / tot, 0)
    on_edge_ends <- .on_along(network, boundary, node_on, pu, pv)
    sp <- .with_point_nodes(network, pu, pv, frac, extensive)
    network <- sp$net
    pnode <- sp$node
    node_on <- boundary_distance(boundary, network$coords) <= boundary$tol
    node_on[pnode[on_edge_ends]] <- TRUE
    node_T <- boundary_is_treated(boundary, network$coords) & !node_on
    bnodes <- which(node_on)
    u <- v <- pnode[seq_len(k_obs)]
    cu <- cv <- numeric(k_obs)
    snapped <- .snap(u, v, cu, cv, depth)
    e <- pnode[-seq_len(k_obs)]
    eval_snap <- .snap(e, e, numeric(length(e)), numeric(length(e)), numeric(length(e)))
  }

  ns <- nearest_source(network, bnodes)
  node_dist <- ns$dist
  node_src <- ns$src
  via_u <- cu + node_dist[u]
  via_v <- cv + node_dist[v]
  dist <- depth_cost + pmin(via_u, via_v)
  crossing <- ifelse(via_u <= via_v, node_src[u], node_src[v])
  # A location snapped onto a segment that runs along the boundary is already there.
  on_edge <- .on_along(network, boundary, node_on, u, v)
  dist[on_edge] <- depth_cost[on_edge]
  crossing[on_edge] <- ifelse(cu <= cv, u, v)[on_edge]
  unreachable <- !is.finite(dist) | is.na(crossing)
  if (any(unreachable)) {
    warning(sprintf(paste0("%d observations (%d locations) have no network path to the boundary and are ",
                           "dropped (disconnected parts of the network, or edge directions)"),
                    sum(unreachable[loc]), sum(unreachable)), call. = FALSE)
  }
  s <- ifelse(unreachable, NA_real_, ifelse(T_l, dist, -dist))
  crossing[unreachable] <- NA_integer_
  safe <- ifelse(is.na(crossing), 1L, crossing)
  cross_loc <- .snap(ifelse(on_edge, u, safe), ifelse(on_edge, v, safe), ifelse(on_edge, cu, 0),
                     ifelse(on_edge, cv, 0), numeric(n))
  cross_arc <- ifelse(on_edge, boundary_project(boundary, xy_l),
                      boundary_project(boundary, network$coords[safe, , drop = FALSE]))
  cross_arc[unreachable] <- NA_real_

  masked <- boundary_is_masked(boundary, eval_arc)
  spot <- paste(pmin(eval_snap$u, eval_snap$v), pmax(eval_snap$u, eval_snap$v),
                round(ifelse(eval_snap$u <= eval_snap$v, eval_snap$cost_u, eval_snap$cost_v), 6))
  n_spots <- length(unique(spot[!masked]))
  if (n_spots < 0.75 * sum(!masked)) { # a few repeats near the ends are harmless
    warning(sprintf(paste0("the %d evaluation points fall on only %d distinct places (the network meets the ",
                           "boundary at few points); their estimates repeat, so pass one evaluation point ",
                           "per crossing, e.g. eval_points = <crossing coordinates>"), sum(!masked), n_spots),
            call. = FALSE)
  }
  ecu <- eval_snap$cost_u
  ecv <- eval_snap$cost_v
  tot <- ecu + ecv
  frac <- ifelse(tot > 0, ecu / tot, 0)
  eval_xy <- network$coords[eval_snap$u, , drop = FALSE] * (1 - frac) +
    network$coords[eval_snap$v, , drop = FALSE] * frac

  structure(
    list(
      network = network,
      boundary = boundary,
      xy = xy,
      treated = T_l[loc],
      snap = .snap_subset(snapped, loc),
      depth_cost = depth_cost[loc],
      s = s[loc],
      crossing_node = crossing[loc],
      crossing_arc = cross_arc[loc],
      crossing = .snap_subset(cross_loc, loc),
      node_treated = node_T,
      node_on_boundary = node_on,
      node_dist = node_dist,
      node_source = node_src,
      eval_arc = eval_arc,
      eval_snap = eval_snap,
      eval_node = ifelse(ecu <= ecv, eval_snap$u, eval_snap$v),
      eval_xy = eval_xy,
      eval_node_arc = boundary_project(boundary, eval_xy),
      eval_masked = masked,
      loc = loc,
      first = ul$first,
      lsnap = snapped,
      ldepth = depth_cost,
      lcross = cross_loc,
      lcross_arc = cross_arc,
      lcross_node = crossing,
      n_obs = nrow(xy),
      n_eval = length(eval_arc),
      n_locations = n,
      cache = new.env(parent = emptyenv())
    ),
    class = "rd2dnet_coords"
  )
}

#' @export
print.rd2dnet_coords <- function(x, ...) {
  cat(sprintf("<rd2dnet coordinates> n=%d, treated=%d, eval_points=%d, masked=%d, locations=%d\n",
              x$n_obs, sum(x$treated), x$n_eval, sum(x$eval_masked), x$n_locations))
  invisible(x)
}

.check_coords <- function(coords) {
  if (!inherits(coords, "rd2dnet_coords")) {
    stop("expected coordinates created with network_coordinates()", call. = FALSE)
  }
  coords
}

# Logical treatment indicator from logical or 0/1 input, with checks.
.as_treated <- function(treated, n) {
  t <- if (is.factor(treated)) as.character(treated) else treated
  if (!is.logical(t)) {
    v <- suppressWarnings(as.numeric(t))
    if (anyNA(v) || !all(v %in% c(0, 1))) {
      stop("treated must be logical or 0/1 without missing values", call. = FALSE)
    }
    t <- v == 1
  }
  if (anyNA(t)) stop("treated must be logical or 0/1 without missing values", call. = FALSE)
  if (length(t) != n) {
    stop(sprintf("treated must have one entry per point (%d), got %d", n, length(t)), call. = FALSE)
  }
  as.vector(t)
}

# `(loc, first)`: location id of each point and one point per location (by
# x, y and side), numbered in order of first appearance.
.unique_locations <- function(xy, T = NULL) {
  n <- nrow(xy)
  if (n == 0) {
    return(list(loc = integer(), first = integer()))
  }
  ord <- if (is.null(T)) order(xy[, 1], xy[, 2]) else order(xy[, 1], xy[, 2], T)
  a <- xy[ord, , drop = FALSE]
  new <- rep(TRUE, n)
  if (n > 1) {
    new[-1] <- a[-1, 1] != a[-n, 1] | a[-1, 2] != a[-n, 2]
    if (!is.null(T)) new[-1] <- new[-1] | T[ord][-1] != T[ord][-n]
  }
  loc <- integer(n)
  loc[ord] <- cumsum(new)
  first <- ord[new]
  srt <- order(first)
  rank <- integer(length(first))
  rank[srt] <- seq_along(first)
  list(loc = rank[loc], first = first[srt])
}

# Place evaluation points on boundary street segments (or boundary nodes).
# Edge-frame rows of street segments that run along the boundary. Both ends on
# the boundary is not enough: a road that leaves the boundary and comes back
# (a crescent, a loop into one side) is not on it. Its quarter points must lie
# within the boundary's tolerance too.
.along_boundary <- function(net, b, node_on) {
  E <- edge_frame(net)
  out <- node_on[E$u] & node_on[E$v] & (E$u != E$v)
  for (k in which(out)) {
    sg <- .segments(E$geometry[[k]])
    q <- .line_interpolate(sg, c(0.25, 0.5, 0.75) * sg$total)
    out[k] <- all(boundary_distance(b, q) <= b$tol)
  }
  out
}

# Whether points snapped between nodes u and v lie on a segment along the boundary.
.on_along <- function(net, b, node_on, u, v) {
  on <- node_on[u] & node_on[v] & (u != v)
  if (!any(on)) return(on)
  E <- edge_frame(net)
  keys <- paste(pmin(E$u, E$v), pmax(E$u, E$v))[.along_boundary(net, b, node_on)]
  on & paste(pmin(u, v), pmax(u, v)) %in% keys
}

.snap_eval_points <- function(net, b, eval_arc, node_on, bnodes, bnode_arc) {
  k <- length(eval_arc)
  L <- b$length
  closed <- isTRUE(b$closed)
  nearest <- bnodes[vapply(eval_arc, function(a) {
    gap <- abs(bnode_arc - a)
    if (closed) gap <- pmin(gap, L - gap) # the shorter way round the ring
    which.min(gap)
  }, integer(1))]
  out <- .snap(nearest, nearest, numeric(k), numeric(k), numeric(k))
  if (!isTRUE(net$metric$supports_edge_snapping)) {
    return(out)
  }
  E <- edge_frame(net)
  along <- .along_boundary(net, b, node_on) & is.finite(E$cost)
  if (!any(along)) {
    return(out)
  }
  E <- E[along, , drop = FALSE]
  au <- boundary_project(b, net$coords[E$u, , drop = FALSE])
  av <- boundary_project(b, net$coords[E$v, , drop = FALSE])
  # signed arc length from u to v; on a ring the shorter way, so a street across
  # the ring's start (arcs near 0 and near L) spans the seam, not the ring
  span <- av - au
  if (closed) span <- (span + L / 2) %% L - L / 2
  for (j in seq_len(k)) {
    along_u <- eval_arc[j] - au
    if (closed) along_u <- ifelse(span >= 0, along_u %% L, -((-along_u) %% L))
    frac_all <- ifelse(span != 0, along_u / span, -1)
    inside <- which(frac_all >= 0 & frac_all <= 1 & span != 0)
    if (!length(inside)) next
    r <- inside[1]
    frac <- frac_all[r]
    out$u[j] <- E$u[r]
    out$v[j] <- E$v[r]
    out$cost_u[j] <- frac * E$cost[r]
    out$cost_v[j] <- (1 - frac) * E$cost[r]
  }
  out
}

# Metric distance from evaluation point j to every network node (cached).
.node_row <- function(coords, j) {
  key <- as.character(j)
  row <- coords$cache[[key]]
  if (is.null(row)) {
    .node_rows(coords, j)
    row <- coords$cache[[key]]
  }
  row
}

# Compute and cache the rows of several evaluation points in one search.
.node_rows <- function(coords, js) {
  js <- js[!vapply(as.character(js), exists, logical(1), envir = coords$cache, inherits = FALSE)]
  if (!length(js)) {
    return(invisible())
  }
  e <- coords$eval_snap
  D <- distances_from(coords$network, c(e$u[js], e$v[js]))
  k <- length(js)
  for (r in seq_len(k)) {
    j <- js[r]
    assign(as.character(j), pmin(e$cost_u[j] + D[r, ], e$cost_v[j] + D[k + r, ]), envir = coords$cache)
  }
  invisible()
}

# Direct distance for locations on the same segment as evaluation point j (else Inf).
.along_eval_edge <- function(coords, loc, j) {
  e <- coords$eval_snap
  u0 <- e$u[j]
  v0 <- e$v[j]
  out <- rep(Inf, length(loc$u))
  if (u0 == v0) {
    return(out)
  }
  same <- loc$u == u0 & loc$v == v0
  rev <- loc$u == v0 & loc$v == u0
  out[same] <- abs(loc$cost_u[same] - e$cost_u[j])
  out[rev] <- abs(loc$cost_u[rev] - e$cost_v[j])
  out
}

#' Kernel distance and boundary offset at an evaluation point
#'
#' @param coords Output of [network_coordinates()].
#' @param j Evaluation point (1-based index into `coords$eval_arc`).
#' @return A list with `d`, the network distance from evaluation point `j`
#'   to every observation (for kernel weights), and `o`, each observation's
#'   signed offset along the boundary from `j` (`NA` when the boundary cannot
#'   be reached).
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' str(coords_at(coords, 3))
#' @export
coords_at <- function(coords, j) {
  .check_coords(coords)
  j <- as.integer(j)
  if (length(j) != 1 || is.na(j) || j < 1 || j > coords$n_eval) {
    stop("j must be a single evaluation point index", call. = FALSE)
  }
  D <- .node_row(coords, j)
  sn <- coords$lsnap
  d <- coords$ldepth + pmin(sn$cost_u + D[sn$u], sn$cost_v + D[sn$v])
  d <- pmin(d, coords$ldepth + .along_eval_edge(coords, sn, j))
  diff <- coords$lcross_arc - coords$eval_node_arc[j]
  if (isTRUE(coords$boundary$closed)) { # the shorter way round the ring
    L <- coords$boundary$length
    diff <- (diff + L / 2) %% L - L / 2
  }
  sgn <- sign(diff)
  cr <- coords$lcross
  reach <- pmin(cr$cost_u + D[cr$u], cr$cost_v + D[cr$v])
  reach <- pmin(reach, .along_eval_edge(coords, cr, j))
  o <- ifelse(is.na(coords$lcross_node), NA_real_, sgn * reach)
  list(d = d[coords$loc], o = o[coords$loc])
}

#' Evaluation points of network coordinates
#'
#' @param coords Output of [network_coordinates()].
#' @return A data frame with one row per evaluation point: requested arc
#'   position `arc`, nearest network node id `node`, arc position after
#'   snapping onto the network `node_arc`, coordinates `x`, `y`, and `masked`.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' eval_table(coords)
#' @export
eval_table <- function(coords) {
  .check_coords(coords)
  data.frame(arc = coords$eval_arc, node = coords$network$nodes$name[coords$eval_node],
             node_arc = coords$eval_node_arc, x = coords$eval_xy[, 1], y = coords$eval_xy[, 2],
             masked = coords$eval_masked)
}

#' @rdname network_coordinates
#' @param x Output of [network_coordinates()].
#' @param row.names,optional Ignored (for compatibility with the generic).
#' @param j Optional evaluation point index; adds the columns `d` and `o` at that point.
#' @param ... Ignored.
#' @export
as.data.frame.rd2dnet_coords <- function(x, row.names = NULL, optional = FALSE, j = NULL, ...) {
  out <- data.frame(x = x$xy[, 1], y = x$xy[, 2], treated = x$treated, s = x$s,
                    crossing_arc = x$crossing_arc, depth = x$snap$depth)
  if (!is.null(j)) {
    at <- coords_at(x, j)
    out$d <- at$d
    out$o <- at$o
  }
  out
}

#' Network distances between observations
#'
#' @param coords Output of [network_coordinates()].
#' @param idx Observation indices for the rows.
#' @param cols Observation indices for the columns (default `idx`).
#' @param cutoff Distances beyond `cutoff` may be reported as `Inf`.
#' @return A `length(idx) x length(cols)` matrix of network distances
#'   (including both snapping depths); an observation is at distance 0 from itself.
#' @examples
#' city <- load_mock_midtown()
#' coords <- network_coordinates(city$network(), city$boundary, city$obs, n_eval = 5)
#' round(pairwise_distances(coords, 1:4), 1)
#' @export
pairwise_distances <- function(coords, idx, cols = NULL, cutoff = Inf) {
  .check_coords(coords)
  idx <- as.integer(idx)
  cols <- if (is.null(cols)) idx else as.integer(cols)
  sn <- coords$snap
  row_nodes <- sort(unique(c(sn$u[idx], sn$v[idx])))
  Dn <- distances_from(coords$network, row_nodes, limit = cutoff)
  iu <- match(sn$u[idx], row_nodes)
  iv <- match(sn$v[idx], row_nodes)
  best <- matrix(Inf, length(idx), length(cols))
  for (rr in list(list(iu, sn$cost_u[idx]), list(iv, sn$cost_v[idx]))) {
    Drow <- Dn[rr[[1]], , drop = FALSE]
    for (cc in list(list(sn$u[cols], sn$cost_u[cols]), list(sn$v[cols], sn$cost_v[cols]))) {
      best <- pmin(best, outer(rr[[2]], cc[[2]], "+") + Drow[, cc[[1]], drop = FALSE])
    }
  }
  dep <- outer(coords$depth_cost[idx], coords$depth_cost[cols], "+")
  out <- best + dep
  same_edge <- outer(sn$u[idx], sn$u[cols], "==") & outer(sn$v[idx], sn$v[cols], "==")
  along <- abs(outer(sn$cost_u[idx], sn$cost_u[cols], "-")) + dep
  out[same_edge] <- pmin(out, along)[same_edge]
  out[outer(idx, cols, "==")] <- 0
  out
}

# Conley (Bartlett) HAC quadratic form Psi' K Psi for observations `idx`, with
# K_ik = max(0, 1 - d_ik / cutoff). Rows of Psi are first summed by location
# (K depends on location only); only pairs closer than the cutoff are visited.
.hac_quadratic <- function(coords, idx, Psi, cutoff, distance = "network") {
  if (!distance %in% c("network", "euclidean")) {
    stop("distance must be 'network' or 'euclidean'", call. = FALSE)
  }
  Psi <- matrix(as.numeric(Psi), nrow = length(idx))
  lid <- coords$loc[idx]
  locs <- sort(unique(lid))
  P <- rowsum(Psi, match(lid, locs), reorder = TRUE)
  if (distance == "euclidean") {
    pts <- coords$xy[coords$first[locs], , drop = FALSE]
    out <- .hac_euclidean(pts[, 1], pts[, 2], P, cutoff)
  } else {
    sn <- coords$lsnap
    eu <- sn$u[locs]
    ev <- sn$v[locs]
    nodes <- sort(unique(c(eu, ev)))
    pr <- .engine_pairs(.engine(coords$network), nodes, cutoff)
    tgt <- match(pr$node, nodes)
    keep <- !is.na(tgt)
    out <- .hac_network(match(eu, nodes), match(ev, nodes), sn$cost_u[locs], sn$cost_v[locs],
                        coords$ldepth[locs], pr$i[keep], tgt[keep], pr$dist[keep], length(nodes),
                        P, cutoff)
  }
  (out + t(out)) / 2
}
