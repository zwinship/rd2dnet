#' Boundary discontinuity estimates with the data entry of rd2d
#'
#' `rd2dnet(Y, X, t, b, network)` mirrors `rd2d(Y, X, t, b)` from the rd2d
#' package: the same outcome, locations, treatment indicator and boundary
#' points, plus the network on which distance is measured. It builds the
#' boundary and the network coordinates and calls [rdnet()].
#'
#' @param Y Outcome, length `n`.
#' @param X Locations: an `n x 2` matrix, a data frame with `x`/`y` (or
#'   exactly two columns), or sf points, in the network's projected coordinates.
#' @param t Treatment indicator (0/1 or logical), length `n`.
#' @param b Evaluation points along the boundary, a `J x 2` matrix in order
#'   along it.
#' @param network A network from [network_graph()] or one of its
#'   constructors, an igraph graph with `x`/`y` vertex attributes, or line
#'   geometries accepted by [network_graph_from_lines()] (an sf layer of
#'   LINESTRINGs or a data frame of segments with `x1, y1, x2, y2`).
#' @param boundary The boundary line. By default it is the line through `b`,
#'   with the treated side taken from where the treated observations lie
#'   ([boundary_from_points()]). Pass a [boundary()] when the line extends
#'   beyond the first and last evaluation points or a treated region is known.
#' @param treated_side `"left"` or `"right"` of the line through `b`, instead
#'   of inferring it from `X` and `t`.
#' @param metric Distance on the network when `network` is not already an
#'   `rd2dnet_network`.
#' @param snap `"edge"` or `"node"`, passed to [network_coordinates()].
#' @param lonlat `X` and `b` are longitude/latitude; they are projected into
#'   the network's CRS (requires sf).
#' @param ... Passed to [rdnet()]: `h`, `p`, `kernel`, `bwselect`, `vce`,
#'   `cluster`, `fuzzy`, `covs`, `level`, ...
#'
#' @return An `rd2dnet_result` (see [rdnet()]) with estimates at each row of
#'   `b`; its `coords` element holds the network coordinates for [density()]
#'   and [rdbwnet()].
#' @seealso [rdnet()] for the step-by-step interface.
#' @examples
#' city <- load_mock_midtown()
#' b <- boundary_interpolate(city$boundary, boundary_grid(city$boundary, 7, start = 200, end = 2800))
#' fit <- rd2dnet(city$obs$outcome, city$obs[, c("x", "y")], city$obs$treated, b,
#'                city$network(), boundary = city$boundary, h = 900, cbands = FALSE)
#' fit
#' @export
rd2dnet <- function(Y, X, t, b, network, boundary = NULL, treated_side = NULL, metric = "length",
                    snap = "edge", lonlat = FALSE, ...) {
  net <- .as_network(network, metric)
  xy <- .as_points(X, net$crs)
  b <- .as_points(b)
  if (lonlat) {
    if (is.null(net$crs)) stop("network has no crs to project lon/lat into", call. = FALSE)
    xy <- project_lonlat(xy[, 1], xy[, 2], net$crs)
    b <- project_lonlat(b[, 1], b[, 2], net$crs)
  }
  t <- .as_treatment(t, nrow(xy))
  if (is.null(boundary)) boundary <- boundary_from_points(b, xy, t, treated_side = treated_side)
  coords <- network_coordinates(net, boundary, xy, treated = t, eval_points = b, snap = snap)
  rdnet(Y, coords, ...)
}

.as_treatment <- function(t, n) {
  if (is.factor(t)) t <- as.character(t)
  v <- suppressWarnings(as.numeric(t))
  if (anyNA(v) || !all(v %in% c(0, 1))) {
    stop("t must be a 0/1 (or logical) treatment indicator without missing values", call. = FALSE)
  }
  if (length(v) != n) stop("X and t must have the same number of observations", call. = FALSE)
  v == 1
}

.as_network <- function(network, metric = "length") {
  if (inherits(network, "rd2dnet_network")) return(network)
  if (inherits(network, "igraph")) return(network_graph_from_igraph(network, metric = metric))
  if (is.function(network)) stop("network is a function; call it, e.g. city$network()", call. = FALSE)
  network_graph_from_lines(network, metric = metric)
}

#' Boundaries from evaluation points or from a treated region
#'
#' `boundary_from_points()` draws the boundary through the ordered evaluation
#' points `b` of rd2d and infers the treated side from the observations.
#' `boundary_from_region()` uses the outline of the treated region (a school
#' zone, a city, a low-emission zone) as the boundary.
#'
#' @param b A `J x 2` matrix of points along the boundary, in order.
#' @param X,t Observation locations and treatment indicators, used to infer
#'   the treated side: the side of the line on which most treated observations
#'   lie. A warning is given when fewer than 90% of observations agree.
#' @param treated_side `"left"` or `"right"`, instead of inferring it.
#' @param region The treated region: a ring matrix, a list of rings, or an sf
#'   polygon. Points inside it are treated.
#' @param within Keep only the part of the outline inside this polygon (same
#'   formats), for example the stretch where the zone meets its neighbour
#'   rather than the coast or the edge of the study area. The kept part must be
#'   one connected line.
#' @param ... Passed to [boundary()] (`mask`, `tol`).
#' @return A [boundary()].
#' @examples
#' X <- cbind(runif(200, 0, 100), runif(200, -50, 50))
#' t <- X[, 2] < 0
#' boundary_from_points(rbind(c(0, 0), c(50, 0), c(100, 0)), X, t)
#'
#' zone <- rect_polygon(1000, -1000, 3000, 1000)
#' boundary_from_region(zone, within = rect_polygon(-1e6, -1e6, 1500, 1e6))
#' @export
boundary_from_points <- function(b, X = NULL, t = NULL, treated_side = NULL, ...) {
  b <- .as_points(b)
  keep <- c(TRUE, rowSums(abs(diff(b))) > 0)
  b <- b[keep, , drop = FALSE]
  if (nrow(b) < 2) stop("b needs at least two distinct points to define a boundary line", call. = FALSE)
  if (is.null(treated_side)) {
    if (is.null(X) || is.null(t)) stop("give treated_side, or X and t to infer it", call. = FALSE)
    treated_side <- .infer_side(b, .as_points(X), t)
  }
  boundary(b, treated_side = treated_side, ...)
}

.infer_side <- function(L, X, t) {
  t <- .as_treatment(t, nrow(X))
  if (all(t) || !any(t)) stop("t must contain both treated and control observations", call. = FALSE)
  left <- boundary_is_treated(boundary(L, treated_side = "left"), X)
  agree <- mean(left == t)
  side <- if (agree >= 0.5) "left" else "right"
  share <- max(agree, 1 - agree)
  if (share < 0.9) {
    warning(sprintf(paste0("only %.0f%% of observations lie on the side of the boundary that matches t; ",
                           "check that b is ordered along the boundary or pass treated_side/boundary"),
                    100 * share), call. = FALSE)
  }
  side
}

#' @rdname boundary_from_points
#' @export
boundary_from_region <- function(region, within = NULL, ...) {
  g <- .as_geom(region)
  if (g$type != "polygon") stop("region must be a polygon", call. = FALSE)
  if (length(g$parts) != 1 && is.null(within)) {
    stop("the region has several rings; pick one stretch of its outline with within =", call. = FALSE)
  }
  pieces <- list()
  for (ring in g$parts) {
    if (ring[1, 1] != ring[nrow(ring), 1] || ring[1, 2] != ring[nrow(ring), 2]) ring <- rbind(ring, ring[1, ])
    pieces <- c(pieces, if (is.null(within)) list(ring) else .clip_line(ring, .as_geom(within)))
  }
  pieces <- Filter(function(p) nrow(p) >= 2, pieces)
  if (!length(pieces)) stop("no part of the outline lies within `within`", call. = FALSE)
  line <- .join_lines(pieces)
  if (is.null(line)) stop("the boundary must be a single connected line; narrow it with within =", call. = FALSE)
  boundary(line, treated_region = region, ...)
}

# Pieces of a closed ring inside a polygon: split every segment at its
# crossings with the polygon outline and keep the parts whose midpoint is inside.
.clip_line <- function(L, poly) {
  ps <- .geom_segments(poly)
  pts <- list(L[1, , drop = FALSE])
  for (k in seq_len(nrow(L) - 1)) {
    p0 <- L[k, ]
    p1 <- L[k + 1, ]
    d <- p1 - p0
    den <- d[1] * (ps$y1 - ps$y0) - d[2] * (ps$x1 - ps$x0)
    ok <- abs(den) > 1e-12
    s <- ((ps$x0 - p0[1]) * (ps$y1 - ps$y0) - (ps$y0 - p0[2]) * (ps$x1 - ps$x0)) / den
    u <- ((ps$x0 - p0[1]) * d[2] - (ps$y0 - p0[2]) * d[1]) / den
    s <- sort(unique(s[ok & s > 0 & s < 1 & u >= 0 & u <= 1]))
    for (si in s) pts[[length(pts) + 1]] <- matrix(p0 + si * d, 1)
    pts[[length(pts) + 1]] <- matrix(p1, 1)
  }
  P <- do.call(rbind, pts)
  mid <- (P[-1, , drop = FALSE] + P[-nrow(P), , drop = FALSE]) / 2
  inside <- .geom_contains(poly, mid) | .geom_distance(poly, mid) <= 1e-9 * max(1, abs(P))
  r <- rle(inside)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  lapply(which(r$values), function(i) P[starts[i]:(ends[i] + 1), , drop = FALSE])
}

# Join polylines that meet end to end into one line (NULL if they do not).
.join_lines <- function(pieces) {
  line <- pieces[[1]]
  rest <- pieces[-1]
  same <- function(a, b) sum(abs(a - b)) <= 1e-9 * max(1, abs(a))
  while (length(rest)) {
    hit <- FALSE
    for (i in seq_along(rest)) {
      p <- rest[[i]]
      n <- nrow(line)
      m <- nrow(p)
      if (same(line[n, ], p[1, ])) line <- rbind(line, p[-1, , drop = FALSE])
      else if (same(line[n, ], p[m, ])) line <- rbind(line, p[rev(seq_len(m - 1)), , drop = FALSE])
      else if (same(line[1, ], p[m, ])) line <- rbind(p[-m, , drop = FALSE], line)
      else if (same(line[1, ], p[1, ])) line <- rbind(p[rev(seq_len(m))[-m], , drop = FALSE], line)
      else next
      rest <- rest[-i]
      hit <- TRUE
      break
    }
    if (!hit) return(NULL)
  }
  line
}

#' Networks from line geometries
#'
#' Builds a network from a layer of lines, such as GIS roads, rivers, rail or
#' pipes. Line ends that meet (within `tol`) become shared nodes; each line
#' becomes an edge with its geometry and length.
#'
#' @param lines An sf object (or sfc) of `LINESTRING`/`MULTILINESTRING`, a list
#'   of `n x 2` vertex matrices, or a data frame of straight segments with
#'   endpoint columns `coords`. Other columns (speed, lanes, flow, ...) become
#'   edge attributes, so they can be used in a metric.
#' @param metric How distance is measured (see [network_graph()]).
#' @param tol Line ends closer than this are joined into one node. Increase it
#'   for layers whose lines do not quite meet.
#' @param split_at_intersections Also split lines where they cross each
#'   other, for layers that are not split at junctions (requires sf). Leave it
#'   off when crossings without a junction (bridges, tunnels) matter.
#' @param directed Lines are drawn in the direction of travel or flow
#'   (one-way streets, rivers). Combine with `direction = "out"` (downstream)
#'   or `"in"`.
#' @param coords Endpoint columns of a segment data frame.
#' @param to_crs For an sf layer in longitude/latitude: `"utm"` (default)
#'   projects it to its local UTM zone so distances are in metres; any other
#'   CRS projects to that; `NULL` keeps degrees.
#' @param crs CRS label of plain (non-sf) coordinates, if known.
#' @param ... Passed to [network_graph()] (`direction`, `snap_scale`).
#'
#' @return An `rd2dnet_network`. Each edge has `geometry` and `length` (the
#'   line's length, unless a `length` column is supplied). Closed rings are
#'   split in two so they can be routed through.
#' @examples
#' seg <- data.frame(x1 = c(0, 100, 100), y1 = c(0, 0, 0),
#'                   x2 = c(100, 200, 100), y2 = c(0, 0, 100), speed_kph = c(30, 30, 50))
#' net <- network_graph_from_lines(seg)
#' net
#' node_distance(net, 1, 3)
#' @export
network_graph_from_lines <- function(lines, metric = "length", tol = 1e-6, split_at_intersections = FALSE,
                                     directed = FALSE, coords = c("x1", "y1", "x2", "y2"), to_crs = "utm",
                                     crs = NULL, ...) {
  inp <- .lines_input(lines, coords, to_crs, crs, split_at_intersections)
  geoms <- inp$geoms
  attrs <- inp$attrs
  len <- vapply(geoms, function(m) .segments(m)$total, numeric(1))
  ok <- len > 0
  if (!all(ok)) {
    warning(sprintf("%d empty or zero-length lines were dropped", sum(!ok)), call. = FALSE)
    geoms <- geoms[ok]
    attrs <- attrs[ok, , drop = FALSE]
    len <- len[ok]
  }
  # closed rings: split at the midpoint so both halves are routable
  ring <- vapply(geoms, function(m) all(m[1, ] == m[nrow(m), ]), logical(1))
  if (any(ring)) {
    idx <- rep(seq_along(geoms), ifelse(ring, 2L, 1L))
    geoms <- unlist(lapply(seq_along(geoms), function(k) {
      if (ring[k]) list(.line_substring(geoms[[k]], 0, 0.5), .line_substring(geoms[[k]], 0.5, 1))
      else geoms[k]
    }), recursive = FALSE)
    attrs <- attrs[idx, , drop = FALSE]
    len <- vapply(geoms, function(m) .segments(m)$total, numeric(1))
  }
  m <- length(geoms)
  ends <- rbind(t(vapply(geoms, function(g) g[1, ], numeric(2))),
                t(vapply(geoms, function(g) g[nrow(g), ], numeric(2))))
  lab <- .merge_close(ends, tol)
  first <- !duplicated(lab)
  nodes <- data.frame(name = lab[first], x = ends[first, 1], y = ends[first, 2])
  nodes <- nodes[order(nodes$name), , drop = FALSE]
  edges <- data.frame(from = lab[seq_len(m)], to = lab[m + seq_len(m)])
  if (ncol(attrs)) edges <- cbind(edges, attrs[setdiff(names(attrs), c("from", "to", "geometry"))])
  if (!"length" %in% names(edges)) edges$length <- len
  edges$geometry <- geoms
  rownames(edges) <- NULL
  network_graph(edges, nodes, metric = metric, directed = directed, crs = inp$crs, ...)
}

.lines_input <- function(lines, coords, to_crs, crs, split) {
  if (.is_sf(lines)) {
    .need_sf()
    if (inherits(lines, "sfg")) lines <- sf::st_sfc(lines)
    g <- sf::st_geometry(lines)
    attrs <- if (inherits(lines, "sf")) sf::st_drop_geometry(lines) else data.frame(row.names = seq_along(g))
    if (!all(as.character(sf::st_geometry_type(g)) %in% c("LINESTRING", "MULTILINESTRING"))) {
      stop("every geometry must be a LINESTRING or MULTILINESTRING", call. = FALSE)
    }
    if (isTRUE(sf::st_is_longlat(g)) && !is.null(to_crs)) {
      target <- if (identical(to_crs, "utm")) .utm_crs(g) else to_crs
      g <- sf::st_transform(g, target)
    }
    lay <- sf::st_sf(attrs, geometry = g)
    lay <- suppressWarnings(sf::st_cast(lay, "MULTILINESTRING"))
    lay <- suppressWarnings(sf::st_cast(lay, "LINESTRING"))
    if (split) {
      noded <- sf::st_cast(sf::st_sf(geometry = sf::st_node(sf::st_union(sf::st_geometry(lay)))), "LINESTRING")
      mid <- suppressWarnings(sf::st_line_sample(sf::st_geometry(noded), sample = 0.5))
      src <- sf::st_nearest_feature(mid, lay)
      piece_len <- as.numeric(sf::st_length(noded))
      new_attrs <- sf::st_drop_geometry(lay)[src, , drop = FALSE]
      if ("length" %in% names(new_attrs)) {
        new_attrs$length <- new_attrs$length * piece_len / as.numeric(sf::st_length(lay))[src]
      }
      lay <- sf::st_sf(new_attrs, geometry = sf::st_geometry(noded))
    }
    crs_out <- sf::st_crs(lay)
    geoms <- lapply(sf::st_geometry(lay), function(l) unname(matrix(unclass(l)[, 1:2], ncol = 2)))
    attrs <- sf::st_drop_geometry(lay)
    rownames(attrs) <- NULL
    crs_lab <- if (is.na(crs_out)) crs else (crs_out$input %||% crs)
    if (!is.na(crs_out) && !is.null(crs_out$epsg) && !is.na(crs_out$epsg)) crs_lab <- paste0("EPSG:", crs_out$epsg)
    return(list(geoms = geoms, attrs = attrs, crs = crs_lab))
  }
  if (split) stop("split_at_intersections needs the lines as an sf object", call. = FALSE)
  if (is.data.frame(lines)) {
    miss <- setdiff(coords, names(lines))
    if (length(miss)) {
      stop(sprintf("segment table has no %s columns (set coords = to the endpoint columns)",
                   paste(miss, collapse = ", ")), call. = FALSE)
    }
    xy <- vapply(coords, function(cl) as.numeric(lines[[cl]]), numeric(nrow(lines)))
    xy <- matrix(xy, ncol = 4)
    geoms <- lapply(seq_len(nrow(xy)), function(k) matrix(xy[k, ], 2, 2, byrow = TRUE))
    attrs <- lines[setdiff(names(lines), coords)]
    return(list(geoms = geoms, attrs = as.data.frame(attrs), crs = crs))
  }
  if (is.list(lines)) {
    geoms <- lapply(lines, .as_points)
    return(list(geoms = geoms, attrs = data.frame(row.names = seq_along(geoms)), crs = crs))
  }
  stop("lines must be an sf layer, a list of vertex matrices or a data frame of segments", call. = FALSE)
}

# Local UTM zone of an sf geometry in longitude/latitude.
.utm_crs <- function(g) {
  bb <- sf::st_bbox(g)
  lon <- (bb[["xmin"]] + bb[["xmax"]]) / 2
  lat <- (bb[["ymin"]] + bb[["ymax"]]) / 2
  zone <- floor((lon + 180) / 6) %% 60 + 1
  paste0("EPSG:", (if (lat >= 0) 32600 else 32700) + zone)
}

# Label points so that points within `tol` of each other share a label
# (numbered in order of first appearance).
.merge_close <- function(xy, tol) {
  n <- nrow(xy)
  cell <- max(tol, 1e-12)
  cx <- floor(xy[, 1] / cell)
  cy <- floor(xy[, 2] / cell)
  cell_key <- function(a, b) paste(sprintf("%.0f", a), sprintf("%.0f", b))
  key <- cell_key(cx, cy)
  pts <- data.frame(i = seq_len(n), key = key)
  pairs <- list()
  for (o in list(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(1, -1))) {
    other <- data.frame(j = seq_len(n), key = cell_key(cx - o[1], cy - o[2]))
    pr <- merge(pts, other, by = "key")
    pairs[[length(pairs) + 1]] <- pr[pr$i != pr$j, c("i", "j")]
  }
  pr <- do.call(rbind, pairs)
  if (nrow(pr)) {
    d <- sqrt(rowSums((xy[pr$i, , drop = FALSE] - xy[pr$j, , drop = FALSE])^2))
    pr <- pr[d <= tol, , drop = FALSE]
  }
  lab <- seq_len(n)
  if (nrow(pr)) {
    repeat {
      new <- lab
      m <- pmin(new[pr$i], new[pr$j])
      idx <- c(pr$i, pr$j)
      o <- order(idx, c(m, m))
      first <- o[!duplicated(idx[o])]  # smallest neighbouring label of each point
      new[idx[first]] <- pmin(new[idx[first]], c(m, m)[first])
      new <- new[new]
      if (identical(new, lab)) break
      lab <- new
    }
  }
  match(lab, unique(lab))
}
