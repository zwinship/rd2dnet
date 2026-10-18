#' Spatial network with a distance metric
#'
#' A street (or any spatial) network together with the metric that defines
#' distance on it.
#'
#' @param edges Data frame with one row per edge: columns `from` and `to`
#'   (node ids, as in the node id column of `nodes`) plus any edge attributes
#'   (`length`, `highway`, `speed_kph`, ...). An optional `geometry` column
#'   holds each edge's shape as a list of `n x 2` matrices or an sf
#'   `LINESTRING` column; edges without one are straight segments. Repeated
#'   node pairs (multigraphs) are allowed; self-loops are ignored.
#' @param nodes Data frame with one row per node: a node id column (`name` or
#'   `id`; otherwise ids are the row numbers), projected coordinates `x` and
#'   `y`, and any node attributes.
#' @param metric How distance is measured: an edge column name (`"length"`),
#'   a vectorised cost function, a metric from [metrics] such as
#'   [transition_cost()], or a square node-by-node matrix.
#' @param directed Whether edges are one-way (`from` -> `to`).
#' @param direction For directed networks: `"both"` treats every edge as
#'   two-way (pedestrian access, the default); `"out"` follows edge direction;
#'   `"in"` reverses it.
#' @param x,y Names of the coordinate columns in `nodes`.
#' @param crs Coordinate reference system label of the projected coordinates
#'   (used to project longitude/latitude points).
#' @param snap_scale Metric cost per unit of off-network snapping distance.
#'   Defaults to the metric's own default (1 for length, 0 otherwise).
#'
#' @return An object of class `rd2dnet_network`.
#' @seealso [network_graph_from_adjacency()], [network_graph_from_igraph()],
#'   [network_graph_from_place()], [network_coordinates()].
#' @examples
#' nodes <- data.frame(name = c("a", "b", "c"), x = c(0, 100, 200), y = 0)
#' edges <- data.frame(from = c("a", "b"), to = c("b", "c"), length = c(100, 100))
#' net <- network_graph(edges, nodes)
#' net
#' node_distance(net, "a", "c")
#' @export
network_graph <- function(edges, nodes, metric = "length", directed = FALSE, direction = "both",
                          x = "x", y = "y", crs = NULL, snap_scale = NULL) {
  if (!direction %in% c("both", "out", "in")) {
    stop("direction must be one of 'both', 'out', 'in'", call. = FALSE)
  }
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  if (!all(c(x, y) %in% names(nodes))) {
    stop(sprintf("every node needs '%s' and '%s' coordinates", x, y), call. = FALSE)
  }
  if (!all(c("from", "to") %in% names(edges))) {
    stop("edges need `from` and `to` columns", call. = FALSE)
  }
  id_col <- intersect(c("name", "id"), names(nodes))
  ids <- if (length(id_col)) nodes[[id_col[1]]] else seq_len(nrow(nodes))
  if (anyDuplicated(ids)) stop("node ids must be unique", call. = FALSE)
  if (is.factor(ids)) ids <- as.character(ids)
  coords <- cbind(as.numeric(nodes[[x]]), as.numeric(nodes[[y]]))
  if (anyNA(coords)) stop(sprintf("every node needs '%s' and '%s' coordinates", x, y), call. = FALSE)
  attrs <- nodes[setdiff(names(nodes), c(id_col[1], x, y))]
  nodes <- data.frame(name = ids, x = coords[, 1], y = coords[, 2], stringsAsFactors = FALSE)
  if (ncol(attrs)) nodes <- cbind(nodes, attrs)
  from <- match(.as_id(edges$from, ids), ids)
  to <- match(.as_id(edges$to, ids), ids)
  if (anyNA(from) || anyNA(to)) {
    stop("every edge endpoint must be a node id", call. = FALSE)
  }
  edges$from <- from
  edges$to <- to
  if ("geometry" %in% names(edges)) {
    edges$geometry <- .edge_geometries(edges$geometry)
  }
  metric <- as_metric(metric)
  structure(
    list(
      nodes = nodes,
      edges = edges,
      directed = isTRUE(directed),
      direction = direction,
      metric = metric,
      crs = crs,
      snap_scale = if (is.null(snap_scale)) .default_snap_scale(metric) else as.numeric(snap_scale),
      coords = unname(coords),
      n_nodes = nrow(nodes),
      cache = new.env(parent = emptyenv())
    ),
    class = "rd2dnet_network"
  )
}

.as_id <- function(v, ids) {
  if (is.factor(v)) v <- as.character(v)
  if (is.character(ids) && !is.character(v)) v <- as.character(v)
  if (is.numeric(ids) && is.character(v)) v <- suppressWarnings(as.numeric(v))
  v
}

.edge_geometries <- function(g) {
  if (.is_sf(g)) {
    .need_sf()
    return(lapply(seq_along(g), function(k) unname(sf::st_coordinates(g[k])[, c("X", "Y"), drop = FALSE])))
  }
  lapply(g, function(m) if (is.null(m) || length(m) == 0) NULL else .as_points(m))
}

#' @export
print.rd2dnet_network <- function(x, ...) {
  cat(sprintf("<rd2dnet network> nodes=%d, edges=%d, metric=%s, direction=%s%s\n", x$n_nodes,
              nrow(x$edges), format(x$metric), if (x$directed) x$direction else "undirected",
              if (!is.null(x$crs)) paste0(", crs=", x$crs) else ""))
  invisible(x)
}

.check_network <- function(net) {
  if (!inherits(net, "rd2dnet_network")) {
    stop("expected a network created with network_graph()", call. = FALSE)
  }
  net
}

# Edges without self-loops (from/to are node indices).
.stored_edges <- function(net) {
  E <- net$edges
  E[E$from != E$to, , drop = FALSE]
}

# Shortest-path engine compiled from the metric (built on first use).
.engine <- function(net) {
  if (is.null(net$cache$engine)) {
    net$cache$engine <- .compile_metric(net$metric, net)
  }
  net$cache$engine
}

#' Network distances
#'
#' Metric distances between network nodes.
#'
#' @param net A [network_graph()].
#' @param sources Node indices (row numbers of `net$nodes`) to search from.
#' @param limit Distances beyond `limit` are reported as `Inf`.
#' @param a,b Node ids.
#' @return `distances_from()`: a `length(sources) x n_nodes` matrix.
#'   `nearest_source()`: a list with, for every node, the distance to the
#'   nearest source (`dist`) and that source's node index (`src`, `NA` if
#'   unreachable). `node_distance()`: a single distance.
#' @examples
#' nodes <- data.frame(name = 1:3, x = c(0, 1, 2), y = 0)
#' net <- network_graph(data.frame(from = 1:2, to = 2:3, length = c(5, 7)), nodes)
#' distances_from(net, 1)
#' nearest_source(net, c(1, 3))
#' node_distance(net, 1, 3)
#' @name network_distances
NULL

#' @rdname network_distances
#' @export
distances_from <- function(net, sources, limit = Inf) {
  .engine_distances(.engine(.check_network(net)), sources, limit)
}

#' @rdname network_distances
#' @export
nearest_source <- function(net, sources) {
  .engine_nearest(.engine(.check_network(net)), sources)
}

#' @rdname network_distances
#' @export
node_distance <- function(net, a, b) {
  .check_network(net)
  ia <- match(.as_id(a, net$nodes$name), net$nodes$name)
  ib <- match(.as_id(b, net$nodes$name), net$nodes$name)
  if (is.na(ia) || is.na(ib)) stop("unknown node id", call. = FALSE)
  distances_from(net, ia)[1, ib]
}

#' Street segments of a network
#'
#' One row per undirected street segment (the first edge of each node pair):
#' endpoint indices `u` and `v`, geometric `length` (the `length` column, or
#' the geometry's length), the metric `cost` of traversing it (`NA` for a
#' precomputed metric), and its `geometry` oriented from `u` to `v`.
#'
#' @param net A [network_graph()].
#' @return A data frame with a list column `geometry`.
#' @examples
#' city <- load_mock_midtown()
#' head(edge_frame(city$network())[, c("u", "v", "length", "cost")])
#' @export
edge_frame <- function(net) {
  .check_network(net)
  if (!is.null(net$cache$edges)) {
    return(net$cache$edges)
  }
  E <- .stored_edges(net)
  key <- paste(pmin(E$from, E$to), pmax(E$from, E$to))
  E <- E[!duplicated(key), , drop = FALSE]
  u <- E$from
  v <- E$to
  geom <- if ("geometry" %in% names(E)) E$geometry else vector("list", nrow(E))
  for (k in seq_along(geom)) {
    g <- geom[[k]]
    pu <- net$coords[u[k], ]
    if (is.null(g)) {
      geom[[k]] <- rbind(pu, net$coords[v[k], ], deparse.level = 0)
    } else {
      d0 <- sqrt(sum((g[1, ] - pu)^2))
      d1 <- sqrt(sum((g[nrow(g), ] - pu)^2))
      if (d0 > d1) geom[[k]] <- g[rev(seq_len(nrow(g))), , drop = FALSE]
    }
  }
  geom_len <- vapply(geom, function(g) sum(sqrt(rowSums(diff(g)^2))), numeric(1))
  length <- if ("length" %in% names(E)) .first_value(E$length) else rep(NA_real_, nrow(E))
  length[is.na(length)] <- geom_len[is.na(length)]
  cost <- if (isTRUE(net$metric$supports_edge_snapping)) {
    .edge_costs(net$metric, E, net)
  } else {
    rep(NA_real_, nrow(E))
  }
  out <- data.frame(u = u, v = v, length = length, cost = cost)
  out$geometry <- geom
  out$geom_length <- geom_len
  net$cache$edges <- out
  out
}

# Flattened segments of all street geometries, with their edge row.
.edge_segments <- function(net) {
  if (!is.null(net$cache$segments)) {
    return(net$cache$segments)
  }
  G <- edge_frame(net)$geometry
  nv <- vapply(G, nrow, integer(1))
  XY <- do.call(rbind, G)
  last <- cumsum(nv)
  i0 <- seq_len(nrow(XY))[-last]
  edge <- rep(seq_along(G), nv - 1L)
  x0 <- XY[i0, 1]
  y0 <- XY[i0, 2]
  x1 <- XY[i0 + 1L, 1]
  y1 <- XY[i0 + 1L, 2]
  len <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  cs <- cumsum(len)
  first <- c(1L, cumsum(nv - 1L) + 1L)[seq_along(G)]
  out <- data.frame(edge = edge, x0 = x0, y0 = y0, x1 = x1, y1 = y1, len = len,
                    cum = cs - len - c(0, cs)[first][edge])
  net$cache$segments <- out
  out
}

.snap <- function(u, v, cost_u, cost_v, depth) {
  list(u = as.integer(u), v = as.integer(v), cost_u = as.numeric(cost_u),
       cost_v = as.numeric(cost_v), depth = as.numeric(depth))
}

.snap_subset <- function(sn, i) {
  lapply(sn, function(a) a[i])
}

#' Attach points to a network
#'
#' Each point reaches node `u` at cost `cost_u` and node `v` at cost `cost_v`
#' along its street segment, after an off-network `depth` (Euclidean units).
#' Points snapped to a node have `u == v` and zero edge costs.
#'
#' @param net A [network_graph()].
#' @param xy Points: an `n x 2` matrix, a data frame with `x`/`y`, or sf points.
#' @param allowed Optional logical vector over nodes. With `mode = "edge"` a
#'   point may only snap to segments whose two endpoints are allowed.
#' @param mode `"edge"` (nearest street segment) or `"node"` (nearest node).
#'   Precomputed metrics always snap to nodes.
#' @return A list with integer vectors `u`, `v` (node indices) and numeric
#'   vectors `cost_u`, `cost_v`, `depth`.
#' @examples
#' nodes <- data.frame(name = 1:2, x = c(0, 100), y = 0)
#' net <- network_graph(data.frame(from = 1, to = 2, length = 100), nodes)
#' snap_points(net, rbind(c(30, 10)))
#' @export
snap_points <- function(net, xy, allowed = NULL, mode = "edge") {
  .check_network(net)
  xy <- .as_points(xy)
  allowed <- if (is.null(allowed)) rep(TRUE, net$n_nodes) else as.logical(allowed)
  if (mode == "edge" && !isTRUE(net$metric$supports_edge_snapping)) mode <- "node"
  if (mode == "node") {
    cand <- which(allowed)
    if (!length(cand)) stop("no network nodes available for snapping", call. = FALSE)
    cx <- net$coords[cand, 1]
    cy <- net$coords[cand, 2]
    r <- .nearest_segment(xy[, 1], xy[, 2], cx, cy, cx, cy)
    node <- cand[r$seg]
    z <- rep(0, nrow(xy))
    return(.snap(node, node, z, z, r$dist))
  }
  if (mode != "edge") stop("mode must be 'edge' or 'node'", call. = FALSE)
  E <- edge_frame(net)
  ok <- allowed[E$u] & allowed[E$v] & is.finite(E$cost)
  if (!any(ok)) stop("no network edges available for snapping", call. = FALSE)
  S <- .edge_segments(net)
  S <- S[ok[S$edge], , drop = FALSE]
  r <- .nearest_segment(xy[, 1], xy[, 2], S$x0, S$y0, S$x1, S$y1)
  ei <- S$edge[r$seg]
  glen <- E$geom_length[ei]
  frac <- ifelse(glen > 0, (S$cum[r$seg] + r$t * S$len[r$seg]) / glen, 0)
  cost <- E$cost[ei]
  .snap(E$u[ei], E$v[ei], frac * cost, (1 - frac) * cost, r$dist)
}

#' Split edges where they cross the boundary
#'
#' Inserts a node wherever an edge crosses the boundary line, so that every
#' path between the two sides passes through a node on the boundary. Columns
#' listed in `extensive` are divided in proportion to the split; all others
#' are copied. [network_coordinates()] does this automatically.
#'
#' @param net A [network_graph()].
#' @param b A [boundary()].
#' @param extensive Edge columns that scale with length.
#' @return A network (the same one when nothing needs splitting), with
#'   attribute `"n_crossings_inserted"`.
#' @examples
#' nodes <- data.frame(name = 1:2, x = 0, y = c(0, 100))
#' net <- network_graph(data.frame(from = 1, to = 2, length = 100), nodes)
#' b <- boundary(rbind(c(-10, 40), c(10, 40)), treated_region = rect_polygon(-50, -50, 50, 40))
#' split <- with_boundary_crossings(net, b)
#' split$edges
#' @export
with_boundary_crossings <- function(net, b, extensive = c("length", "travel_time")) {
  .check_network(net)
  .check_boundary(b)
  if (!isTRUE(net$metric$supports_edge_splitting)) {
    warning(sprintf("%s cannot split edges; make sure boundary nodes exist in the graph",
                    format(net$metric)), call. = FALSE)
    return(net)
  }
  if (!is.null(net$cache$split_boundary) && identical(net$cache$split_boundary, b)) {
    return(net$cache$split_net)
  }
  Tn <- boundary_is_treated(b, net$coords)
  on <- boundary_distance(b, net$coords) <= b$tol
  E <- net$edges
  iu <- E$from
  iv <- E$to
  cand <- which(!on[iu] & !on[iv] & (Tn[iu] != Tn[iv]))
  if (!length(cand)) {
    attr(net, "n_crossings_inserted") <- 0L
    net$cache$split_boundary <- b
    net$cache$split_net <- net
    return(net)
  }
  nodes <- net$nodes
  has_geom <- "geometry" %in% names(E)
  new_nodes <- list()
  crossing_ids <- character()
  new_edges <- list()
  drop <- integer()
  for (k in cand) {
    geom <- if (has_geom) E$geometry[[k]] else NULL
    if (is.null(geom)) geom <- rbind(net$coords[iu[k], ], net$coords[iv[k], ], deparse.level = 0)
    hit <- .first_crossing(geom, b$line)
    if (is.null(hit)) next # crosses the side partition away from this boundary line
    pt <- hit$point
    t <- hit$t
    key <- sprintf("%.6f %.6f", round(pt[1], 6), round(pt[2], 6))
    if (!key %in% names(crossing_ids)) {
      cid <- sprintf("rdnet_crossing_%d", length(crossing_ids))
      crossing_ids[key] <- cid
      new_nodes[[length(new_nodes) + 1]] <- data.frame(name = cid, x = pt[1], y = pt[2])
    }
    cid <- crossing_ids[[key]]
    first <- E[k, , drop = FALSE]
    second <- E[k, , drop = FALSE]
    for (a in intersect(extensive, names(E))) {
      val <- .first_value(E[[a]][k])
      first[[a]] <- val * t
      second[[a]] <- val * (1 - t)
    }
    if (has_geom && !is.null(E$geometry[[k]])) {
      first$geometry <- list(.line_substring(geom, 0, t))
      second$geometry <- list(.line_substring(geom, t, 1))
    }
    first$to <- -match(cid, crossing_ids)
    second$from <- -match(cid, crossing_ids)
    new_edges[[length(new_edges) + 1]] <- first
    new_edges[[length(new_edges) + 1]] <- second
    drop <- c(drop, k)
  }
  if (!length(drop)) {
    attr(net, "n_crossings_inserted") <- 0L
    net$cache$split_boundary <- b
    net$cache$split_net <- net
    return(net)
  }
  added <- do.call(rbind, new_nodes)
  ids <- nodes$name
  if (!is.character(ids)) ids <- as.character(ids)
  all_nodes <- nodes
  all_nodes$name <- ids
  extra <- as.data.frame(matrix(NA, nrow(added), ncol(nodes) - 3, dimnames = list(NULL, names(nodes)[-(1:3)])))
  added <- cbind(added, extra)
  all_nodes <- rbind(all_nodes, added[names(all_nodes)])
  all_nodes$rdnet_crossing <- c(rep(FALSE, nrow(nodes)), rep(TRUE, nrow(added)))
  NE <- do.call(rbind, new_edges)
  n0 <- nrow(nodes)
  NE$from <- ifelse(NE$from < 0, n0 - NE$from, NE$from)
  NE$to <- ifelse(NE$to < 0, n0 - NE$to, NE$to)
  E2 <- rbind(E[-drop, , drop = FALSE], NE)
  E2$from <- all_nodes$name[E2$from]
  E2$to <- all_nodes$name[E2$to]
  rownames(E2) <- NULL
  out <- network_graph(E2, all_nodes, metric = net$metric, directed = net$directed,
                       direction = net$direction, crs = net$crs, snap_scale = net$snap_scale)
  attr(out, "n_crossings_inserted") <- length(drop)
  net$cache$split_boundary <- b
  net$cache$split_net <- out
  out
}

#' Build a network from an adjacency matrix
#'
#' Each nonzero entry `A[i, j]` becomes an edge with `weight` (the entry) and
#' `length` (the straight-line distance between the two nodes).
#'
#' @param adjacency Square weighted adjacency matrix (base or `Matrix`).
#' @param coords `n x 2` node coordinates.
#' @param nodes Optional node ids (default `1:n`).
#' @param directed Whether entries are one-way.
#' @param metric Metric specification, as in [network_graph()]; defaults to the
#'   matrix entries.
#' @param ... Further arguments to [network_graph()].
#' @return A [network_graph()].
#' @examples
#' A <- rbind(c(0, 2, 0), c(2, 0, 3), c(0, 3, 0))
#' net <- network_graph_from_adjacency(A, rbind(c(0, 0), c(1, 0), c(2, 0)))
#' node_distance(net, 1, 3)
#' @export
network_graph_from_adjacency <- function(adjacency, coords, nodes = NULL, directed = FALSE,
                                         metric = "weight", ...) {
  if (inherits(adjacency, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("package 'Matrix' is needed", call. = FALSE)
    s <- Matrix::summary(methods::as(adjacency, "TsparseMatrix"))
    ij <- cbind(s$i, s$j)
    w <- s$x
  } else {
    A <- as.matrix(adjacency)
    ij <- which(A != 0, arr.ind = TRUE)
    ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
    w <- A[ij]
  }
  coords <- .as_points(coords)
  n <- nrow(coords)
  ids <- nodes %||% seq_len(n)
  keep <- ij[, 1] != ij[, 2]
  ij <- ij[keep, , drop = FALSE]
  w <- w[keep]
  if (!directed) {
    # later entries of an undirected pair overwrite earlier ones
    key <- paste(pmin(ij[, 1], ij[, 2]), pmax(ij[, 1], ij[, 2]))
    last <- !duplicated(key, fromLast = TRUE)
    w_last <- stats::setNames(w[last], key[last])
    sel <- which(!duplicated(key))
    w <- unname(w_last[key[sel]])
    ij <- ij[sel, , drop = FALSE]
  }
  len <- sqrt(rowSums((coords[ij[, 1], , drop = FALSE] - coords[ij[, 2], , drop = FALSE])^2))
  edges <- data.frame(from = ids[ij[, 1]], to = ids[ij[, 2]], weight = as.numeric(w), length = len)
  nodes_df <- data.frame(name = ids, x = coords[, 1], y = coords[, 2])
  network_graph(edges, nodes_df, metric = metric, directed = directed, ...)
}

#' Build a network from an igraph graph
#'
#' Vertices must carry projected coordinates; vertex names (if any) become
#' node ids and edge attributes become edge columns.
#'
#' @param g An `igraph` graph.
#' @param metric Metric specification, as in [network_graph()].
#' @param x,y Vertex attributes holding the coordinates.
#' @param ... Further arguments to [network_graph()].
#' @return A [network_graph()].
#' @examplesIf requireNamespace("igraph", quietly = TRUE)
#' g <- igraph::make_ring(4)
#' igraph::V(g)$x <- c(0, 1, 1, 0)
#' igraph::V(g)$y <- c(0, 0, 1, 1)
#' igraph::E(g)$weight <- 1
#' network_graph_from_igraph(g)
#' @export
network_graph_from_igraph <- function(g, metric = "weight", x = "x", y = "y", ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("package 'igraph' is needed", call. = FALSE)
  }
  df <- igraph::as_data_frame(g, what = "both")
  v <- df$vertices
  if (!"name" %in% names(v)) v <- cbind(name = seq_len(igraph::vcount(g)), v)
  network_graph(df$edges, v, metric = metric, directed = igraph::is_directed(g), x = x, y = y, ...)
}

#' Download an OpenStreetMap street network
#'
#' Downloads streets with \pkg{dodgr} (which uses \pkg{osmdata}) and projects
#' them to the local UTM zone with \pkg{sf}. If a package is missing or the
#' download fails and `offline_fallback` is `TRUE`, the bundled mock Midtown
#' network is returned with a warning, so scripts keep working offline.
#'
#' @param query Place name, e.g. `"Midtown Manhattan, New York"`.
#' @param bbox Alternatively a bounding box `c(xmin, ymin, xmax, ymax)` in
#'   longitude/latitude.
#' @param network_type `"walk"`, `"bike"` or `"drive"` (or a \pkg{dodgr}
#'   weighting profile).
#' @param metric Metric specification, as in [network_graph()].
#' @param offline_fallback Return the bundled mock network when the download fails.
#' @param dist,min_area `query` is geocoded with Nominatim and the bounding box
#'   of its top match is downloaded. When that match is a point (many
#'   neighbourhoods, such as Midtown Manhattan, are mapped as points) or covers
#'   less than `min_area` square metres (a single building), the streets
#'   within `dist` metres of it are downloaded instead, with a warning.
#' @param ... Further arguments to [network_graph()].
#' @return A [network_graph()] with `length` (m), `highway`, and
#'   `travel_time` (s) edge columns.
#' @examples
#' \donttest{
#' net <- suppressWarnings(network_graph_from_place("Midtown Manhattan, New York"))
#' net
#' }
#' @export
network_graph_from_place <- function(query = NULL, bbox = NULL, network_type = "walk",
                                     metric = "length", offline_fallback = TRUE, dist = 1000,
                                     min_area = 250000, ...) {
  if (is.null(query) && is.null(bbox)) stop("give a place query or a bbox", call. = FALSE)
  profile <- switch(network_type, walk = "foot", bike = "bicycle", drive = "motorcar", network_type)
  res <- tryCatch({
    where <- if (!is.null(bbox)) as.numeric(bbox) else .place_bbox(query, dist, min_area)
    .network_from_dodgr(.osm_download(where, profile), metric = metric, ...)
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (!offline_fallback) stop(res)
    warning(sprintf("OSM download failed (%s); using the bundled mock Midtown network",
                    conditionMessage(res)), call. = FALSE)
    return(load_mock_midtown()$network(metric, ...))
  }
  res
}

# Bounding box c(xmin, ymin, xmax, ymax) of the top geocoding match for `query`,
# or a square of half-width `dist` metres around it when the match is a point or
# covers less than `min_area` square metres.
.place_bbox <- function(query, dist, min_area) {
  bb <- .geocode_bbox(query)
  r <- .bbox_or_radius(bb, dist, min_area)
  if (r$small) {
    warning(sprintf("'%s' geocodes to an area of %s m2; using the streets within %s m of it", query,
                    format(round(r$area), big.mark = ","), format(dist)), call. = FALSE)
  }
  r$bbox
}

# Top Nominatim match as c(xmin, ymin, xmax, ymax) in longitude/latitude.
.geocode_bbox <- function(query) {
  if (!requireNamespace("osmdata", quietly = TRUE)) stop("package 'osmdata' is needed", call. = FALSE)
  m <- osmdata::getbb(query)
  bb <- c(m[1, 1], m[2, 1], m[1, 2], m[2, 2])
  if (anyNA(bb)) stop(sprintf("Nominatim found no match for '%s'", query), call. = FALSE)
  bb
}

.bbox_or_radius <- function(bb, dist, min_area) {
  lon0 <- (bb[1] + bb[3]) / 2
  lat0 <- (bb[2] + bb[4]) / 2
  mx <- 111320 * cos(lat0 * pi / 180) # metres per degree of longitude
  my <- 110574 # metres per degree of latitude
  area <- (bb[3] - bb[1]) * mx * (bb[4] - bb[2]) * my
  if (area >= min_area) {
    return(list(bbox = bb, small = FALSE, area = area))
  }
  list(bbox = c(lon0 - dist / mx, lat0 - dist / my, lon0 + dist / mx, lat0 + dist / my),
       small = TRUE, area = area)
}

# Street network as a dodgr edge table (needs dodgr, osmdata and internet access).
.osm_download <- function(bbox, profile) {
  for (p in c("dodgr", "sf")) {
    if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("package '%s' is needed", p), call. = FALSE)
  }
  streets <- dodgr::dodgr_streetnet(bbox, expand = 0, quiet = TRUE)
  dodgr::weight_streetnet(streets, wt_profile = profile)
}

.network_from_dodgr <- function(graph, metric = "length", ...) {
  graph <- as.data.frame(graph)
  lon <- c(graph$from_lon, graph$to_lon)
  lat <- c(graph$from_lat, graph$to_lat)
  ids <- c(as.character(graph$from_id), as.character(graph$to_id))
  keep <- !duplicated(ids)
  zone <- floor((mean(lon) + 180) / 6) + 1
  crs <- sprintf("EPSG:%d", if (mean(lat) >= 0) 32600 + zone else 32700 + zone)
  xy <- project_lonlat(lon[keep], lat[keep], crs)
  nodes <- data.frame(name = ids[keep], x = xy[, 1], y = xy[, 2])
  edges <- data.frame(from = as.character(graph$from_id), to = as.character(graph$to_id), length = graph$d)
  if ("highway" %in% names(graph)) edges$highway <- graph$highway
  if ("time" %in% names(graph)) edges$travel_time <- graph$time
  network_graph(edges, nodes, metric = metric, directed = TRUE, crs = crs, ...)
}

#' Project longitude/latitude to a planar coordinate system
#'
#' @param lon,lat Longitude and latitude (WGS84).
#' @param to_crs Target coordinate reference system (e.g. `"EPSG:32618"`).
#' @return An `n x 2` matrix of projected coordinates. Requires \pkg{sf}.
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' project_lonlat(-73.98, 40.76, "EPSG:32618")
#' @export
project_lonlat <- function(lon, lat, to_crs) {
  .need_sf()
  out <- sf::sf_project("EPSG:4326", to_crs, cbind(as.numeric(lon), as.numeric(lat)))
  unname(out[, 1:2, drop = FALSE])
}
