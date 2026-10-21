# Network metrics: how distance (dissimilarity) is measured on a network.
#
# A metric turns a graph into a shortest-path engine. Three families exist:
# edge costs (sum of per-edge costs), transition costs (edge costs plus a
# penalty for changing edge class, with transparent connector classes), and a
# precomputed node-by-node matrix. All costs must be non-negative; Inf marks an
# untraversable edge.

#' Network metrics
#'
#' A metric defines what distance means on a [network_graph()]. Distance is
#' dissimilarity: it does not have to be metres. The coordinates `s` and `o`,
#' the kernel distances and the bandwidths are all in the metric's units.
#'
#' * `edge_cost()`: path cost is the sum of edge costs. `cost` is an edge
#'   column name (`"length"`, `"travel_time"`) or a function. A function is
#'   called once, vectorised over all edges, as `cost(edges)` or, if it takes
#'   two arguments, `cost(edges, nodes)`; `edges` has columns `from` and `to`
#'   (node ids) plus every edge attribute, `nodes` has `name`, `x`, `y` and the
#'   node attributes. It returns one non-negative cost per edge (`Inf` blocks
#'   an edge).
#' * `transition_cost()`: edge costs plus a penalty for *changing* edge class
#'   (e.g. highway to local street). Edges whose class is in `connectors`
#'   (ramps, links) are transparent: the penalty is charged between the last
#'   non-connector class and the next one, so highway -> ramp -> highway stays
#'   cheap while highway -> ramp -> residential is expensive.
#' * `precomputed()`: a user-supplied node-by-node dissimilarity matrix.
#'   Observations then snap to nodes, since per-edge costs are unknown.
#' * `travel_time()`: seconds to traverse each edge from a length (m) and a
#'   speed (km/h) column.
#' * `node_contrast()`: `base + weight * |attr(u) - attr(v)|`, for land-use or
#'   other node-level dissimilarity.
#'
#' @param cost Edge column name or a vectorised function (see Details). For
#'   `precomputed_from_od()`: the column with the origin-destination cost.
#' @param snap_scale Cost per unit of off-network snapping distance. Defaults
#'   to 1 when the cost is `"length"` and 0 otherwise.
#' @param name Label used when printing.
#' @return An object of class `rd2dnet_metric`.
#' @examples
#' edge_cost("length")
#' edge_cost(function(edges) 2 * edges$length, name = "double length")
#' travel_time()
#' node_contrast("far", 50)
#' transition_cost("length", class_attr = "highway",
#'                 penalty = data.frame(from = "motorway", to = "residential", penalty = 2000),
#'                 connectors = c("motorway_link", "trunk_link"))
#' @name metrics
NULL

.new_metric <- function(x, subclass) {
  structure(x, class = c(subclass, "rd2dnet_metric"))
}

#' @rdname metrics
#' @export
edge_cost <- function(cost = "length", snap_scale = NULL, name = NULL) {
  if (!(is.character(cost) && length(cost) == 1) && !is.function(cost)) {
    stop("cost must be an edge column name or a function", call. = FALSE)
  }
  .new_metric(list(
    cost = cost,
    snap_scale = snap_scale,
    name = name %||% (if (is.character(cost)) cost else "custom"),
    supports_edge_splitting = TRUE,
    supports_edge_snapping = TRUE
  ), "rd2dnet_edge_cost")
}

#' @rdname metrics
#' @param edge_cost Base per-edge cost (column name or function), as in `edge_cost()`.
#' @param class_attr Edge column holding the edge class (OSM's `"highway"` by default).
#' @param penalty Transition penalties: a function `f(from_class, to_class)`
#'   returning a cost, a data frame with columns `from`, `to` and `penalty`,
#'   or a matrix with class names as row (from) and column (to) names. Pairs
#'   not given cost `default_penalty` when the classes differ and 0 otherwise.
#' @param default_penalty Penalty for unlisted pairs of different classes.
#' @param connectors Classes that do not reset the remembered class (ramps,
#'   links, crosswalks).
#' @param default_class Class of edges with a missing class value.
#' @export
transition_cost <- function(edge_cost = "length", class_attr = "highway", penalty = NULL,
                            default_penalty = 0, connectors = character(),
                            default_class = "unclassified", snap_scale = NULL) {
  base <- if (inherits(edge_cost, "rd2dnet_edge_cost")) edge_cost else rd2dnet::edge_cost(edge_cost, snap_scale = snap_scale)
  .new_metric(list(
    base = base,
    class_attr = class_attr,
    penalty = penalty,
    default_penalty = as.numeric(default_penalty),
    connectors = as.character(connectors),
    default_class = default_class,
    supports_edge_splitting = TRUE,
    supports_edge_snapping = TRUE
  ), "rd2dnet_transition_cost")
}

#' @rdname metrics
#' @param matrix Square matrix; `matrix[i, j]` is the distance from node
#'   `nodes[i]` to node `nodes[j]`.
#' @param nodes Node ids for the rows and columns. Defaults to the network's
#'   node order.
#' @export
precomputed <- function(matrix, nodes = NULL) {
  M <- as.matrix(matrix)
  storage.mode(M) <- "double"
  if (length(dim(M)) != 2 || nrow(M) != ncol(M)) {
    stop("precomputed metric must be a square matrix", call. = FALSE)
  }
  if (anyNA(M) || any(M < 0)) {
    stop("precomputed distances must be non-negative and not NA", call. = FALSE)
  }
  .new_metric(list(
    matrix = unname(M),
    nodes = nodes,
    supports_edge_splitting = FALSE,
    supports_edge_snapping = FALSE
  ), "rd2dnet_precomputed")
}

#' @rdname metrics
#' @param table For `precomputed_from_od()`: a long origin-destination table,
#'   one row per pair. Pairs not listed are unreachable (`Inf`), the diagonal
#'   is 0 and repeated pairs keep their smallest cost.
#' @param origin,destination Columns with the origin and destination node ids.
#' @param symmetric Whether a listed pair `a -> b` also gives `b -> a` when
#'   that pair is not listed.
#' @export
precomputed_from_od <- function(table, origin = "origin", destination = "destination", cost = "cost",
                                nodes = NULL, symmetric = FALSE) {
  table <- as.data.frame(table)
  o <- table[[origin]]
  d <- table[[destination]]
  c <- as.numeric(table[[cost]])
  if (anyNA(c) || any(c < 0)) stop("OD costs must be non-negative and not missing", call. = FALSE)
  ids <- unique(c(o, d))
  if (is.null(nodes)) nodes <- sort(ids)
  if (!all(ids %in% nodes)) stop("the OD table has ids that are not in nodes", call. = FALSE)
  k <- length(nodes)
  M <- matrix(Inf, k, k)
  i <- match(o, nodes)
  j <- match(d, nodes)
  for (r in seq_along(c)) M[i[r], j[r]] <- min(M[i[r], j[r]], c[r])
  if (symmetric) {
    fill <- !is.finite(t(M)) & is.finite(M)
    M[t(fill)] <- t(M)[t(fill)]
  }
  diag(M) <- 0
  precomputed(M, nodes)
}

#' @rdname metrics
#' @param speed_attr,length_attr Edge columns with the speed (km/h) and length (m).
#' @param default_kph Speed used where `speed_attr` is absent or `NA`. For
#'   data in feet and miles per hour, convert first (`length = 0.3048 *
#'   length_ft`, `speed_kph = 1.609344 * speed_mph`).
#' @export
travel_time <- function(speed_attr = "speed_kph", length_attr = "length", default_kph = NULL) {
  warned <- FALSE
  fn <- function(edges) {
    speed <- if (speed_attr %in% names(edges)) .first_value(edges[[speed_attr]]) else rep(NA_real_, nrow(edges))
    if (!is.null(default_kph)) {
      if (!warned && nrow(edges) && mean(is.na(speed)) > 0.5) {
        warning(sprintf(paste0("%.0f%% of edges have no '%s' and use default_kph = %g: travel times are then ",
                               "mostly length / one speed. OpenStreetMap downloads carry speeds only where ",
                               "mapped"), 100 * mean(is.na(speed)), speed_attr, default_kph), call. = FALSE)
        warned <<- TRUE
      }
      speed[is.na(speed)] <- default_kph
    }
    miss <- which(is.na(speed))
    if (length(miss)) {
      stop(sprintf("edge (%s, %s) has no '%s' and no default_kph", edges$from[miss[1]],
                   edges$to[miss[1]], speed_attr), call. = FALSE)
    }
    len <- .first_value(edges[[length_attr]])
    ifelse(speed <= 0, Inf, len / (speed / 3.6))
  }
  edge_cost(fn, name = "travel_time")
}

#' @rdname metrics
#' @param attr Node column compared across each edge.
#' @param weight Cost per unit of absolute difference in `attr`.
#' @param base Edge column with the base cost.
#' @export
node_contrast <- function(attr, weight, base = "length", snap_scale = NULL) {
  force(attr)
  force(weight)
  force(base)
  fn <- function(edges, nodes) {
    a <- if (attr %in% names(nodes)) as.numeric(nodes[[attr]]) else rep(0, nrow(nodes))
    a[is.na(a)] <- 0
    i <- match(edges$from, nodes$name)
    j <- match(edges$to, nodes$name)
    .first_value(edges[[base]]) + weight * abs(a[i] - a[j])
  }
  edge_cost(fn, snap_scale = snap_scale %||% (if (identical(base, "length")) 1 else NULL),
            name = sprintf("%s+%s*|d%s|", base, format(weight), attr))
}

#' @export
format.rd2dnet_metric <- function(x, ...) {
  if (inherits(x, "rd2dnet_edge_cost")) {
    sprintf("edge_cost(\"%s\")", x$name)
  } else if (inherits(x, "rd2dnet_transition_cost")) {
    sprintf("transition_cost(\"%s\", class_attr = \"%s\")", x$base$name, x$class_attr)
  } else if (inherits(x, "rd2dnet_precomputed")) {
    sprintf("precomputed(%d x %d)", nrow(x$matrix), ncol(x$matrix))
  } else {
    "rd2dnet_metric"
  }
}

#' @export
print.rd2dnet_metric <- function(x, ...) {
  cat("<rd2dnet metric>", format(x), "\n")
  invisible(x)
}

# Coerce a metric specification (name, function, matrix or metric) to a metric.
as_metric <- function(spec) {
  if (inherits(spec, "rd2dnet_metric")) {
    return(spec)
  }
  if ((is.character(spec) && length(spec) == 1) || is.function(spec)) {
    return(edge_cost(spec))
  }
  if (is.matrix(spec)) {
    return(precomputed(spec))
  }
  stop("cannot interpret this object as a network metric", call. = FALSE)
}

# OSM attributes are sometimes list-valued (merged ways): use the first entry.
.first_value <- function(x) {
  if (is.list(x)) {
    x <- vapply(x, function(v) if (length(v)) as.numeric(v[[1]]) else NA_real_, numeric(1))
  }
  as.numeric(x)
}

.default_snap_scale <- function(metric) {
  if (inherits(metric, "rd2dnet_transition_cost")) {
    return(.default_snap_scale(metric$base))
  }
  if (inherits(metric, "rd2dnet_edge_cost")) {
    if (!is.null(metric$snap_scale)) {
      return(as.numeric(metric$snap_scale))
    }
    # NA: the snapped edge's cost per unit length (see .snap_rates())
    return(if (identical(metric$cost, "length")) 1 else NA_real_)
  }
  0
}

# Per-edge costs for the edge table `E` (columns from/to are node indices).
.edge_costs <- function(metric, E, net) {
  if (inherits(metric, "rd2dnet_transition_cost")) {
    return(.edge_costs(metric$base, E, net))
  }
  if (!inherits(metric, "rd2dnet_edge_cost")) {
    return(rep(NA_real_, nrow(E)))
  }
  labels <- net$nodes$name
  if (is.character(metric$cost)) {
    attr <- metric$cost
    if (!attr %in% names(E)) {
      stop(sprintf("edges have no attribute '%s'", attr), call. = FALSE)
    }
    out <- .first_value(E[[attr]])
    miss <- which(is.na(out))
    if (length(miss)) {
      k <- miss[1]
      stop(sprintf("edge (%s, %s) has no value for attribute '%s'", labels[E$from[k]],
                   labels[E$to[k]], attr), call. = FALSE)
    }
    return(out)
  }
  ed <- E
  ed$from <- labels[E$from]
  ed$to <- labels[E$to]
  fn <- metric$cost
  out <- if (length(formals(fn)) >= 2) fn(ed, net$nodes) else fn(ed)
  out <- as.numeric(out)
  if (length(out) != nrow(E)) {
    stop("a cost function must return one value per edge", call. = FALSE)
  }
  out
}

.check_costs <- function(c, E, labels) {
  bad <- which(is.na(c) | c < 0)
  if (length(bad)) {
    k <- bad[1]
    stop(sprintf("metric produced invalid cost %s for edge (%s, %s); costs must be >= 0",
                 format(c[k]), labels[E$from[k]], labels[E$to[k]]), call. = FALSE)
  }
  c
}

.check_cost <- function(c, where) {
  c <- as.numeric(c)
  if (length(c) != 1 || is.na(c) || c < 0) {
    stop(sprintf("metric produced invalid cost %s for %s; costs must be >= 0", format(c), where),
         call. = FALSE)
  }
  c
}

# Stored edges -> traversable arcs, list(tail, head, idx) with idx the edge row.
.orient <- function(net, from, to) {
  idx <- seq_along(from)
  if (!net$directed || net$direction == "both") {
    return(list(tail = c(from, to), head = c(to, from), idx = c(idx, idx)))
  }
  if (net$direction == "out") {
    return(list(tail = from, head = to, idx = idx))
  }
  list(tail = to, head = from, idx = idx)
}

# Compressed sparse row graph (0-based) from 1-based arcs.
.csr <- function(tail, head, w, n) {
  o <- order(tail, method = "radix")
  list(indptr = as.integer(c(0L, cumsum(tabulate(tail, n)))),
       indices = as.integer(head[o] - 1L), data = as.numeric(w[o]), n = as.integer(n))
}

.new_engine <- function(csr, n_nodes, start_state, state_node, identity) {
  structure(list(type = "csr", csr = csr, n_nodes = as.integer(n_nodes),
                 start_state = as.integer(start_state), state_node = as.integer(state_node),
                 identity = identity), class = "rd2dnet_engine")
}

# Build the shortest-path engine for a network (called once, lazily).
.compile_metric <- function(metric, net) {
  n <- nrow(net$nodes)
  E <- .stored_edges(net)
  labels <- net$nodes$name
  if (inherits(metric, "rd2dnet_precomputed")) {
    M <- metric$matrix
    if (nrow(M) != n) {
      stop(sprintf("precomputed matrix is %dx%d but the graph has %d nodes", nrow(M), nrow(M), n),
           call. = FALSE)
    }
    if (!is.null(metric$nodes)) {
      perm <- match(labels, metric$nodes)
      if (anyNA(perm)) stop("precomputed `nodes` do not cover the network's nodes", call. = FALSE)
      M <- M[perm, perm, drop = FALSE]
    }
    return(structure(list(type = "dense", D = M, n_nodes = n), class = "rd2dnet_engine"))
  }
  cost <- .check_costs(.edge_costs(metric, E, net), E, labels)
  fin <- is.finite(cost)
  arcs <- .orient(net, E$from[fin], E$to[fin])
  acost <- cost[fin][arcs$idx]
  if (inherits(metric, "rd2dnet_edge_cost")) {
    idx <- seq_len(n) - 1L
    return(.new_engine(.csr(arcs$tail, arcs$head, acost, n), n, idx, idx, TRUE))
  }
  # transition costs: state-expanded graph
  cls_raw <- if (metric$class_attr %in% names(E)) E[[metric$class_attr]] else rep(NA, nrow(E))
  if (is.list(cls_raw)) {
    cls_raw <- vapply(cls_raw, function(v) if (length(v)) as.character(v[[1]]) else NA_character_, "")
  }
  cls_raw <- as.character(cls_raw)
  cls_raw[is.na(cls_raw)] <- as.character(metric$default_class)
  ecls <- cls_raw[fin][arcs$idx]
  tail <- arcs$tail
  head <- arcs$head
  m <- length(tail)
  is_conn <- ecls %in% metric$connectors
  classes <- sort(unique(ecls[!is_conn]))
  n_mem <- length(classes) + 1L # memory 1 = no class yet
  pen <- matrix(0, n_mem, n_mem)
  for (a in seq_along(classes)) {
    for (b in seq_along(classes)) {
      pen[a + 1L, b + 1L] <- .transition_penalty(metric, classes[a], classes[b])
    }
  }
  if (net$directed && net$direction == "in") {
    # arcs are reversed, so the search meets each turn from its far side
    pen <- t(pen)
  }
  cid <- match(ecls, classes) + 1L
  cid[is.na(cid) | is_conn] <- 1L
  # States: one "standing at node" state per node, one per non-connector arc,
  # and one per (connector arc, remembered class).
  nc <- which(!is_conn)
  cn <- which(is_conn)
  nc_state <- rep(NA_integer_, m)
  nc_state[nc] <- n + seq_along(nc)
  conn_base <- rep(NA_integer_, m)
  conn_base[cn] <- n + length(nc) + (seq_along(cn) - 1L) * n_mem
  n_states <- n + length(nc) + length(cn) * n_mem
  state_node <- c(seq_len(n), head[nc], rep(head[cn], each = n_mem))
  # Every state (at a node, with a remembered class) leaves along each out-arc of its node.
  from_state <- c(seq_len(n), nc_state[nc], rep(conn_base[cn], each = n_mem) + rep(seq_len(n_mem), times = length(cn)))
  from_node <- state_node[from_state]
  from_mem <- c(rep(1L, n), cid[nc], rep(seq_len(n_mem), times = length(cn)))
  ord <- order(tail, method = "radix")
  indptr <- c(0L, cumsum(tabulate(tail, n)))
  deg <- indptr[from_node + 1L] - indptr[from_node]
  rows <- rep(from_state, deg)
  mem <- rep(from_mem, deg)
  start <- rep(indptr[from_node], deg)
  within <- seq_len(sum(deg)) - 1L - rep(cumsum(deg) - deg, deg)
  k <- ord[start + within + 1L]
  conn_k <- is_conn[k]
  cols <- ifelse(conn_k, conn_base[k] + mem, nc_state[k])
  vals <- acost[k] + ifelse(conn_k, 0, pen[cbind(mem, cid[k])])
  .new_engine(.csr(rows, cols, vals, n_states), n, seq_len(n) - 1L, state_node - 1L, FALSE)
}

.transition_penalty <- function(metric, a, b) {
  pen <- metric$penalty
  where <- sprintf("transition %s -> %s", a, b)
  if (is.function(pen)) {
    return(.check_cost(pen(a, b), where))
  }
  if (is.data.frame(pen)) {
    hit <- which(as.character(pen$from) == a & as.character(pen$to) == b)
    if (length(hit)) return(.check_cost(pen$penalty[hit[1]], where))
  } else if (is.matrix(pen)) {
    if (a %in% rownames(pen) && b %in% colnames(pen) && !is.na(pen[a, b])) {
      return(.check_cost(pen[a, b], where))
    }
  } else if (!is.null(pen)) {
    stop("penalty must be a function, a data frame (from, to, penalty) or a matrix", call. = FALSE)
  }
  .check_cost(if (a == b) 0 else metric$default_penalty, where)
}

# --- engine queries (node indices are 1-based) ---------------------------------

.engine_distances <- function(eng, sources, limit = Inf) {
  sources <- as.integer(sources)
  if (eng$type == "flow") {
    return(pmin(.engine_distances(eng$out, sources, limit), .engine_distances(eng$inn, sources, limit)))
  }
  if (eng$type == "dense") {
    out <- eng$D[sources, , drop = FALSE]
    out[out > limit] <- Inf
    return(out)
  }
  g <- eng$csr
  .dijkstra_distances(g$indptr, g$indices, g$data, g$n, eng$start_state[sources], eng$state_node,
                      eng$n_nodes, as.numeric(limit))
}

# For every node: distance to the nearest source and that source's node index (NA if none).
.engine_nearest <- function(eng, sources) {
  sources <- as.integer(sources)
  if (eng$type == "flow") {
    a <- .engine_nearest(eng$out, sources)
    b <- .engine_nearest(eng$inn, sources)
    take_b <- b$dist < a$dist
    return(list(dist = ifelse(take_b, b$dist, a$dist), src = ifelse(take_b, b$src, a$src)))
  }
  if (eng$type == "dense") {
    sub <- eng$D[sources, , drop = FALSE]
    arg <- apply(sub, 2, which.min)
    dist <- sub[cbind(arg, seq_len(ncol(sub)))]
    return(list(dist = dist, src = ifelse(is.finite(dist), sources[arg], NA_integer_)))
  }
  g <- eng$csr
  r <- .dijkstra_nearest(g$indptr, g$indices, g$data, g$n, eng$start_state[sources], eng$state_node,
                         eng$n_nodes)
  list(dist = r$dist, src = ifelse(r$src > 0, sources[pmax(r$src, 1L)], NA_integer_))
}

# (i, node, dist) for every node closer than `limit` to sources[i].
.engine_pairs <- function(eng, sources, limit) {
  sources <- as.integer(sources)
  if (eng$type == "flow") {
    a <- .engine_pairs(eng$out, sources, limit)
    b <- .engine_pairs(eng$inn, sources, limit)
    i <- c(a$i, b$i)
    node <- c(a$node, b$node)
    dist <- c(a$dist, b$dist)
    o <- order(i, node, dist)
    keep <- o[!duplicated(cbind(i, node)[o, , drop = FALSE])]
    keep <- keep[order(i[keep], node[keep])]
    return(list(i = i[keep], node = node[keep], dist = dist[keep]))
  }
  if (eng$type == "dense") {
    D <- eng$D[sources, , drop = FALSE]
    w <- which(D < limit, arr.ind = TRUE)
    w <- w[order(w[, 1], w[, 2]), , drop = FALSE]
    return(list(i = unname(w[, 1]), node = unname(w[, 2]), dist = D[w]))
  }
  g <- eng$csr
  .dijkstra_pairs(g$indptr, g$indices, g$data, g$n, eng$start_state[sources], eng$state_node,
                  eng$n_nodes, as.numeric(limit))
}
