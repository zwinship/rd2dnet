# Internal planar geometry: polylines, polygons and conversions from sf.
#
# A polyline is an n x 2 numeric matrix of vertices. A polygon is a list of
# ring matrices evaluated with the even-odd rule, so holes and multipolygons
# need no special handling. Geometries used as masks are tagged lists
# `list(type = "polygon" | "line" | "point", parts = list(<matrix>, ...))`.

.is_sf <- function(x) inherits(x, c("sf", "sfc", "sfg"))

.need_sf <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("package 'sf' is needed to use sf geometries", call. = FALSE)
  }
}

# sf geometry -> tagged list with one matrix per ring / line part / point.
.sf_parts <- function(x) {
  .need_sf()
  type <- unique(as.character(sf::st_geometry_type(x)))
  cc <- sf::st_coordinates(x)
  xy <- unname(cc[, c("X", "Y"), drop = FALSE])
  lcols <- grep("^L[0-9]+$", colnames(cc), value = TRUE)
  kind <- if (all(type %in% c("POLYGON", "MULTIPOLYGON"))) {
    "polygon"
  } else if (all(type %in% c("LINESTRING", "MULTILINESTRING"))) {
    "line"
  } else if (all(type %in% c("POINT", "MULTIPOINT"))) {
    "point"
  } else {
    stop(sprintf("unsupported geometry type: %s", paste(type, collapse = ", ")), call. = FALSE)
  }
  if (kind == "point" || length(lcols) == 0) {
    return(list(type = kind, parts = list(xy)))
  }
  key <- do.call(paste, c(unname(as.list(as.data.frame(cc[, lcols, drop = FALSE]))), sep = "/"))
  parts <- unname(split.data.frame(xy, factor(key, levels = unique(key))))
  list(type = kind, parts = parts)
}

# Observation points -> n x 2 matrix.
.as_points <- function(points) {
  if (.is_sf(points)) {
    return(.sf_parts(points)$parts[[1]])
  }
  if (is.data.frame(points)) {
    if (!all(c("x", "y") %in% names(points))) {
      stop("a data frame of points needs `x` and `y` columns", call. = FALSE)
    }
    return(cbind(as.numeric(points$x), as.numeric(points$y)))
  }
  if (is.numeric(points) && is.null(dim(points)) && length(points) == 2) {
    return(matrix(as.numeric(points), 1, 2))
  }
  m <- as.matrix(points)
  if (!is.numeric(m) || length(dim(m)) != 2 || ncol(m) != 2) {
    stop("points must be an (n x 2) matrix, a data frame with x/y columns, or sf points", call. = FALSE)
  }
  storage.mode(m) <- "double"
  unname(m)
}

# A boundary line -> polyline matrix.
.as_line <- function(line) {
  if (.is_sf(line)) {
    g <- .sf_parts(line)
    if (g$type != "line" || length(g$parts) != 1) {
      stop("the boundary must be a single LINESTRING", call. = FALSE)
    }
    return(g$parts[[1]])
  }
  m <- .as_points(line)
  if (nrow(m) < 2) stop("a line needs at least two vertices", call. = FALSE)
  m
}

# A region or mask geometry -> tagged list. Plain matrices are polygon rings.
.as_geom <- function(g) {
  if (.is_sf(g)) {
    return(.sf_parts(g))
  }
  if (is.list(g) && !is.data.frame(g) && !is.null(g$type) && !is.null(g$parts)) {
    return(g)
  }
  if (is.list(g) && !is.data.frame(g)) {
    return(list(type = "polygon", parts = lapply(g, .as_points)))
  }
  list(type = "polygon", parts = list(.as_points(g)))
}

#' Rectangle polygon
#'
#' A convenience constructor for an axis-aligned rectangle, usable wherever a
#' polygon is accepted (treated regions, masks).
#'
#' @param xmin,ymin,xmax,ymax Corners of the rectangle.
#' @return A 5 x 2 matrix of ring vertices (closed).
#' @examples
#' rect_polygon(0, 0, 10, 5)
#' @export
rect_polygon <- function(xmin, ymin, xmax, ymax) {
  matrix(c(xmin, xmax, xmax, xmin, xmin, ymin, ymin, ymax, ymax, ymin), ncol = 2)
}

# Segments of a polyline with their lengths and starting arc positions.
.segments <- function(L) {
  n <- nrow(L)
  x0 <- L[-n, 1]
  y0 <- L[-n, 2]
  x1 <- L[-1, 1]
  y1 <- L[-1, 2]
  len <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1, len = len, cum = c(0, cumsum(len))[seq_along(len)],
       total = sum(len))
}

# Segments of every part of a tagged geometry (rings are closed).
.geom_segments <- function(g) {
  segs <- lapply(g$parts, function(m) {
    if (g$type == "polygon" && (m[1, 1] != m[nrow(m), 1] || m[1, 2] != m[nrow(m), 2])) {
      m <- rbind(m, m[1, ])
    }
    if (g$type == "point" || nrow(m) == 1) {
      return(list(x0 = m[, 1], y0 = m[, 2], x1 = m[, 1], y1 = m[, 2]))
    }
    s <- .segments(m)
    s[c("x0", "y0", "x1", "y1")]
  })
  list(x0 = unlist(lapply(segs, `[[`, "x0")), y0 = unlist(lapply(segs, `[[`, "y0")),
       x1 = unlist(lapply(segs, `[[`, "x1")), y1 = unlist(lapply(segs, `[[`, "y1")))
}

.nearest <- function(xy, s) {
  .nearest_segment(xy[, 1], xy[, 2], s$x0, s$y0, s$x1, s$y1)
}

# Arc-length position of each point's projection onto a polyline.
.line_locate <- function(segs, xy) {
  r <- .nearest(xy, segs)
  segs$cum[r$seg] + r$t * segs$len[r$seg]
}

.line_distance <- function(segs, xy) .nearest(xy, segs)$dist

# Points at the given arc positions (clamped to the line).
.line_interpolate <- function(segs, arc) {
  a <- pmin(pmax(as.numeric(arc), 0), segs$total)
  k <- pmax(findInterval(a, segs$cum), 1L)
  t <- ifelse(segs$len[k] > 0, pmin(pmax((a - segs$cum[k]) / segs$len[k], 0), 1), 0)
  cbind(segs$x0[k] + t * (segs$x1[k] - segs$x0[k]), segs$y0[k] + t * (segs$y1[k] - segs$y0[k]))
}

# Part of a polyline between normalized positions t0 < t1.
.line_substring <- function(L, t0, t1) {
  s <- .segments(L)
  a0 <- t0 * s$total
  a1 <- t1 * s$total
  inner <- which(c(s$cum, s$total) > a0 & c(s$cum, s$total) < a1)
  rbind(.line_interpolate(s, a0), L[inner, , drop = FALSE], .line_interpolate(s, a1))
}

.geom_contains <- function(g, xy) {
  if (g$type != "polygon") {
    return(rep(FALSE, nrow(xy)))
  }
  rings <- g$parts
  .points_in_rings(xy[, 1], xy[, 2], unlist(lapply(rings, function(m) m[, 1])),
                   unlist(lapply(rings, function(m) m[, 2])),
                   rep(seq_along(rings), vapply(rings, nrow, integer(1))))
}

.geom_distance <- function(g, xy) {
  d <- .nearest(xy, .geom_segments(g))$dist
  if (g$type == "polygon") d[.geom_contains(g, xy)] <- 0
  d
}

.cross2 <- function(ax, ay, bx, by) ax * by - ay * bx

# Where polyline G first meets polyline B, as list(point, t) with t the
# normalized position along G; NULL when they do not meet.
.first_crossing <- function(G, B) {
  sg <- .segments(G)
  sb <- .segments(B)
  ng <- length(sg$len)
  nb <- length(sb$len)
  i <- rep(seq_len(ng), times = nb)
  j <- rep(seq_len(nb), each = ng)
  px <- sg$x0[i]
  py <- sg$y0[i]
  rx <- sg$x1[i] - px
  ry <- sg$y1[i] - py
  qx <- sb$x0[j]
  qy <- sb$y0[j]
  sx <- sb$x1[j] - qx
  sy <- sb$y1[j] - qy
  den <- .cross2(rx, ry, sx, sy)
  t <- .cross2(qx - px, qy - py, sx, sy) / den
  u <- .cross2(qx - px, qy - py, rx, ry) / den
  hit <- den != 0 & t >= 0 & t <= 1 & u >= 0 & u <= 1
  arc <- sg$cum[i] + t * sg$len[i]
  if (!any(hit)) {
    # collinear overlaps: take the start of the overlap along G
    rr <- rx * rx + ry * ry
    col <- den == 0 & .cross2(qx - px, qy - py, rx, ry) == 0 & rr > 0
    if (!any(col)) {
      return(NULL)
    }
    ta <- ((qx - px) * rx + (qy - py) * ry) / rr
    tb <- ta + (sx * rx + sy * ry) / rr
    lo <- pmax(0, pmin(ta, tb))
    hi <- pmin(1, pmax(ta, tb))
    col <- col & lo <= hi
    if (!any(col)) {
      return(NULL)
    }
    t <- lo
    arc <- sg$cum[i] + t * sg$len[i]
    hit <- col
  }
  k <- which(hit)[which.min(arc[hit])]
  point <- c(px[k] + t[k] * rx[k], py[k] + t[k] * ry[k])
  list(point = point, t = if (sg$total > 0) arc[k] / sg$total else 0.5)
}
