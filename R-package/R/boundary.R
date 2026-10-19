#' Boundary separating treated from control locations
#'
#' A boundary line, the side that is treated, and optional masked stretches.
#'
#' @param line The boundary: an `n x 2` matrix (or data frame with `x`/`y`) of
#'   vertices, or an sf `LINESTRING`, in the same projected coordinates as the
#'   network.
#' @param treated_region Polygon containing the treated side: a ring matrix, a
#'   list of ring matrices (holes and several parts use the even-odd rule), or
#'   an sf polygon. Points inside it, or within `tol` of it, are treated.
#' @param treated_side Alternatively `"left"` or `"right"` of the line's direction.
#' @param mask A list of geometries (polygons as above, or sf geometries such as
#'   a park) and/or numeric `c(start, end)` arc-length intervals along the line.
#'   Evaluation points touching them are skipped: use this for stretches with
#'   no mass on one side (Central Park on 60th Street).
#' @param tol Distance (in the network's coordinate units) at which a network
#'   node counts as lying on the boundary. Increase it when the boundary
#'   follows a street but its vertices are not exactly on the street's nodes.
#'
#' @details Points exactly on the line count as treated with `treated_region`
#'   (and are assigned by `treated_side` otherwise); pass `treated` to
#'   [network_coordinates()] to assign them yourself. Everything outside
#'   `treated_region` is control, so restrict the sample to the two
#'   neighbouring areas first. A closed line (a zone outline, see
#'   [boundary_from_region()]) is treated as a ring: offsets along it wrap
#'   around the shorter way.
#'
#' @return An object of class `rd2dnet_boundary`.
#' @seealso [boundary_project()], [boundary_is_treated()], [boundary_is_masked()].
#' @examples
#' b <- boundary(rbind(c(0, 0), c(100, 0)), treated_side = "left",
#'               mask = list(rect_polygon(20, 0, 40, 10), c(80, 90)))
#' b
#' boundary_is_treated(b, rbind(c(50, 1), c(50, -1)))
#' boundary_is_masked(b, c(10, 30, 50, 85))
#' @export
boundary <- function(line, treated_region = NULL, treated_side = NULL, mask = NULL, tol = 1e-6) {
  gtype <- if (.is_sf(line)) as.character(sf::st_geometry_type(line)) else character()
  if (length(gtype) && all(gtype %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("line is a polygon: use boundary_from_region(polygon) to make its outline the boundary with ",
         "the inside treated", call. = FALSE)
  }
  if (any(gtype == "MULTILINESTRING")) {
    L <- .join_lines(.sf_parts(line)$parts)
    if (is.null(L)) stop("the boundary has several disconnected pieces; pass one connected line", call. = FALSE)
  } else {
    L <- .as_line(line)
  }
  seglen <- sqrt(rowSums(diff(L)^2))
  if (!all(is.finite(L)) || !(sum(seglen) > 0)) {
    stop("the boundary line has zero length (its vertices coincide) or non-finite vertices", call. = FALSE)
  }
  if (is.null(treated_region) == is.null(treated_side)) {
    stop("give exactly one of treated_region or treated_side", call. = FALSE)
  }
  if (!is.null(treated_side) && !(identical(treated_side, "left") || identical(treated_side, "right"))) {
    stop("treated_side must be 'left' or 'right'", call. = FALSE)
  }
  if (.is_sf(mask) || is.matrix(mask) || (is.numeric(mask) && is.null(dim(mask)))) mask <- list(mask)
  mask_geoms <- list()
  mask_arcs <- list()
  for (m in mask) {
    if (is.numeric(m) && is.null(dim(m)) && length(m) == 2) {
      mask_arcs[[length(mask_arcs) + 1]] <- c(min(m), max(m))
    } else {
      mask_geoms[[length(mask_geoms) + 1]] <- .as_geom(m)
    }
  }
  structure(
    list(
      line = L,
      segs = .segments(L),
      treated_region = if (is.null(treated_region)) NULL else .as_geom(treated_region),
      treated_side = treated_side,
      mask_geoms = mask_geoms,
      mask_arcs = mask_arcs,
      tol = as.numeric(tol),
      closed = nrow(L) > 2 && all(L[1, ] == L[nrow(L), ]),
      length = .segments(L)$total
    ),
    class = "rd2dnet_boundary"
  )
}

#' @export
print.rd2dnet_boundary <- function(x, ...) {
  side <- if (!is.null(x$treated_region)) "region" else x$treated_side
  cat(sprintf("<rd2dnet boundary> length=%.1f, treated=%s, masks=%d\n", x$length, side,
              length(x$mask_geoms) + length(x$mask_arcs)))
  invisible(x)
}

.check_boundary <- function(b) {
  if (!inherits(b, "rd2dnet_boundary")) {
    stop("expected a boundary created with boundary()", call. = FALSE)
  }
  b
}

#' Boundary geometry operations
#'
#' Distances to the boundary line, arc-length positions along it, and points
#' at given positions.
#'
#' @param b A [boundary()].
#' @param xy Points: an `n x 2` matrix, a data frame with `x`/`y`, or sf points.
#' @param arc Arc-length positions along the boundary line.
#' @return `boundary_distance()`: Euclidean distance from each point to the
#'   line. `boundary_project()`: arc-length position of each point's
#'   projection onto the line. `boundary_interpolate()`: an `n x 2` matrix of
#'   coordinates.
#' @examples
#' b <- boundary(rbind(c(0, 0), c(100, 0)), treated_side = "left")
#' boundary_distance(b, rbind(c(10, 5), c(120, 0)))
#' boundary_project(b, rbind(c(10, 5), c(120, 0)))
#' boundary_interpolate(b, c(0, 50))
#' @name boundary_geometry
NULL

#' @rdname boundary_geometry
#' @export
boundary_distance <- function(b, xy) {
  .line_distance(.check_boundary(b)$segs, .as_points(xy))
}

#' @rdname boundary_geometry
#' @export
boundary_project <- function(b, xy) {
  .line_locate(.check_boundary(b)$segs, .as_points(xy))
}

#' @rdname boundary_geometry
#' @export
boundary_interpolate <- function(b, arc) {
  .line_interpolate(.check_boundary(b)$segs, arc)
}

#' Treatment side of points
#'
#' @inheritParams boundary_geometry
#' @return Logical: whether each point lies on the treated side (inside
#'   `treated_region`, within `tol` of it, or on `treated_side` of the line).
#' @examples
#' b <- boundary(rbind(c(0, 0), c(10, 0)), treated_side = "right")
#' boundary_is_treated(b, rbind(c(5, 1), c(5, -1)))
#' @export
boundary_is_treated <- function(b, xy) {
  .check_boundary(b)
  xy <- .as_points(xy)
  if (!is.null(b$treated_region)) {
    g <- b$treated_region
    return(.geom_contains(g, xy) | .geom_distance(g, xy) <= b$tol)
  }
  arc <- .line_locate(b$segs, xy)
  eps <- min(1e-3, b$length * 1e-6)
  a <- .line_interpolate(b$segs, pmin(pmax(arc - eps, 0), b$length))
  z <- .line_interpolate(b$segs, pmin(pmax(arc + eps, 0), b$length))
  foot <- .line_interpolate(b$segs, arc)
  tx <- z[, 1] - a[, 1]
  ty <- z[, 2] - a[, 2]
  left <- tx * (xy[, 2] - foot[, 2]) - ty * (xy[, 1] - foot[, 1]) > 0
  if (identical(b$treated_side, "left")) left else !left
}

#' Masked stretches of the boundary
#'
#' @inheritParams boundary_geometry
#' @return Logical: whether each arc-length position falls in a masked stretch.
#' @examples
#' b <- boundary(rbind(c(0, 0), c(100, 0)), treated_side = "left", mask = list(c(80, 90)))
#' boundary_is_masked(b, c(10, 85))
#' @export
boundary_is_masked <- function(b, arc) {
  .check_boundary(b)
  arc <- as.numeric(arc)
  out <- rep(FALSE, length(arc))
  for (m in b$mask_arcs) {
    out <- out | (!is.na(arc) & arc >= m[1] & arc <= m[2])
  }
  if (length(b$mask_geoms)) {
    pts <- .line_interpolate(b$segs, arc)
    for (g in b$mask_geoms) {
      d <- .geom_distance(g, pts)
      out <- out | (!is.na(d) & d <= b$tol)
    }
  }
  out
}

#' Equally spaced positions along the boundary
#'
#' @inheritParams boundary_geometry
#' @param n Number of positions (default 25) when `spacing` is not given.
#' @param spacing Distance between consecutive positions.
#' @param start,end Arc-length range (default: the whole line).
#' @return Numeric vector of arc-length positions.
#' @examples
#' b <- boundary(rbind(c(0, 0), c(100, 0)), treated_side = "left")
#' boundary_grid(b, n = 5)
#' boundary_grid(b, spacing = 30)
#' @export
boundary_grid <- function(b, n = NULL, spacing = NULL, start = 0, end = NULL) {
  .check_boundary(b)
  # on a closed ring the end is the start again: leave it out
  ring <- isTRUE(b$closed) && is.null(end) && start == 0
  end <- end %||% b$length
  if (!is.null(spacing)) {
    out <- seq(start, end + 1e-9, by = spacing)
    return(if (ring) out[out < end - 1e-9] else out)
  }
  k <- n %||% 25
  if (ring) return(start + (end - start) * (seq_len(k) - 1) / k)
  seq(start, end, length.out = k)
}
