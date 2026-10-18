#' Distance-decay kernels
#'
#' Kernel functions evaluated at `u = d / h`, with `d` a network distance.
#'
#' @param kernel `"triangular"`, `"epanechnikov"`, `"uniform"`, `"gaussian"`,
#'   their three-letter abbreviations (`"tri"`, `"epa"`, `"uni"`, `"gau"`), or
#'   a function, which is returned unchanged.
#' @return A vectorised function of `u`.
#' @examples
#' k <- get_kernel("epanechnikov")
#' k(c(0, 0.5, 1.5))
#' @export
get_kernel <- function(kernel) {
  if (is.function(kernel)) {
    return(kernel)
  }
  name <- .kernel_name(kernel)
  switch(name,
    triangular = function(u) {
      u <- abs(as.numeric(u))
      ifelse(!is.na(u) & u <= 1, 1 - u, 0)
    },
    epanechnikov = function(u) {
      u <- abs(as.numeric(u))
      ifelse(!is.na(u) & u <= 1, 0.75 * (1 - u^2), 0)
    },
    uniform = function(u) {
      u <- abs(as.numeric(u))
      ifelse(!is.na(u) & u <= 1, 0.5, 0)
    },
    gaussian = function(u) {
      u <- as.numeric(u)
      ifelse(is.finite(u), exp(-0.5 * u^2) / sqrt(2 * pi), 0)
    }
  )
}

.KERNELS <- c(triangular = "triangular", tri = "triangular", epanechnikov = "epanechnikov",
              epa = "epanechnikov", uniform = "uniform", uni = "uniform", gaussian = "gaussian",
              gau = "gaussian")

.kernel_name <- function(kernel) {
  if (!is.character(kernel) || length(kernel) != 1 || !tolower(kernel) %in% names(.KERNELS)) {
    stop(sprintf("unknown kernel; choose from %s", paste(sort(names(.KERNELS)), collapse = ", ")),
         call. = FALSE)
  }
  unname(.KERNELS[tolower(kernel)])
}
