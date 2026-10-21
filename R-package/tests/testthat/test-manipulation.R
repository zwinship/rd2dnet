# Reference values from the Python package (same clusters and the same portable draws).

test_that("manipulation test matches the Python package", {
  cc <- network_coordinates(mock_city()$network("length"), mock_city()$boundary, mock_city()$obs, n_eval = 9)
  r <- manipulation_test(cc, q = c(50, 200), reps = 2000)
  expect_equal(r$path$statistic, c(1.1666666666666667, 1.1666666666666667), tolerance = 1e-12)
  expect_equal(r$path$pvalue, c(0.2375, 0.689))
  expect_equal(r$path$n_treated, c(23, 100))
  expect_equal(r$q, 50L)
})

test_that("the portable draws match the Python package", {
  expect_equal(rd2dnet:::.lcg_labels(0, 2, 10),
               rbind(c(1, 1, 0, 1, 0, 1, 1, 0, 0, 0), c(1, 0, 0, 1, 1, 0, 0, 1, 1, 1)), ignore_attr = TRUE)
})

test_that("the 1-D k-means partitions are optimal", {
  set.seed(1)
  for (trial in 1:50) {
    pos <- sort(sample(0:11, sample(3:8, 1), replace = TRUE))
    vals <- unique(pos)
    parts <- rd2dnet:::.kmeans_1d_partitions(pos, length(vals))
    sse <- function(firsts) {
      ends <- c(firsts[-1], length(pos) + 1)
      sum(mapply(function(a, e) sum((pos[a:(e - 1)] - mean(pos[a:(e - 1)]))^2), firsts, ends))
    }
    cands <- which(diff(pos) != 0) + 1
    for (k in seq_along(parts)) {
      best <- if (k == 1) sse(1) else min(apply(utils::combn(cands, k - 1, simplify = TRUE), 2,
                                                function(cmb) sse(c(1, cmb))))
      if (k > 1 && length(cands) == k - 1) best <- sse(c(1, cands))
      expect_equal(sse(parts[[k]]), best, tolerance = 1e-9)
    }
  }
})

test_that("density uses the straight-line distance by default", {
  d <- density(mock_coords())
  expect_equal(d$normalize, "euclidean")
  expect_equal(d$diff_bc, -1.16812813223094e-05, tolerance = 1e-8)
  expect_equal(d$se, 0.00010308279614672041, tolerance = 1e-8)
})
