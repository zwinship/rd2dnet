# Build data/*.rda from the CSVs written by data-raw/mock_midtown.py
# (run from R-package/: Rscript data-raw/mock_midtown.R).
mock_midtown_nodes <- read.csv("data-raw/mock_midtown_nodes.csv")
mock_midtown_edges <- read.csv("data-raw/mock_midtown_edges.csv", stringsAsFactors = FALSE)
mock_midtown_obs <- read.csv("data-raw/mock_midtown_obs.csv")
mock_midtown_obs$treated <- as.logical(mock_midtown_obs$treated)
mock_midtown_nodes$name <- as.integer(mock_midtown_nodes$name)
mock_midtown_edges$from <- as.integer(mock_midtown_edges$from)
mock_midtown_edges$to <- as.integer(mock_midtown_edges$to)
mock_midtown_obs$block <- as.integer(mock_midtown_obs$block)
save(mock_midtown_nodes, file = "data/mock_midtown_nodes.rda", compress = "xz")
save(mock_midtown_edges, file = "data/mock_midtown_edges.rda", compress = "xz")
save(mock_midtown_obs, file = "data/mock_midtown_obs.rda", compress = "xz")
