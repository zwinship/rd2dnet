# Runs the R package on every parity case written by make_cases.py.
#   Rscript parity/run_r.R [parity/cases]
suppressPackageStartupMessages(library(rd2dnet))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else file.path("parity", "cases")

status_code <- function(s) {
  codes <- c(ok = 0, masked = 1, insufficient = 2, singular = 3, "no bandwidth" = 4,
             "weak first stage" = 5, "few clusters" = 6)
  unname(codes[as.character(s)])
}

write_res <- function(t, file) {
  basis <- if ("basis" %in% names(t)) c("s,o" = 1, s = 2)[as.character(t$basis)] else NA
  out <- data.frame(h0 = t$h0, h1 = t$h1, n0 = t$n0, n1 = t$n1, tau = t$tau, se = t$se,
                    tau_bc = t$tau_bc, se_rb = t$se_rb, status = status_code(t$status),
                    basis = unname(basis))
  write.csv(out, file, row.names = FALSE)
}

cases <- readLines(file.path(root, "cases.txt"))
failed <- character()
for (cs in cases) {
  d <- file.path(root, cs)
  ok <- tryCatch({
    eval(parse(file.path(d, "run.R")), envir = environment())
    TRUE
  }, error = function(e) {
    message("case ", cs, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) failed <- c(failed, cs)
  cat("case", cs, if (ok) "done" else "ERROR", "\n")
}
if (length(failed)) {
  cat("R errors in cases:", failed, "\n")
  quit(status = 1)
}
