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

write_agg <- function(res, file) {
  agg <- data.frame(wbate = NA_real_, lbate_max = NA_real_, lbate_min = NA_real_, lbate_absmax = NA_real_)
  if (any(res$table$status == "ok") && !is.null(res$crit)) {
    agg <- data.frame(wbate = unname(wbate(res)["tau_bc"]), lbate_max = unname(lbate(res)["tau_bc"]),
                      lbate_min = unname(lbate(res, "min")["tau_bc"]),
                      lbate_absmax = unname(lbate(res, "absmax")["tau_bc"]))
  }
  write.csv(agg, file, row.names = FALSE)
}

write_res <- function(t, file) {
  basis <- if ("basis" %in% names(t)) c("s,o" = 1, s = 2)[as.character(t$basis)] else NA
  out <- data.frame(h0 = t$h0, h1 = t$h1, n0 = t$n0, n1 = t$n1, tau = t$tau, se = t$se,
                    tau_bc = t$tau_bc, se_rb = t$se_rb, status = status_code(t$status),
                    basis = unname(basis))
  write.csv(out, file, row.names = FALSE)
}

# Nearest-crossing bins (mirrors run_bins() in make_cases.py).
run_bins <- function(obs, cc, d, h) {
  b <- bin_coordinates(cc, features = "lanes")
  br <- b$branch
  canon <- ifelse(br > 0, match(br, unique(br[br > 0])), 0L)
  write.csv(data.frame(bin = b$bin, branch = canon, S = b$similarity), file.path(d, "r_bins_obs.csv"),
            row.names = FALSE)
  tab <- bins_table(b)
  write.csv(data.frame(arc = tab$arc, n0 = tab$n0, n1 = tab$n1, h_S = b$h_S), file.path(d, "r_bins_tab.csv"),
            row.names = FALSE)
  agg <- list()
  for (k in 1:2) {
    res <- if (k == 1) rdnet(obs$y_out, b, h = h, min_obs = 10, seed = 0) else rdnet(obs$y_out, b, min_obs = 10, seed = 0)
    write_res(res$table, file.path(d, sprintf("r_bins_res_%d.csv", k)))
    if (k == 1) {
      for (rule in c("density", "gls")) {
        a <- tryCatch(suppressWarnings(wbate(res, rule)), error = function(e) c(tau_bc = NA_real_, se_rb = NA_real_))
        agg[[paste0(rule, "_tau_bc")]] <- unname(a["tau_bc"])
        agg[[paste0(rule, "_se")]] <- unname(a["se_rb"])
      }
    }
  }
  fe <- tryCatch(rdnet_fe(obs$y_out, b, h = h), error = function(e) NULL)
  if (is.null(fe)) {
    agg[c("fe_tau", "fe_se", "fe_tau_bc", "fe_se_rb")] <- NA_real_
  } else {
    agg[c("fe_tau", "fe_se", "fe_tau_bc", "fe_se_rb")] <- list(fe$tau, fe$se, fe$tau_bc, fe$se_rb)
  }
  write.csv(as.data.frame(agg), file.path(d, "r_bins_agg.csv"), row.names = FALSE)
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
