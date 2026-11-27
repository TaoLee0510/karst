# The paper's Fig. 6 (08_misseg_rate_analysis.R) on the Salehi fits of one
# arm (R/paper_fig456.R): for each sample (lineage start), the landscape with
# the highest xv (xv > 0, a test passage), and its steady-state karyotype
# distribution under missegregation rates p = 2e-5 * 2^(1:8), from
# alfakR::find_steady_state(lscape, p, Nmax = 1).
#
#   Rscript scripts/fig6_steady_states.R --index N     (arm N of fig456_arms)
#   Rscript scripts/fig6_steady_states.R ARM_KEY
#
# FIG456_EXCLUDE=0 keeps the trajectories that are not analysable (the paper
# used all). Writes BENCH_RESULTS/analysis/fig456/fig6{,_all_inputs}/<arm>.Rds:
# a list named by trajectory of karyotype x rate matrices (rows ordered by
# their sum, attributes fi and min_obs), as the paper's
# steadyStatePredictions.Rds.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
arm_i <- if (identical(args[1], "--index")) as.integer(args[2]) else match(args[1], fig456_arms$key)
stopifnot(!is.na(arm_i), arm_i <= nrow(fig456_arms))
arm <- fig456_arms[arm_i, ]
exclude <- Sys.getenv("FIG456_EXCLUDE", "1") != "0"
out <- file.path(analysis_root(paths), "fig456", if (exclude) "fig6" else "fig6_all_inputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cat("Fig. 6, arm", arm$key, " exclusion", exclude, "\n")

rows <- fig456_fit_rows(arm, paths, exclude)
x <- do.call(rbind, lapply(rows, function(r) r[c("fi", "min_obs", "xv", "test_treat")]))
x <- x[!is.na(x$xv) & x$xv > 0 & x$test_treat != "NaN", ]
# the best minobs of each trajectory, then the best trajectory of each start,
# ties resolved (R/select.R)
x <- select_by_cv(x, "fi", "min_obs", "xv")
x <- x[!duplicated(x$fi), ]
rownames(x) <- NULL
x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- select_by_cv(x, "clean_id", "fi", "xv")
x <- x[!duplicated(x$fi), ]
cat("landscapes:", nrow(x), "\n")
p <- 0.00002 * 2^(1:8)
# version B (ALIGN=karst): landscapes on E, too long for alfakR's matrix
# builder; the same steady state by align_find_steady_state() (R/align.R),
# and only the karyotypes with a share above zero at some rate are kept
align <- align_mode() == "karst"
steady_state <- if (align) align_find_steady_state else
  function(lscape, p) alfakR::find_steady_state(lscape = lscape, p = p, Nmax = 1)
started <- proc.time()[["elapsed"]]
df <- lapply(seq_len(nrow(x)), function(i) {
  lscape <- fig456_landscape(fig456_fit_dir(arm, x$fi[i], x$min_obs[i], paths))
  if (is.null(lscape)) return(NULL)
  res <- do.call(cbind, lapply(p, function(p0) steady_state(lscape, p0)))
  colnames(res) <- p
  res <- res[order(rowSums(res), decreasing = TRUE), ]
  if (align) res <- res[rowSums(res) > 0, , drop = FALSE]
  attr(res, "min_obs") <- x$min_obs[i]
  attr(res, "fi") <- x$fi[i]
  cat(sprintf("%d %s: %d karyotypes, %d dominant, %.0f s\n", i, x$fi[i], nrow(res),
              length(unique(apply(res, 2, which.max))), proc.time()[["elapsed"]] - started))
  res
})
names(df) <- x$fi
df <- Filter(Negate(is.null), df)
saveRDS(df, file.path(out, paste0(arm$key, ".Rds")))
cat("Done.\n")
