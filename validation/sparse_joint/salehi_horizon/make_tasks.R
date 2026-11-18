#!/usr/bin/env Rscript
# Task table of forecast_score.R: the Salehi fits of a round of the paper
# benchmark (benchmark/paper) at minobs 5, 10 and 20. Design "ends" (the
# first and last passage of every trajectory): KARST at its default
# (initial_state = "balance") and alfak, from the round's results
# (BENCH_RESULTS; alfak's fits are linked in from the alfak2 era), and KARST
# with initial_state = "free" (method karst_free, fitted by the same round
# into a folder of its own, BENCH_RESULTS_FREE, default BENCH_RESULTS
# followed by "_free"). Design "all" (every passage of the trajectories of
# three or more passages): KARST anchored (the default), pairwise and
# state_space, alfak, and KARST "free". A task without a fit or without
# descendants records that and stops.
#   Rscript make_tasks.R OUT.tsv [ends | all] [METHOD ...]   (from benchmark/paper)
# METHOD keeps those methods' tasks only (e.g. KARST's, when alfak's rows are
# carried over from an earlier run with reuse_rows.R).
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
design <- if (length(args) > 1L) args[[2]] else "ends"
stopifnot(design %in% c("ends", "all"))
ids <- sort(salehi_trajectories(paths)$traj_id)
round <- file.path(paths$results, "fits", "salehi")
free <- file.path(Sys.getenv("BENCH_RESULTS_FREE", paste0(paths$results, "_free")), "fits", "salehi")
roots <- c(karst = round, alfak = round, karst_free = free)
if (design == "all") {
  roots <- c(karst = round, karst_pairwise = round, karst_state_space = round, alfak = round,
             karst_free = free)
  ids <- ids[as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", ids)) >= 3L]
}
if (length(args) > 2L) {
  keep <- args[-(1:2)]
  stopifnot(all(keep %in% names(roots)))
  roots <- roots[keep]
}
tasks <- expand.grid(id = ids, method = names(roots), minobs = c(5, 10, 20), stringsAsFactors = FALSE)
tasks$design <- design
tasks$root <- roots[tasks$method]
tasks$pm <- bench_settings$salehi$pm
tasks$nb <- bench_settings$salehi$nb
utils::write.table(tasks, args[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
cat("tasks", nrow(tasks), "\n")
print(table(tasks$method))
