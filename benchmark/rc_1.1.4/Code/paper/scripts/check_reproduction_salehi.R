# Check that the Salehi scoring code reproduces the paper's scores.
#
#   Rscript scripts/check_reproduction_salehi.R [--out DIR] [--index N | TRAJ_ID ...]
#
# For each of the paper's Salehi fits (PAPER_DATA/processed/salehi/
# alfak_outputs/<traj>/minobs_<m>) the rows of 03_process_alfak_fits.R are
# computed again and compared with the paper's alfak_outputs_proc/
# <traj>_minobs_<m>.Rds. `--index N` takes the N-th trajectory in sorted
# order. With `--out`, a summary per fit is saved to DIR/<traj>.Rds. Stops if
# a number differs by more than 1e-8 or the rows do not line up.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
out <- NULL
if (identical(args[1], "--out")) {
  out <- args[2]
  args <- args[-(1:2)]
}
all_traj <- sort(salehi_trajectories(paths)$traj_id)
traj <- if (identical(args[1], "--index")) all_traj[as.integer(args[2])] else if (length(args)) args else all_traj
lineages <- readRDS(paths$salehi_lineages)
meta <- utils::read.csv(paths$salehi_meta)
proc <- file.path(paths$ref, "processed", "salehi", "alfak_outputs_proc")
numbers <- c("xv", "pred0", "pred", "base", "angle")
labels <- c("fi", "min_obs", "ntrain", "t_train", "train_treat", "metric", "tt", "test_treat")

rows <- list()
for (id in traj) {
  for (d in Sys.glob(file.path(paths$paper_salehi, id, "minobs_*"))) {
    saved_path <- file.path(proc, paste0(id, "_", basename(d), ".Rds"))
    if (!file.exists(saved_path)) next
    saved <- readRDS(saved_path)
    mine <- salehi_fit_scores(d, id, lineages, meta, paths$salehi_inputs)
    aligned <- nrow(mine) == nrow(saved) &&
      all(vapply(labels, function(n) identical(as.character(mine[[n]]), as.character(saved[[n]])), TRUE))
    diff <- if (aligned) max(vapply(numbers, function(n) {
      a <- as.numeric(mine[[n]])
      b <- as.numeric(saved[[n]])
      if (!identical(is.na(a), is.na(b))) return(Inf)
      if (all(is.na(a))) 0 else max(abs(a - b), na.rm = TRUE)
    }, 0)) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(traj = id, fit = basename(d), rows = nrow(saved), aligned = aligned,
                                            max_diff = diff)
  }
}
rows <- do.call(rbind, rows)
if (is.null(rows)) {
  cat("No paper fits for", paste(traj, collapse = ", "), "\n")
  quit(save = "no", status = 0)
}
if (!is.null(out)) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(rows, file.path(out, paste0(if (length(traj) == 1L) traj else "all", ".Rds")))
}
print(rows)
bad <- rows[!rows$aligned | rows$max_diff > 1e-8, ]
if (nrow(bad)) stop("Scores differ from the paper's.")
cat("All Salehi scores match the paper's.\n")
