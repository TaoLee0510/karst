# Scores of the Salehi fits of one trajectory (03_process_alfak_fits.R).
#
#   Rscript scripts/score_salehi.R --index N | TRAJ_ID ...
#
# `--index N` takes the N-th trajectory in sorted order (one array task per
# trajectory). For every fit directory of the trajectory (each method,
# design and minobs) the rows of 03 are computed and saved together, with
# the method and design added, to BENCH_RESULTS/scores/salehi/<traj>.Rds.
# SCORE_METHODS (comma-separated, e.g. karst_auto) scores the folders of those
# methods only and merges them into the trajectory's existing file, replacing
# that method's rows (e.g. karst's "auto" arm, linked in after the rest was
# scored).

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
all_traj <- sort(salehi_trajectories(paths)$traj_id)
traj <- if (identical(args[1], "--index")) all_traj[as.integer(args[2])] else args
lineages <- readRDS(paths$salehi_lineages)
meta <- utils::read.csv(paths$salehi_meta)
out <- file.path(paths$results, "scores", "salehi")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

only <- strsplit(Sys.getenv("SCORE_METHODS", ""), ",", fixed = TRUE)[[1]]
for (id in traj) {
  base <- file.path(paths$results, "fits", "salehi", id)
  dirs <- Sys.glob(file.path(base, "*", "*", "minobs_*"))
  if (length(only)) dirs <- dirs[basename(dirname(dirs)) %in% only]
  rows <- lapply(dirs, function(d) {
    if (!file.exists(file.path(d, "run.Rds"))) return(NULL)
    parts <- strsplit(sub(paste0(base, "/"), "", d, fixed = TRUE), "/", fixed = TRUE)[[1]]
    s <- salehi_fit_scores(d, id, lineages, meta, paths$salehi_inputs)
    cbind(design = parts[1], method = parts[2], s, stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  path <- file.path(out, paste0(id, ".Rds"))
  if (length(only) && file.exists(path)) {
    old <- readRDS(path)
    old <- old[!old$method %in% only, , drop = FALSE]
    rows <- if (is.null(rows)) old else rbind(old, rows[names(old)])
  }
  saveRDS(rows, path)
  cat(id, ":", length(dirs), "fit directories,", if (is.null(rows)) 0 else nrow(rows), "rows\n")
}
