# Compose karst's arms of the missegregation-rate round (ac36fd1) from
# BENCH_RESULTS/tables/beta_picks.tsv (scripts/beta_picks.R), and write the
# task table of the picked candidates' forecasts.
#
#   Rscript scripts/compose_beta.R TASK_DIR
#
# For each input, fits/<source>/<dataset>/<design>/:
#   karst_beta_auto   relative link to the picked candidate,
#                      ../../../../candidates/<...>/<label>, or to the sibling
#                      karst where the pick is the reused 5e-5 fit; with no
#                      pick (every candidate failed) to the 5e-5 candidate,
#                      whose failed fit it then reports
#   karst_beta_fixed  missegregation simulations only: the 5e-5 candidate
#                      (their karst fits take the true rate)
# Existing links of these names are replaced; anything else is refused. An
# input with a missing candidate (no run.Rds: rerun its task) stops the
# script, so that no arm is composed from part of the grid.
#
# TASK_DIR/beta_forecast_mem<tier>.tsv: the ABM and Salehi inputs whose pick
# is not the reused 5e-5 fit (that one has its forecasts), as tasks of method
# karst_beta_auto with pm = the picked rate (beta_multiple), at the
# candidate's memory tier: run_task.R finds the fit done and runs the
# forecasts and the ABM scores in the picked candidate's folder. Salehi
# trajectories without a daughter passage have no forecast and no task.
# TASK_DIR/beta_forecast_tiers.tsv lists them.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
task_dir <- args[[1]]
picks <- utils::read.delim(file.path(paths$results, "tables", "beta_picks.tsv"), stringsAsFactors = FALSE)
tasks <- utils::read.delim(file.path(task_dir, "beta_tasks.tsv"), stringsAsFactors = FALSE,
                           colClasses = c(minobs = "character"))
if (any(picks$candidates_missing > 0)) {
  print(utils::head(picks[picks$candidates_missing > 0, c("source", "dataset", "design", "candidates_missing")], 30))
  stop(sum(picks$candidates_missing > 0), " inputs have a missing candidate; rerun their tasks first.")
}

# (Sys.readlink() gives "" for a file that is not a link and NA, on some
# systems, for a path that does not exist)
is_link <- function(path) {
  l <- Sys.readlink(path)
  !is.na(l) && nzchar(l)
}
link_to <- function(base, name, target) {
  link <- file.path(base, name)
  if (is_link(link)) file.remove(link) else if (file.exists(link)) stop("Not a link, refusing to replace: ", link)
  stopifnot(file.symlink(target, link))
}
target_of <- function(p, label) {
  if (label == "x1" && p$reused_x1) return("karst")
  file.path("..", "..", "..", "..", "candidates", p$source, p$dataset, p$design, label)
}
made <- c(auto = 0L, fixed = 0L)
for (i in seq_len(nrow(picks))) {
  p <- picks[i, ]
  base <- file.path(paths$results, "fits", p$source, p$dataset, p$design)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  target <- target_of(p, if (is.na(p$pick_label)) "x1" else p$pick_label)
  stopifnot(dir.exists(file.path(base, target)))
  link_to(base, "karst_beta_auto", target)
  made["auto"] <- made["auto"] + 1L
  if (p$source == "misseg") {
    link_to(base, "karst_beta_fixed", target_of(p, "x1"))
    made["fixed"] <- made["fixed"] + 1L
  }
}
cat("links: karst_beta_auto", made[["auto"]], " karst_beta_fixed", made[["fixed"]], "\n")

## Forecasts of the picked candidates ---------------------------------------------------------
has_forecast <- function(p) {
  # Salehi: a trajectory with a daughter passage; the reference fit has a
  # forecast.Rds (ok or failed) then
  if (p$source != "salehi") return(TRUE)
  length(Sys.glob(file.path(paths$results, "fits", "salehi", p$dataset, p$design, "karst", "minobs_*",
                            "forecast.Rds"))) > 0
}
need <- picks[picks$source %in% c("abm", "salehi") & !is.na(picks$pick_label) &
                !(picks$pick_label == "x1" & picks$reused_x1) & picks$pick_outputs_ok > 0, ]
need <- need[vapply(seq_len(nrow(need)), function(i) has_forecast(need[i, ]), logical(1)), ]
info <- unique(tasks[c("source", "dataset", "design", "minobs", "n0", "nb", "mem")])
fc <- merge(need[c("source", "dataset", "design", "pick_multiple")], info, by = c("source", "dataset", "design"))
fc <- data.frame(source = fc$source, dataset = fc$dataset, design = fc$design, method = "karst_beta_auto",
                 minobs = fc$minobs, n0 = fc$n0, nb = fc$nb, pm = karst_beta_default * fc$pick_multiple,
                 beta_multiple = fc$pick_multiple, mem = fc$mem, stringsAsFactors = FALSE)
fc <- fc[order(fc$source, fc$dataset, fc$design), ]
old <- list.files(task_dir, pattern = "^beta_forecast_.*\\.tsv$", full.names = TRUE)
if (length(old)) invisible(file.remove(old))
listing <- list()
for (m in unique(fc$mem)) {
  part <- fc[fc$mem == m, ]
  part <- cbind(task = seq_len(nrow(part)), part)
  file <- file.path(task_dir, paste0("beta_forecast_mem", m, ".tsv"))
  utils::write.table(part, file, sep = "\t", quote = FALSE, row.names = FALSE)
  listing[[length(listing) + 1L]] <- data.frame(file = file, tasks = nrow(part), mem = m, time = "06:00:00")
}
if (length(listing)) {
  listing <- do.call(rbind, listing)
  utils::write.table(listing, file.path(task_dir, "beta_forecast_tiers.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE)
  print(transform(listing, file = basename(file)), row.names = FALSE)
}
cat("forecast tasks:", nrow(fc), "\n")
print(table(fc$source, fc$design))
