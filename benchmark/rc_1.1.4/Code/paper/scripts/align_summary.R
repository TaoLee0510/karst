# Sizes, times and checks of version B's tasks (scripts/align_task.R), from
# their task.Rds:
#
#   ALIGN=karst Rscript scripts/align_summary.R     (ALIGN_RESULTS: another folder)
#
# Prints, by source and design: inputs done, the extended list's size, the
# time and peak memory of a task, and per fit kind the time of its landscape,
# cross-validation and forecasts; the share of alfak's landscape that Kriging
# fills and the agreement of the Kriging with alfak's mean on its own list;
# the largest difference of KARST's landscape on its own list from the round's
# (check_original). Writes align_tasks.Rds into the folder.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
if (nzchar(Sys.getenv("ALIGN_RESULTS"))) paths$results <- Sys.getenv("ALIGN_RESULTS")
# one record per input (the pilot's first layout) or per landscape
files <- c(Sys.glob(file.path(paths$results, "fits", "*", "*", "*", "task.Rds")),
           Sys.glob(file.path(paths$results, "fits", "*", "*", "*", "*", "task_*.Rds")))
recs <- lapply(files, readRDS)
inputs <- do.call(rbind, lapply(recs, function(r) data.frame(source = r$source, dataset = r$dataset, design = r$design,
  E = r$E, observed = r$observed, seconds = r$seconds, peak_mb = r$peak_mb, stringsAsFactors = FALSE)))
fits <- do.call(rbind, lapply(recs, function(r) if (is.null(r$fits)) NULL else
  cbind(source = r$source, dataset = r$dataset, design = r$design, E = r$E, r$fits, stringsAsFactors = FALSE)))
options(width = 200)
cat("tasks done:", length(files), "\n")
q <- function(x) paste(round(stats::quantile(x, c(0, 0.5, 1), na.rm = TRUE)), collapse = " / ")
print(do.call(rbind, lapply(split(inputs, inputs[c("source", "design")], drop = TRUE), function(x) data.frame(
  source = x$source[1], design = x$design[1], inputs = nrow(x), E_min_med_max = q(x$E),
  seconds_min_med_max = q(x$seconds), peak_mb_max = round(max(x$peak_mb, na.rm = TRUE))))), row.names = FALSE)
if (!is.null(fits)) {
  fits$kind <- ifelse(fits$method == "alfak", "alfak (per minobs)", fits$method)
  print(do.call(rbind, lapply(split(fits, fits[c("source", "kind")], drop = TRUE), function(x) data.frame(
    source = x$source[1], fit = x$kind[1], fits = nrow(x), landscape_s = q(x$landscape_seconds),
    xval_s = q(x$xval_seconds), forecast_s = q(x$forecast_seconds),
    kriged_share = round(stats::median(x$kriged / x$E, na.rm = TRUE), 3),
    kriging_agreement_min = round(suppressWarnings(min(x$kriging_agreement, na.rm = TRUE)), 3),
    check_original_max = signif(suppressWarnings(max(x$check_original, na.rm = TRUE)), 3)))), row.names = FALSE)
}
saveRDS(list(inputs = inputs, fits = fits), file.path(paths$results, "align_tasks.Rds"))
