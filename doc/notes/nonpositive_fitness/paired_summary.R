
# Run from the repository root:  Rscript doc/notes/nonpositive_fitness/paired_summary.R
# Notes: doc/notes/nonpositive_fitness.md. KARST_OWN_ROUND names the own benchmark's
# round (default 6ae2a24); results are read from and written to round_<round>/ of this
# folder. OWN_RESULTS is that round's held-out results (datasets/, nodes/, fits/,
# dataset_index_heldout.csv); the default is the copy made by
# doc/report/build/fetch_own.sh. NF_HERE is this folder (default
# doc/notes/nonpositive_fitness), NF_CORES the workers (default SLURM_CPUS_PER_TASK, else 3).
here <- Sys.getenv("NF_HERE", "doc/notes/nonpositive_fitness")
round <- Sys.getenv("KARST_OWN_ROUND", "6ae2a24")
own <- Sys.getenv("OWN_RESULTS", file.path("doc/report/data", paste0("own_", round)))
out_dir <- file.path(here, paste0("round_", round))
dir.create(out_dir, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("NF_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "3")))
out <- readRDS(file.path(out_dir, "paired_experiment.rds"))
id <- do.call(rbind, lapply(out, `[[`, "id"))
pop <- do.call(rbind, lapply(out, function(o) cbind(o$id, o$pop)))
ds <- do.call(rbind, lapply(out, function(o) cbind(o$id, o$ds)))
cat("Deaths under the paper rule over 2520 days, by wavelength (median, max over the runs):\n")
print(aggregate(died_total ~ wavelength, id, function(x) c(median = median(x), max = max(x), zero = sum(x == 0))))
cat("\nPopulations, persist vs die (same seeds): TV distance, by wavelength and day\n")
print(aggregate(tv_rules ~ wavelength + day, pop, function(x) signif(c(median = median(x), max = max(x), n_gt_0.01 = sum(x > 0.01)), 3)))
cat("\nShare of f<=0 cells in the persist population\n")
print(aggregate(nonpos_persist ~ wavelength, pop, function(x) signif(c(median = median(x), max = max(x)), 3)))
ds$core <- !grepl("^turned", ds$scenario)
cat("\nDatasets with identical counts under both rules\n")
print(aggregate(counts_identical ~ wavelength + core, ds, function(x) c(identical = sum(x), of = length(x))))
cat("\nDatasets whose counts differ: cells that differ (of 100-6000)\n")
x <- ds[!ds$counts_identical, ]; print(x[order(-x$cells_differing), c("landscape_id","patient","scenario","cells_differing","tv_truth_t1")][1:min(25, nrow(x)), ], row.names = FALSE)
cat("\nTV of the t1 population between rules, core scenarios\n"); print(aggregate(tv_truth_t1 ~ wavelength, ds[ds$core, ], function(x) signif(c(median = median(x), max = max(x)), 3)))
cat("\nturned scenarios\n"); print(aggregate(cbind(tv_truth_t1, cells_differing) ~ wavelength, ds[!ds$core, ], function(x) signif(c(median = median(x), max = max(x)), 3)))
