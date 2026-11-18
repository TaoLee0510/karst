#!/usr/bin/env Rscript
# Task table of run_task.R for a rerun of the development benchmark with
# karst: the fits of an earlier run (FITS_FROM, a fits/ folder of files
# <dataset>__<method>.rds) whose method karst still has (config.R) are refitted;
# alfak's fits (alfak_minobs5, alfak_minobs20) and the null predictor are
# copied, not refitted. The earlier run is the alfak2 era's, whose methods were
# named alfak2_* and alfakR_*; they map to karst_* and alfak_*.
#
#   Rscript benchmark/sparse_two_timepoint/make_tasks.R FITS_FROM TASKS.tsv
# Run from the repository root; copies into benchmark/results/sparse_two_timepoint/fits.
for (part in c("landscape.R", "config.R")) source(file.path("benchmark", "sparse_two_timepoint", part))
args <- commandArgs(TRUE)
from <- args[[1]]
config <- stt_config()
files <- list.files(from, pattern = "__.*[.]rds$")
dataset <- sub("__.*$", "", files)
old <- sub("[.]rds$", "", sub("^.*__", "", files))
new <- sub("^alfak2_", "karst_", old)
new <- sub("^alfakR_", "alfak_", new)
split <- ifelse(startsWith(dataset, "heldout_"), "heldout", "dev")
known <- new %in% names(config$methods)
cat("Earlier fits:", length(files), " of methods karst no longer has (left out):", sum(!known), "\n")
print(table(old[!known]))
copy <- known & (startsWith(new, "alfak_") | new == "zero")
refit <- known & !copy
out <- file.path(config$output_dir, "fits")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
copied <- 0L
for (i in which(copy)) {
  target <- file.path(out, sprintf("%s__%s.rds", dataset[i], new[i]))
  if (file.exists(target)) next
  fit <- readRDS(file.path(from, files[i]))
  fit$method <- new[i]
  saveRDS(fit, target)
  copied <- copied + 1L
}
tasks <- data.frame(split = split[refit], dataset_id = dataset[refit], method = new[refit], stringsAsFactors = FALSE)
tasks <- tasks[order(tasks$method, tasks$dataset_id), ]
utils::write.table(tasks, args[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
cat("Copied (alfak and zero):", copied, " tasks (karst):", nrow(tasks), "\n")
print(table(tasks$method, tasks$split))
