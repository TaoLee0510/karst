#!/usr/bin/env Rscript
# Add landscapes to a split whose datasets already exist, without touching the
# others: the truth population, the datasets of the scoreboard patient (every
# scenario) and their scored nodes, as run.R --stage=truth and --stage=datasets
# make them, for the named landscapes only. run.R cannot do this when the
# truth populations of the earlier landscapes are not at hand, because its
# dataset stage reads them for every landscape of the split.
#
#   Rscript benchmark/sparse_two_timepoint/add_landscapes.R --split=heldout --unseal-heldout \
#     --landscapes=heldout_L13,heldout_L14,heldout_L15 [--task=N | --index-only] [--workers=1]
#
# --task=N (or SLURM_ARRAY_TASK_ID) makes only the N-th listed landscape, so that
# one HPC task makes one landscape. --index-only adds the listed landscapes' datasets
# to dataset_index_<split>.csv (replacing earlier rows of theirs), once all are made.
# Run from the repository root.
script_dir <- file.path("benchmark", "sparse_two_timepoint")
# run.R's functions (it runs its own main() only when called as a script)
source(file.path(script_dir, "run.R"))
args <- commandArgs(TRUE)
option <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  if (!grepl("=", hit[[1L]], fixed = TRUE)) return(TRUE)
  sub("^[^=]*=", "", hit[[1L]])
}
config <- stt_config()
split <- option("split")
if (is.null(split)) stop("--split is required", call. = FALSE)
if (identical(split, "heldout") && !isTRUE(option("unseal-heldout"))) {
  stop("Held-out landscapes are reserved for the final confirmation; pass --unseal-heldout to use them.", call. = FALSE)
}
ids <- strsplit(option("landscapes", ""), ",", fixed = TRUE)[[1L]]
landscapes <- config$landscapes[config$landscapes$split == split & config$landscapes$landscape_id %in% ids, , drop = FALSE]
if (!length(ids) || nrow(landscapes) != length(ids)) stop("--landscapes must name landscapes of split ", split, call. = FALSE)
patient <- config$scoreboard_patients
if (length(patient) != 1L) stop("expected one scoreboard patient", call. = FALSE)

# The same seeds as run.R's dataset plan: scenario s of patient p on a landscape.
dataset_seed <- function(landscape, s) as.integer(landscape$landscape_seed %% 100000L * 1000L + patient * 20L + s)
dataset_rows <- function(landscape) {
  do.call(rbind, lapply(seq_len(nrow(config$scenarios)), function(s) {
    scenario <- config$scenarios[s, ]
    dataset <- readRDS(stt_dataset_path(config, sprintf("%s_P%02d_%s", landscape$landscape_id, patient, scenario$scenario)))
    data.frame(dataset_id = dataset$dataset_id, split = dataset$split, landscape_id = dataset$landscape_id,
               wavelength_class = dataset$wavelength_class, patient = dataset$patient,
               scenario = dataset$scenario, role = scenario$role, scoreboard = TRUE,
               t(stt_count_profile(dataset$counts)), stringsAsFactors = FALSE)
  }))
}

if (isTRUE(option("index-only"))) {
  path <- file.path(config$output_dir, sprintf("dataset_index_%s.csv", split))
  index <- utils::read.csv(path, stringsAsFactors = FALSE)
  added <- do.call(rbind, lapply(seq_len(nrow(landscapes)), function(i) dataset_rows(landscapes[i, ])))
  if (!setequal(names(index), names(added))) stop("columns of ", path, " differ from the new rows", call. = FALSE)
  index <- rbind(index[!index$landscape_id %in% ids, , drop = FALSE], added[names(index)])
  temporary <- paste0(path, ".tmp", Sys.getpid())
  utils::write.csv(index, temporary, row.names = FALSE)
  file.rename(temporary, path)
  cat("Index", path, ":", nrow(added), "rows of", paste(ids, collapse = ","), "added;", nrow(index), "in all\n")
  print(added[c("dataset_id", "role", setdiff(names(added), c("dataset_id", "split", "landscape_id", "wavelength_class",
                                                              "patient", "scenario", "role", "scoreboard")))],
        digits = 3, row.names = FALSE)
  quit(save = "no")
}

task <- option("task", Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
if (nzchar(task)) landscapes <- landscapes[as.integer(task), , drop = FALSE]
# Compile the simulator into a folder of this process: tasks may start together.
cache <- file.path(tempdir(), "cpp_cache")
dir.create(cache, showWarnings = FALSE)
Rcpp::sourceCpp(file.path(script_dir, "abm.cpp"), cacheDir = cache, env = globalenv())
for (i in seq_len(nrow(landscapes))) {
  landscape <- landscapes[i, ]
  started <- proc.time()[["elapsed"]]
  truth_path <- stt_generate_truth(config, landscape, patient)
  truth <- readRDS(truth_path)
  cat(landscape$landscape_id, "truth", round(proc.time()[["elapsed"]] - started), "s\n")
  for (s in seq_len(nrow(config$scenarios))) {
    path <- stt_make_dataset(config, truth, config$scenarios[s, ], config$t0, dataset_seed(landscape, s))
    stt_dataset_nodes(config, readRDS(path))
  }
  cat(landscape$landscape_id, "datasets and nodes", round(proc.time()[["elapsed"]] - started), "s\n")
  print(dataset_rows(landscape)[c("dataset_id", "role")], row.names = FALSE)
}
