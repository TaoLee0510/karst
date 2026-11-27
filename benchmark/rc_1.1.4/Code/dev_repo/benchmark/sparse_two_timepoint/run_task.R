#!/usr/bin/env Rscript
# One fit of the development benchmark (one line of a make_tasks.R table), as
# run.R --stage=fit fits it.
#   Rscript benchmark/sparse_two_timepoint/run_task.R TASKS.tsv ID
# Run from the repository root, with karst on the library path.
script_dir <- file.path("benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R", "generate.R", "methods.R", "score.R")) source(file.path(script_dir, part))
args <- commandArgs(TRUE)
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
config <- stt_config()
path <- file.path(config$output_dir, "fits", sprintf("%s__%s.rds", task$dataset_id, task$method))
if (file.exists(path)) quit(save = "no")
cat("karst", as.character(utils::packageVersion("karst")), find.package("karst"), "\n")
dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
nodes <- readRDS(file.path(config$output_dir, "nodes", paste0(task$dataset_id, ".rds")))
work_dir <- file.path(tempdir(), paste0(task$dataset_id, "__", task$method))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(dataset$seed)
result <- stt_run_method(config$methods[[task$method]], dataset, nodes, work_dir)
result$dataset_id <- task$dataset_id
result$method <- task$method
temporary <- paste0(path, ".tmp", Sys.getpid())
saveRDS(result, temporary)
file.rename(temporary, path)
cat(task$dataset_id, task$method, result$status, round(result$seconds), "s\n")
