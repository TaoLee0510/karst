# Part 2 of PLAN.md: one fit of the own benchmark's held-out main table at a
# multiple of the missegregation rate. The round's own code and data
# (benchmark/sparse_two_timepoint of the HPC clone, round folder
# sparse_two_timepoint_6ae2a24, read only); only the rate given to the method
# changes (dataset$beta x multiple, used for KARST's beta and ALFA-K's pm), with
# the round's seed, so x1 reproduces the round. The benchmark code itself is not
# changed (a wrapper instead of the option the plan named).
#   Rscript own_task.R TASK_FILE ID OUT_DIR
# TASK_FILE columns: dataset_id, method, multiple. Writes OUT_DIR/x<mm>/fits/<dataset>__<method>.rds.
args <- commandArgs(TRUE)
# rc freeze: the own benchmark's code is Code/own, its inputs STT_DATA_DIR (data/own_benchmark)
script_dir <- file.path(Sys.getenv("RC_HOME"), "Code", "own", "benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R", "generate.R", "methods.R", "score.R")) source(file.path(script_dir, part))
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
config <- stt_config(output_dir = args[[3]], data_dir = Sys.getenv("STT_DATA_DIR"))
out <- file.path(args[[3]], sprintf("x%02d", task$multiple), "fits")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
path <- file.path(out, sprintf("%s__%s.rds", task$dataset_id, task$method))
if (file.exists(path)) quit(save = "no")
cat("karst", as.character(utils::packageVersion("karst")), find.package("karst"), "\n")
dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
nodes <- readRDS(file.path(config$data_dir, "nodes", paste0(task$dataset_id, ".rds")))
dataset$beta <- dataset$beta * task$multiple
work_dir <- file.path(tempdir(), paste0(task$dataset_id, "__", task$method))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(dataset$seed)
result <- stt_run_method(config$methods[[task$method]], dataset, nodes, work_dir)
result$dataset_id <- task$dataset_id; result$method <- task$method
result$beta_multiple <- task$multiple; result$beta_fitted <- dataset$beta
tmp <- paste0(path, ".tmp", Sys.getpid()); saveRDS(result, tmp); file.rename(tmp, path)
cat(task$dataset_id, task$method, "x", task$multiple, result$status, round(result$seconds), "s\n")
