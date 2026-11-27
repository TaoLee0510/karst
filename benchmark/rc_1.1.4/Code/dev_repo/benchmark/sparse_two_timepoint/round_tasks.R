#!/usr/bin/env Rscript
# Task table of the fits a round still needs: one line per dataset and method
# (columns split, dataset_id, method), for run_task.R (hpc/submit.sh dev) or for
# the distance 3-4 sweep (validation/sparse_joint/far_distance.R; hpc/submit.sh
# far). The datasets are those of the scoreboard patient of the named split,
# listed from the configuration (config.R), so a table can be written before
# the datasets are generated. A dataset and method whose fit file already exists
# in the round's folder is left out.
#
#   Rscript benchmark/sparse_two_timepoint/round_tasks.R --split=heldout --methods=a,b OUT.tsv
#       [--roles=core,stress] [--sweep=main|far] [--landscapes=FROM-TO] [--scenarios=REGEX]
#       [--exclude-scenarios=REGEX]
#
# --roles: the scenario roles (default core). --sweep=far lists the fits of the
# distance 3-4 sweep (core scenarios only; its fits are in sweeps/far/fits).
# --landscapes=16-500 keeps the landscapes with those numbers (L16-L500).
# --scenarios / --exclude-scenarios keep or drop scenarios by name, e.g. to put
# the largest datasets in a task table of their own with more memory.
# Run from the repository root; STT_OUTPUT_DIR names the round's folder.
script_dir <- file.path("benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R")) source(file.path(script_dir, part))
args <- commandArgs(TRUE)
option <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else default
}
out_file <- args[!startsWith(args, "--")]
if (length(out_file) != 1L) stop("give one output file", call. = FALSE)
config <- stt_config()
split <- option("split")
methods <- strsplit(option("methods", ""), ",", fixed = TRUE)[[1L]]
if (is.null(split) || !length(methods)) stop("--split and --methods are required", call. = FALSE)
unknown <- setdiff(methods, names(config$methods))
if (length(unknown)) stop("unknown methods: ", paste(unknown, collapse = ", "), call. = FALSE)
sweep <- option("sweep", "main")
roles <- strsplit(option("roles", "core"), ",", fixed = TRUE)[[1L]]
if (sweep == "far" && !identical(roles, "core")) stop("the distance 3-4 sweep has core scenarios only", call. = FALSE)

landscapes <- config$landscapes[config$landscapes$split == split, ]
range <- option("landscapes")
if (!is.null(range)) {
  bounds <- as.integer(strsplit(range, "-", fixed = TRUE)[[1L]])
  number <- as.integer(sub("^.*_L", "", landscapes$landscape_id))
  landscapes <- landscapes[number >= bounds[1L] & number <= bounds[length(bounds)], ]
}
scenarios <- config$scenarios$scenario[config$scenarios$role %in% roles]
keep <- option("scenarios")
if (!is.null(keep)) scenarios <- grep(keep, scenarios, value = TRUE)
drop <- option("exclude-scenarios")
if (!is.null(drop)) scenarios <- grep(drop, scenarios, value = TRUE, invert = TRUE)
patient <- config$scoreboard_patients
if (length(patient) != 1L) stop("expected one scoreboard patient", call. = FALSE)

grid <- expand.grid(method = methods, scenario = scenarios, landscape_id = landscapes$landscape_id,
                    stringsAsFactors = FALSE)
grid$dataset_id <- sprintf("%s_P%02d_%s", grid$landscape_id, patient, grid$scenario)
fits <- if (sweep == "far") file.path(config$output_dir, "sweeps", "far", "fits") else file.path(config$output_dir, "fits")
done <- file.exists(file.path(fits, sprintf("%s__%s.rds", grid$dataset_id, grid$method)))
tasks <- data.frame(split = split, dataset_id = grid$dataset_id[!done], method = grid$method[!done],
                    stringsAsFactors = FALSE)
utils::write.table(tasks, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d tasks (%d landscapes x %d scenarios x %d methods, %d fits already there) in %s\n",
            sweep, nrow(tasks), nrow(landscapes), length(scenarios), length(methods), sum(done), config$output_dir))
