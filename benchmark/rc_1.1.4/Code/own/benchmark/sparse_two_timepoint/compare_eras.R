#!/usr/bin/env Rscript
# Whether a rerun reproduces the fits of an earlier run: each fit of FITS (a
# fits/ folder of files <dataset>__<method>.rds) against the earlier run's fit
# of the same dataset and method (OLD_FITS; the alfak2 era named the methods
# alfak2_* and alfakR_*, which map to karst_* and alfak_*). Compares status
# and the predicted fitness and its sd per karyotype; prints a table by method
# and writes one row per fit to OUT.csv.
#   Rscript benchmark/sparse_two_timepoint/compare_eras.R FITS OLD_FITS OUT.csv
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
new_dir <- args[[1]]
old_dir <- args[[2]]
files <- list.files(new_dir, pattern = "__.*[.]rds$")
method <- sub("[.]rds$", "", sub("^.*__", "", files))
files <- files[startsWith(method, "karst_")]
old_name <- sub("__karst_", "__alfak2_", files, fixed = TRUE)
max_abs <- function(a, b) if (length(a) && length(a) == length(b)) max(abs(a - b), na.rm = TRUE) else NA_real_
rows <- do.call(rbind, lapply(seq_along(files), function(i) {
  n <- readRDS(file.path(new_dir, files[i]))
  o_path <- file.path(old_dir, old_name[i])
  base <- data.frame(file = files[i], method = sub("[.]rds$", "", sub("^.*__", "", files[i])),
                     status = n$status, old_status = NA_character_, karyotypes_equal = NA,
                     fitness_max_diff = NA_real_, sd_max_diff = NA_real_, seconds = n$seconds,
                     old_seconds = NA_real_, stringsAsFactors = FALSE)
  if (!file.exists(o_path)) return(base)
  o <- readRDS(o_path)
  base$old_status <- o$status
  base$old_seconds <- if (is.null(o$seconds)) NA_real_ else o$seconds
  if (!is.null(n$prediction) && !is.null(o$prediction)) {
    p <- n$prediction
    q <- o$prediction[match(p$karyotype, o$prediction$karyotype), ]
    base$karyotypes_equal <- setequal(p$karyotype, o$prediction$karyotype)
    base$fitness_max_diff <- max_abs(p$fitness, q$fitness)
    base$sd_max_diff <- if (is.null(p$fitness_sd) || all(is.na(p$fitness_sd))) NA_real_ else
      max_abs(p$fitness_sd, q$fitness_sd)
  }
  base
}))
utils::write.csv(rows, args[[3]], row.names = FALSE)
# A comparison that cannot be made (no earlier fit, a missing status or
# prediction) counts as not identical.
same_status <- (rows$status == rows$old_status) %in% TRUE
both_ok <- same_status & rows$status %in% "ok"
rows$identical <- same_status & (!both_ok | (rows$karyotypes_equal %in% TRUE & rows$fitness_max_diff %in% 0 &
                                               (is.na(rows$sd_max_diff) | rows$sd_max_diff %in% 0)))
rows$within_1e8 <- same_status & (!both_ok | (rows$karyotypes_equal %in% TRUE &
                                                (rows$fitness_max_diff < 1e-8) %in% TRUE))
summary <- do.call(rbind, lapply(split(rows, rows$method), function(d) data.frame(
  method = d$method[1], fits = nrow(d), with_old = sum(!is.na(d$old_status)), identical = sum(d$identical, na.rm = TRUE),
  within_1e8 = sum(d$within_1e8, na.rm = TRUE), largest_fitness_diff = suppressWarnings(max(d$fitness_max_diff, na.rm = TRUE)),
  time_ratio = stats::median(d$seconds / d$old_seconds, na.rm = TRUE), stringsAsFactors = FALSE)))
print(summary, row.names = FALSE)
cat("All fits:", nrow(rows), " identical:", sum(rows$identical, na.rm = TRUE), " within 1e-8:",
    sum(rows$within_1e8, na.rm = TRUE), "\n")
