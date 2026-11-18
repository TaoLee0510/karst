#!/usr/bin/env Rscript
# Carry one method's outputs of forecast_score.R over from an earlier run,
# under the method's current name, when its fits did not change: alfak's
# rows of the alfak2-era run (method alfakR there; alfak's fits of a KARST
# round are links to the same folders). Only the method name changes, in the
# file name, task$method and rows$method; everything else is kept as written.
#   Rscript reuse_rows.R FROM_DIR TO_DIR OLD_METHOD NEW_METHOD
# Refuses to overwrite an existing file of TO_DIR.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
from <- args[[1]]
to <- args[[2]]
old <- args[[3]]
new <- args[[4]]
dir.create(to, recursive = TRUE, showWarnings = FALSE)
files <- list.files(from, pattern = paste0("__", old, "__m[0-9]+[.]Rds$"), full.names = TRUE)
n <- 0L
for (f in files) {
  p <- readRDS(f)
  stopifnot(identical(p$task$method, old))
  # named as forecast_score.R names its outputs now (the first "ends" run
  # left the design out of the name)
  design <- if (is.null(p$task$design) || is.na(p$task$design)) "ends" else p$task$design
  out <- file.path(to, sprintf("%s__%s__%s__m%s.Rds", p$task$id, design, new, p$task$minobs))
  if (file.exists(out)) stop("Exists: ", out)
  p$task$method <- new
  if (!is.null(p$rows)) {
    p$rows$method <- rep(new, nrow(p$rows))
    if (is.null(p$rows$design)) p$rows$design <- rep(design, nrow(p$rows))
  }
  saveRDS(p, out)
  n <- n + 1L
}
cat("Carried over", n, "files of", old, "as", new, "from", from, "to", to, "\n")
