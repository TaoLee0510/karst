# Which cross-validation score KARST's fits carry in the analyses, where the
# analyses write, and the check that both methods score the same karyotypes.
#
# XV_SOURCE unset or "own" (the default): every fit carries its own score,
# the R2R of its xval.Rds. For KARST (as_alfak(), sparse_joint_alfak_xval())
# that is its held-out predictions against its own estimates from all the
# data; for alfak (alfakR dev2 f93b2d2, xval()) the Kriging predictions
# against one bootstrap draw of each held-out karyotype's first-stage
# fitness.
# XV_SOURCE = "alfak_target": KARST's fits carry the R2R of their held-out
# predictions against alfak's targets of the same input, design and minobs,
# on the karyotypes both tested (tables/xv_alfak_target.Rds, written by
# scripts/xv_alfak_target.R). A KARST fit whose alfak fit failed has no
# score (NA), so the paper's rule leaves it out.
# XV_SOURCE = "own_target": against targets from each karyotype's own counts,
# resampled once as alfak's bootstrap resamples cells (karst 1.1.4,
# as_alfak(xval_target = "own"); tables/xv_own_target.Rds, written by
# scripts/xv_own_target.R); it does not depend on alfak. "own_point": the
# same targets without resampling (sensitivity analysis).
# alfak's and the paper's fits keep their own scores. The analyses then
# write to analysis_xv_<source>/ and figures_xv_<source>/ instead of
# analysis/ and figures/, so that nothing of the default run is replaced.

xv_source <- function() {
  s <- Sys.getenv("XV_SOURCE", "")
  if (s %in% c("", "own")) return("own")
  if (!s %in% names(xv_source_tables)) {
    stop("XV_SOURCE must be unset or one of \"own\", ", paste0("\"", names(xv_source_tables), "\"", collapse = ", "),
         ", not \"", s, "\"")
  }
  s
}
# Table and column of each source of KARST's scores.
xv_source_tables <- list(alfak_target = c(table = "xv_alfak_target", column = "xv_alfak_target"),
                         own_target = c(table = "xv_own_target", column = "xv_own_target"),
                         own_point = c(table = "xv_own_target", column = "xv_own_point"))
xv_suffix <- function() if (xv_source() == "own") "" else paste0("_xv_", xv_source())
analysis_root <- function(paths = bench_paths()) file.path(paths$results, paste0("analysis", xv_suffix()))
figures_root <- function(paths = bench_paths()) file.path(paths$results, paste0("figures", xv_suffix()))

# KARST's methods whose fits have a cross-validation score of their own.
karst_xv_methods <- c("karst", "karst_pairwise", "karst_state_space")

# A table of BENCH_RESULTS/tables/ (name without ".Rds"). With XV_SOURCE =
# "alfak_target", the cross-validation score of KARST's rows (xv in
# salehi_scores; Rxv, and pos_xv where present, in abm_dfl, abm_dfp and
# misseg_dfl) is replaced by the score against alfak's targets; other tables
# are returned as they are.
read_bench_table <- function(name, paths = bench_paths()) {
  d <- readRDS(file.path(paths$results, "tables", paste0(name, ".Rds")))
  if (xv_source() == "own") return(d)
  col <- switch(name, salehi_scores = "xv", abm_dfl = , abm_dfp = , misseg_dfl = "Rxv", NULL)
  if (is.null(col)) return(d)
  spec <- xv_source_tables[[xv_source()]]
  target <- readRDS(file.path(paths$results, "tables", paste0(spec[["table"]], ".Rds")))
  src <- switch(name, salehi_scores = "salehi", misseg_dfl = "misseg", "abm")
  target <- target[target$source == src, ]
  if (name == "salehi_scores") {
    key_d <- paste(d$fi, d$design, d$method, d$min_obs)
  } else {
    key_d <- paste(d$sim_id, d$design, d$method, d$minobs)
  }
  key_t <- paste(target$dataset, target$design, target$method, target$minobs)
  if (anyDuplicated(key_t)) stop(spec[["table"]], ".Rds has duplicated fits")
  karst <- d$method %in% karst_xv_methods
  # Every KARST fit with a score of its own must be in the target table
  # (with NA where alfak's fit failed); otherwise the table is stale.
  has_own <- karst & is.finite(d[[col]])
  missing <- unique(key_d[has_own & !key_d %in% key_t])
  if (length(missing)) {
    stop(length(missing), " KARST fits of ", name, " with a score are not in ", spec[["table"]], ".Rds, e.g. ",
         missing[1], ": rerun scripts/", spec[["table"]], ".R")
  }
  d[[col]][karst] <- target[[spec[["column"]]]][match(key_d[karst], key_t)]
  if ("pos_xv" %in% names(d)) d$pos_xv[karst] <- d[[col]][karst] > 0
  attr(d, "xv_source") <- xv_source()
  d
}

# Stop unless alfak's and KARST's landscapes of every input, design and
# minobs that both fitted hold the same karyotypes, frequent karyotypes and
# one-step neighbours (tables/ksets.Rds, scripts/check_ksets.R). The paired
# comparisons, the scoring of forecasts on the landscape's karyotypes and
# the cross-validation against alfak's targets assume it. KARST builds the
# sets from the input with alfak's rule (as_alfak()), not from alfak's
# outputs, so this checks the rule's implementation on every fit.
require_same_ksets <- function(sources, paths = bench_paths()) {
  path <- file.path(paths$results, "tables", "ksets.Rds")
  if (!file.exists(path)) stop("tables/ksets.Rds is missing: run scripts/check_ksets.R first")
  k <- readRDS(path)
  k <- k[k$source %in% sources, , drop = FALSE]
  if (!nrow(k)) stop("tables/ksets.Rds has no rows for ", paste(sources, collapse = ", "), ": rerun scripts/check_ksets.R")
  bad <- k[!(k$same_k & k$same_fq & k$same_nn), , drop = FALSE]
  if (nrow(bad)) {
    print(utils::head(bad, 20))
    stop(nrow(bad), " fits of ", paste(sources, collapse = ", "), " have landscapes whose karyotypes differ from alfak's")
  }
  invisible(nrow(k))
}
