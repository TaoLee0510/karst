# KARST's cross-validation scored against alfak's targets
# (docs/xv_target_6ae2a24.md).
#
#   Rscript scripts/xv_alfak_target.R [SOURCE ...]     (default: abm misseg salehi)
#
# Both methods hold out the same folds (each frequent karyotype with its
# one-step neighbours, each karyotype in the first fold that holds it) and
# report alfak's R2R, but they compare the held-out prediction with
# different targets:
#   alfak (alfakR dev2 f93b2d2, xval()): the Kriging prediction from the
#     other karyotypes against one bootstrap draw of the held-out karyotype's
#     first-stage fitness (mostly from its own counts, with their sampling
#     noise);
#   KARST (as_alfak(), sparse_joint_alfak_xval()): the posterior mode
#     without the fold's counts against its own estimate from all the data.
# Both xval.Rds files keep, per karyotype tested, the target
# ("observation") and the held-out prediction ("prediction") in xval_data.
# Here KARST's held-out predictions are compared with alfak's targets
# instead, on the karyotypes both tested, for every fit of a KARST method
# (karst, karst_pairwise, karst_state_space) that succeeded, against alfak's
# fit of the same input, design and minobs (an "all" design: alfak's fit of
# every time point). One row per KARST fit:
#   xv_own                 KARST's own score (its xval.Rds)
#   xv_alfak_target        R2R of alfak's targets and KARST's predictions
#   pearson_alfak_target   their Pearson correlation (scale-free: R2R also
#                          penalises a prediction on another scale than the
#                          targets)
#   alfak_xv_own           alfak's own score
#   alfak_xv_same          alfak's R2R on the same karyotypes
#   alfak_pearson_same     its Pearson correlation there
#   n_karst, n_alfak, n_common   karyotypes tested by each and by both
#   status                 "ok"; "alfak failed" (no alfak fit, or no
#                          targets); "fewer than 2 common karyotypes"
# xv_alfak_target is NA unless status is "ok", and the analyses then leave
# the fit out where the paper selects by the score (R/xv_source.R). R2R is
# alfak's R2R() (both series centred first, as karst's alfak_r2r()); the
# script checks that recomputing each fit's own score from its xval_data
# gives the stored R2R. Writes BENCH_RESULTS/tables/xv_alfak_target.Rds and
# stops first when tables/ksets.Rds shows a fit whose landscapes differ.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
sources <- commandArgs(TRUE)
if (!length(sources)) sources <- c("abm", "misseg", "salehi")
require_same_ksets(sources, paths)
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

r2r <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum(obs^2)
}
run_ok <- function(d) {
  p <- file.path(d, "run.Rds")
  file.exists(p) && identical(readRDS(p)$status, "ok") && file.exists(file.path(d, "xval.Rds"))
}
# Targets and predictions of one fit, one row per karyotype tested.
xval_pairs <- function(d) {
  x <- readRDS(file.path(d, "xval.Rds"))
  v <- x$xval_data
  v <- v[stats::complete.cases(v[c("observation", "prediction")]), c("k", "observation", "prediction")]
  if (anyDuplicated(v$k)) stop("duplicated karyotypes in ", file.path(d, "xval.Rds"))
  list(R2R = as.numeric(x$R2R), v = v)
}
own_check <- function(x) {
  if (nrow(x$v) < 2 || !is.finite(x$R2R)) return(NA_real_)
  abs(r2r(x$v$observation, x$v$prediction) - x$R2R)
}

one_fit <- function(d2, src) {
  p <- rev(strsplit(d2, "/", fixed = TRUE)[[1]])
  head <- data.frame(source = src, dataset = p[4], design = p[3], method = p[2], minobs = sub("minobs_", "", p[1]),
                     stringsAsFactors = FALSE)
  out <- cbind(head, xv_own = NA_real_, xv_alfak_target = NA_real_, pearson_alfak_target = NA_real_,
               alfak_xv_own = NA_real_, alfak_xv_same = NA_real_, alfak_pearson_same = NA_real_,
               n_karst = NA_integer_, n_alfak = NA_integer_, n_common = NA_integer_, status = NA_character_,
               check_karst = NA_real_, check_alfak = NA_real_, stringsAsFactors = FALSE)
  if (!run_ok(d2)) return(NULL)
  k <- xval_pairs(d2)
  out$xv_own <- k$R2R
  out$n_karst <- nrow(k$v)
  out$check_karst <- own_check(k)
  dR <- file.path(dirname(dirname(d2)), "alfak", basename(d2))
  if (!run_ok(dR)) {
    out$status <- "alfak failed"
    return(out)
  }
  a <- xval_pairs(dR)
  out$alfak_xv_own <- a$R2R
  out$n_alfak <- nrow(a$v)
  out$check_alfak <- own_check(a)
  common <- intersect(k$v$k, a$v$k)
  out$n_common <- length(common)
  if (!nrow(a$v)) {
    out$status <- "alfak failed"
    return(out)
  }
  if (length(common) < 2) {
    out$status <- "fewer than 2 common karyotypes"
    return(out)
  }
  obs <- a$v$observation[match(common, a$v$k)]
  pk <- k$v$prediction[match(common, k$v$k)]
  pa <- a$v$prediction[match(common, a$v$k)]
  out$xv_alfak_target <- r2r(obs, pk)
  out$pearson_alfak_target <- suppressWarnings(stats::cor(obs, pk))
  out$alfak_xv_same <- r2r(obs, pa)
  out$alfak_pearson_same <- suppressWarnings(stats::cor(obs, pa))
  out$status <- "ok"
  out
}

rows <- list()
for (src in sources) {
  dirs <- Sys.glob(file.path(paths$results, "fits", src, "*", "*", karst_xv_methods, "minobs_*"))
  cat(src, ":", length(dirs), "KARST fit folders\n")
  res <- parallel::mclapply(dirs, one_fit, src = src, mc.cores = cores, mc.preschedule = TRUE)
  err <- vapply(res, inherits, logical(1), "try-error")
  if (any(err)) stop(src, ": ", sum(err), " fits failed, e.g. ", dirs[err][1], ": ", res[err][[1]])
  rows[[src]] <- do.call(rbind, res)
}
rows <- do.call(rbind, rows)
rownames(rows) <- NULL
# The stored scores must be what alfak's R2R gives on the stored pairs.
worst <- max(c(rows$check_karst, rows$check_alfak), na.rm = TRUE)
cat("largest difference between a stored score and its recomputation:", signif(worst, 3), "\n")
if (worst > 1e-8) {
  print(utils::head(rows[pmax(rows$check_karst, rows$check_alfak, na.rm = TRUE) > 1e-8, ], 10))
  stop("a stored cross-validation score differs from alfak's R2R on its own xval_data")
}
rows$check_karst <- rows$check_alfak <- NULL
dir.create(file.path(paths$results, "tables"), recursive = TRUE, showWarnings = FALSE)
saveRDS(rows, file.path(paths$results, "tables", "xv_alfak_target.Rds"))
cat("\nfits by status:\n")
print(table(source = rows$source, method = rows$method, status = rows$status))
cat("\nmedian score and share > 0 (status ok):\n")
ok <- rows[rows$status == "ok", ]
print(do.call(rbind, lapply(split(ok, ok[c("source", "method")], drop = TRUE), function(x) data.frame(
  source = x$source[1], method = x$method[1], fits = nrow(x),
  own = round(stats::median(x$xv_own, na.rm = TRUE), 3), own_pos = round(mean(x$xv_own > 0, na.rm = TRUE), 3),
  alfak_target = round(stats::median(x$xv_alfak_target), 3), alfak_target_pos = round(mean(x$xv_alfak_target > 0), 3),
  alfak_own = round(stats::median(x$alfak_xv_own, na.rm = TRUE), 3),
  alfak_same = round(stats::median(x$alfak_xv_same), 3), alfak_same_pos = round(mean(x$alfak_xv_same > 0), 3),
  common_share = round(stats::median(x$n_common / x$n_alfak), 3)))), row.names = FALSE)
