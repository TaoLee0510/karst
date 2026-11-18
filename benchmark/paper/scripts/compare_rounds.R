# Compare karst's fits of this round with those of an earlier round, on the
# inputs both fitted with a design matching PATTERN (default: two time points,
# designs ending in "ends"). A commit that leaves those fits unchanged should
# give the same status, landscape (estimate and sd of every karyotype; equal
# within `tol`, since another node may round differently by ~1e-16),
# cross-validation and, for the ABM data, forecasts (identical) and scores
# (Spearman's rho may differ in the 5th decimal: estimates that differ by
# ~1e-16 break ties differently). The bootstrap draws of as_alfak() are not
# compared: they are random.
#
#   [COMPARE_METHOD=karst_pairwise] Rscript scripts/compare_rounds.R OLD_RESULTS [PATTERN]
#
# COMPARE_METHOD: the method folder compared (default karst). Fit folders
# that are links (fits taken over from another round) are skipped.
#
# Writes to BENCH_RESULTS/analysis/compare_rounds/:
#   <old>[_<method>]_fits.tsv     one row per fit (source, input, design, minobs)
#   <old>[_<method>]_summary.tsv  per source, design and the initial state
#                      the new fit used (827c4b5 on: "balance" or "free";
#                      "-" before): fits, same status, equal landscape,
#                      cross-validation and scores, largest difference of a
#                      fitness estimate

for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: compare_rounds.R OLD_RESULTS [PATTERN]")
old_root <- normalizePath(args[[1]])
pattern <- if (length(args) > 1) args[[2]] else "ends$"
method <- Sys.getenv("COMPARE_METHOD", "karst")
# COMPARE_OLD_METHOD: the same method's folder name in OLD_RESULTS, when it
# differs (rounds of the alfak2 era named karst's folders alfak2, alfak2_pairwise
# and alfak2_state_space; default: COMPARE_METHOD).
old_method <- Sys.getenv("COMPARE_OLD_METHOD", method)
old_rel <- function(r) sub(paste0("/", method, "/"), paste0("/", old_method, "/"), r, fixed = TRUE)
paths <- bench_paths()
out <- file.path(paths$results, "analysis", "compare_rounds")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

new_dirs <- Sys.glob(file.path(paths$results, "fits", "*", "*", "*", method, "minobs_*"))
linked <- vapply(dirname(new_dirs), function(d) {
  l <- Sys.readlink(d)
  !is.na(l) && nzchar(l)
}, logical(1))
new_dirs <- new_dirs[!linked]
rel <- substring(new_dirs, nchar(file.path(paths$results, "fits")) + 2L)
parts <- do.call(rbind, strsplit(rel, "/", fixed = TRUE))
keep <- grepl(pattern, parts[, 3])
rel <- rel[keep]
parts <- parts[keep, , drop = FALSE]
cat(method, "fits of this round matching", pattern, ":", length(rel), " (links skipped:", sum(linked), ")\n")

read_or_null <- function(p) if (file.exists(p)) readRDS(p) else NULL
same <- function(a, b) if (is.null(a) || is.null(b)) NA else isTRUE(all.equal(a, b, tolerance = 0))
tol <- 1e-12
max_diff <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA_real_)
  a <- unlist(a)
  b <- unlist(b)
  if (length(a) != length(b)) return(Inf)
  suppressWarnings(max(abs(a - b), na.rm = TRUE))
}
forecasts_same <- function(o, n) {
  fo <- list.files(file.path(o, "abm_preds"))
  fn <- list.files(file.path(n, "abm_preds"))
  if (!length(fo) || !length(fn)) return(NA)
  identical(fo, fn) && all(vapply(fo, function(x) identical(readRDS(file.path(o, "abm_preds", x)),
                                                            readRDS(file.path(n, "abm_preds", x))), logical(1)))
}
compare_one <- function(i) {
  o <- file.path(old_root, "fits", old_rel(rel[i]))
  n <- file.path(paths$results, "fits", rel[i])
  ro <- read_or_null(file.path(o, "run.Rds"))
  rn <- read_or_null(file.path(n, "run.Rds"))
  lo <- read_or_null(file.path(o, "landscape.Rds"))
  ln <- read_or_null(file.path(n, "landscape.Rds"))
  diff_mean <- NA_real_
  diff_sd <- NA_real_
  if (!is.null(lo) && !is.null(ln)) {
    m <- match(lo$k, ln$k)
    if (!anyNA(m) && nrow(lo) == nrow(ln)) {
      diff_mean <- max(abs(lo$mean - ln$mean[m]))
      diff_sd <- max(abs(lo$sd - ln$sd[m]))
    }
  }
  xo <- read_or_null(file.path(o, "xval.Rds"))
  xn <- read_or_null(file.path(n, "xval.Rds"))
  so <- read_or_null(file.path(o, "score.Rds"))
  sn <- read_or_null(file.path(n, "score.Rds"))
  data.frame(source = parts[i, 1], input = parts[i, 2], design = parts[i, 3],
             minobs = sub("minobs_", "", parts[i, 5]),
             status_old = if (is.null(ro)) "missing" else ro$status,
             status_new = if (is.null(rn)) "missing" else rn$status,
             landscape_exact = same(lo, ln),
             landscape_equal = !is.na(diff_mean) && diff_mean <= tol && diff_sd <= tol,
             landscape_rows = if (is.null(ln)) NA_integer_ else nrow(ln),
             max_diff_mean = diff_mean, max_diff_sd = diff_sd,
             xval_identical = same(xo$predictions, xn$predictions),
             r2r_old = if (is.null(xo)) NA_real_ else xo$R2R, r2r_new = if (is.null(xn)) NA_real_ else xn$R2R,
             forecasts_identical = forecasts_same(o, n),
             score_max_diff = max_diff(so[-1], sn[-1]),
             history_old = first_or_na(ro$history), history_new = first_or_na(rn$history),
             state_new = if (is.null(rn$initial_state)) "-" else as.character(rn$initial_state)[[1L]],
             stringsAsFactors = FALSE)
}
rows <- if (cores > 1L) parallel::mclapply(seq_along(rel), compare_one, mc.cores = cores) else
  lapply(seq_along(rel), compare_one)
fits <- do.call(rbind, rows)
name <- basename(old_root)
if (method != "karst") name <- paste0(name, "_", method)
utils::write.table(fits, file.path(out, paste0(name, "_fits.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

count <- function(x) sum(x %in% TRUE)
summary <- do.call(rbind, lapply(split(fits, list(fits$source, fits$design, fits$state_new), drop = TRUE), function(d) {
  ok <- d$status_old == "ok" & d$status_new == "ok"
  data.frame(source = d$source[1], design = d$design[1], state_new = d$state_new[1], fits = nrow(d),
             same_status = sum(d$status_old == d$status_new), both_ok = sum(ok),
             landscape_exact = count(d$landscape_exact[ok]),
             landscape_equal = count(d$landscape_equal[ok]),
             xval_identical = count(d$xval_identical[ok]),
             scored = sum(!is.na(d$forecasts_identical[ok])),
             forecasts_identical = count(d$forecasts_identical[ok]),
             score_max_diff = if (all(is.na(d$score_max_diff[ok]))) NA_real_ else
               max(d$score_max_diff[ok], na.rm = TRUE),
             same_history = sum(d$history_old[ok] == d$history_new[ok], na.rm = TRUE),
             max_diff_mean = if (any(ok)) suppressWarnings(max(d$max_diff_mean[ok], na.rm = TRUE)) else NA_real_,
             max_diff_sd = if (any(ok)) suppressWarnings(max(d$max_diff_sd[ok], na.rm = TRUE)) else NA_real_,
             stringsAsFactors = FALSE)
}))
utils::write.table(summary, file.path(out, paste0(name, "_summary.tsv")), sep = "\t", quote = FALSE,
                   row.names = FALSE)
print(summary, row.names = FALSE)
changed <- fits[fits$status_old != fits$status_new, c("source", "input", "design", "minobs", "status_old",
                                                     "status_new")]
if (nrow(changed)) {
  cat("Fits whose status changed:\n")
  print(utils::head(changed, 30), row.names = FALSE)
}
