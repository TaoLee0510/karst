#!/usr/bin/env Rscript
# Salehi forecasts against later passages (README.md). One array task: one fit
# (trajectory, method, minobs) of a two-time-point input ("ends" design: the
# first and the last passage of a trajectory). Its forecast is run from the
# trajectory's last passage, as the paper benchmark runs it
# (forecast_salehi_dir(): alfakR::predict_evo, ABM, 5 replicates), but up to
# the furthest of the parent's descendants within MAX_DEPTH levels of the
# metadata's parent links, and scored against every one of them. A level is
# a sampled passage, not one passage: adjacent levels are one passage (15
# days) apart in most PDX lines, but 5 to 21 passages (25 to 105 days) in
# SA906 and SA039. Each row has `depth` (levels), `passages_ahead` (passages,
# from the time point labels) and `days`; analyses stratify by days.
#
# Scores for each descendant passage d, on the karyotypes of the fit's
# landscape (as the paper's Fig. 2d axis) with the share of d's cells there:
#   ll_fc, ll_p0   mean log-likelihood per cell of d's counts under the
#                  forecast and under the parent's frequencies (persistence),
#                  each mixed with a uniform distribution over the landscape
#                  with weight EPS (1e-3 and 1e-2): a proper score, in which
#                  sampling noise does not favour a forecast that moves more.
#                  Persistence is the parent's raw frequencies; the forecast
#                  (an ABM run with missegregation) spreads some probability to
#                  one-step neighbours even with no selection, which alone can
#                  score above persistence (README.md).
#   ll_fc_s, ll_p0_s  the same on a common support: the union of the
#                  landscapes of every fit of this trajectory and minobs, both
#                  designs ("ends", "all"), under the task's support roots
#                  (column support_roots, ";"-separated; default the task's
#                  root). The uniform floor, d's counts and the persistence
#                  baseline are then the same for every fit of the trajectory,
#                  so fits of the two designs can be compared.
#   w_fc, w_p0     Wasserstein distance from the forecast and from the parent
#                  to d.
#   noise          Wasserstein distance between two random halves of d's
#                  cells (mean of 20 splits): the error no forecast can go
#                  below.
# Run from the paper benchmark's folder (sources R/), read only.
#   Rscript forecast_score.R TASKS.tsv TASK_ID OUT_DIR
args <- commandArgs(TRUE)
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
tasks <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)
task <- tasks[as.integer(args[[2]]), ]
# rc freeze: the task tables name the as-run dev results root; read the fits
# under DEV_RESULTS instead (Code/dev_repo/dev/hpc/env.sh)
rc_as_run <- "/share/lab_crd/taoli/Project/karst_dev/results/"
rc_root <- function(p) ifelse(startsWith(p, rc_as_run), file.path(Sys.getenv("DEV_RESULTS"), substring(p, nchar(rc_as_run) + 1L)), p)
task$root <- rc_root(task$root)
if (!is.null(task$support_roots) && !is.na(task$support_roots) && nzchar(task$support_roots))
  task$support_roots <- paste(rc_root(strsplit(task$support_roots, ";", fixed = TRUE)[[1]]), collapse = ";")
out_dir <- args[[3]]
max_depth <- 4L
# Design of the fit: "ends" (two time points) unless the task says "all"
# (every passage of the trajectory). The forecast starts from the last passage
# in both.
design <- if (is.null(task$design) || is.na(task$design)) "ends" else task$design
out_file <- file.path(out_dir, sprintf("%s__%s__%s__m%s.Rds", task$id, design, task$method, task$minobs))
if (file.exists(out_file)) quit(save = "no")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fit_dir <- file.path(task$root, task$id, design, task$method, paste0("minobs_", task$minobs))
done <- function(status, rows = NULL) {
  saveRDS(list(task = task, status = status, rows = rows), out_file)
  cat(task$id, task$method, task$minobs, status, "\n")
  quit(save = "no")
}
run <- file.path(fit_dir, "run.Rds")
if (!file.exists(run) || !identical(readRDS(run)$status, "ok") || !file.exists(file.path(fit_dir, "landscape.Rds"))) {
  done("no fit")
}
meta <- utils::read.csv(paths$salehi_meta, stringsAsFactors = FALSE)
lineages <- readRDS(paths$salehi_lineages)
inputs <- list.files(paths$salehi_inputs)
li <- lineages[[task$id]]
parent <- utils::tail(li$ids, 1)
tp <- function(uid) as.numeric(sub("X", "", meta$timepoint[match(uid, meta$uid)]))
days_per_p <- ifelse(meta$PDX_id[match(parent, meta$uid)] %in% c("SA906", "SA039"), 5, 15)
# Descendants of the parent, level by level of the metadata's parent links
# (depth = number of links, not of passages).
desc <- list()
frontier <- parent
for (depth in seq_len(max_depth)) {
  frontier <- meta$uid[!is.na(meta$parent) & meta$parent %in% frontier]
  if (!length(frontier)) break
  desc[[depth]] <- data.frame(uid = frontier, depth = depth, stringsAsFactors = FALSE)
}
if (!length(desc)) done("no descendants")
desc <- do.call(rbind, desc)
desc$passages_ahead <- tp(desc$uid) - tp(parent)
desc$days <- desc$passages_ahead * days_per_p
# A passage's counts: the last column of a trajectory that ends with it (file
# names begin with the passage's data set and time point, then "_l_").
counts_of <- function(uid) {
  stem <- paste(meta$datasetname[match(uid, meta$uid)], meta$timepoint[match(uid, meta$uid)], sep = "_")
  hit <- inputs[startsWith(inputs, paste0(stem, "_l_"))]
  if (!length(hit)) return(NULL)
  x <- readRDS(file.path(paths$salehi_inputs, hit[[1]]))$x
  v <- stats::setNames(x[, ncol(x)], rownames(x))
  v[v > 0]
}
yi <- salehi_input(task$id, ends = TRUE, paths)
lscape <- readRDS(file.path(fit_dir, "landscape.Rds"))
k <- lscape$k
x0 <- stats::setNames(numeric(length(k)), k)
last <- yi$x[rownames(yi$x) %in% k, ncol(yi$x)]
x0[names(last)] <- last
if (sum(x0) <= 0) done("parent outside landscape")
x0 <- x0 / sum(x0)
# Common support: the union of the landscapes of every fit of this trajectory
# and minobs, in both designs, under the support roots (always including this
# fit's own).
support_roots <- if (is.null(task$support_roots) || is.na(task$support_roots) || !nzchar(task$support_roots)) {
  task$root
} else strsplit(task$support_roots, ";", fixed = TRUE)[[1]]
support_globs <- as.vector(outer(support_roots, c("ends", "all"), function(r, d)
  file.path(r, task$id, d, "*", paste0("minobs_", task$minobs), "landscape.Rds")))
support_files <- unique(Sys.glob(support_globs))
ks <- sort(unique(c(k, unlist(lapply(support_files, function(f) as.character(readRDS(f)$k))))))
x0s <- stats::setNames(numeric(length(ks)), ks)
last_s <- yi$x[rownames(yi$x) %in% ks, ncol(yi$x)]
x0s[names(last_s)] <- last_s
x0s <- x0s / sum(x0s)
settings <- bench_settings$forecast_salehi
horizon <- max(desc$days)
started <- proc.time()[["elapsed"]]
fc <- lapply(seq_len(settings$reps), function(rep) {
  y <- alfakR::predict_evo(lscape = lscape, p = task$pm, times = seq(0, horizon, by = 1), x0 = x0,
                           prediction_type = "ABM", abm_pop_size = settings$abm_pop_size,
                           abm_delta_t = settings$abm_delta_t, abm_max_pop = task$nb,
                           abm_record_interval = settings$abm_record_interval,
                           abm_culling_survival = settings$abm_culling_survival, abm_seed = rep)
  m <- as.matrix(y[, -1, drop = FALSE])
  m <- m[, match(k, colnames(m)), drop = FALSE]
  m[is.na(m)] <- 0
  list(time = y$time, m = m)
})
seconds <- proc.time()[["elapsed"]] - started
s2m <- function(kk) do.call(rbind, lapply(strsplit(kk, ".", fixed = TRUE), as.numeric))
wd <- function(ka, ma, kb, mb) {
  a <- is.finite(ma) & ma > 0; b <- is.finite(mb) & mb > 0
  if (!any(a) || !any(b)) return(NA_real_)
  transport::wasserstein(transport::wpp(s2m(ka[a]), mass = ma[a] / sum(ma[a])),
                         transport::wpp(s2m(kb[b]), mass = mb[b] / sum(mb[b])))
}
halves <- function(v, reps = 20L) {
  if (sum(v) < 4) return(NA_real_)
  cells <- rep(names(v), v)
  mean(vapply(seq_len(reps), function(r) {
    set.seed(r)
    i <- sample.int(length(cells), length(cells) %/% 2L)
    a <- table(cells[i]); b <- table(cells[-i])
    wd(names(a), as.numeric(a), names(b), as.numeric(b))
  }, numeric(1)))
}
loglik <- function(y, q, eps) sum(y * log((1 - eps) * q + eps / length(q))) / sum(y)
rows <- lapply(seq_len(nrow(desc)), function(i) {
  y <- counts_of(desc$uid[i])
  if (is.null(y)) return(NULL)
  yl <- stats::setNames(numeric(length(k)), k)
  inl <- names(y)[names(y) %in% k]
  yl[inl] <- y[inl]
  if (sum(yl) <= 0) return(NULL)
  t_index <- vapply(fc, function(f) which.min(abs(f$time - desc$days[i])), integer(1))
  f <- Reduce(`+`, Map(function(f, j) f$m[j, ], fc, t_index)) / length(fc)
  f <- f / sum(f)
  # on the common support: the forecast is zero outside its own landscape
  ys <- stats::setNames(numeric(length(ks)), ks)
  ins <- names(y)[names(y) %in% ks]
  ys[ins] <- y[ins]
  fs <- stats::setNames(numeric(length(ks)), ks)
  fs[k] <- f
  data.frame(id = task$id, design = design, method = task$method, minobs = task$minobs, parent = parent,
             uid = desc$uid[i],
             depth = desc$depth[i], passages_ahead = desc$passages_ahead[i], days = desc$days[i],
             forecast_day = fc[[1]]$time[t_index[[1]]],
             treat_parent = meta$on_treatment[match(parent, meta$uid)],
             treat_desc = meta$on_treatment[match(desc$uid[i], meta$uid)],
             cells_desc = sum(y), cells_parent = sum(yi$x[, ncol(yi$x)]), coverage = sum(yl) / sum(y),
             landscape_size = length(k),
             ll_fc_e3 = loglik(yl, f, 1e-3), ll_p0_e3 = loglik(yl, x0, 1e-3),
             ll_fc_e2 = loglik(yl, f, 1e-2), ll_p0_e2 = loglik(yl, x0, 1e-2),
             support_size = length(ks), coverage_s = sum(ys) / sum(y),
             ll_fc_s_e3 = loglik(ys, fs, 1e-3), ll_p0_s_e3 = loglik(ys, x0s, 1e-3),
             ll_fc_s_e2 = loglik(ys, fs, 1e-2), ll_p0_s_e2 = loglik(ys, x0s, 1e-2),
             w_fc = wd(k, f, k, yl), w_p0 = wd(k, x0, k, yl), noise = halves(yl[yl > 0]),
             forecast_seconds = seconds, stringsAsFactors = FALSE)
})
done("ok", do.call(rbind, rows))
