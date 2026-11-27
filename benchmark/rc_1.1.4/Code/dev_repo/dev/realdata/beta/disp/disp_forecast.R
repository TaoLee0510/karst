#!/usr/bin/env Rscript
# disp/PLAN.md: one forecast of one Salehi trajectory from its x1 fit (alfak
# output, minobs 10), with the prediction and scoring code of
# validation/sparse_joint/salehi_horizon/forecast_score.R, plus
#   variant "smooth"   (I1) parent frequencies with the Good-Turing unseen share
#                      spread over landscape karyotypes the parent sample lacks;
#   variant "drift1e3" / "drift1e4" (I2) ABM bottleneck of about 1e3 / 1e4 cells;
# and, for every later passage, the observed karyotypes' forecast probability,
# class (A-D) and distance to the parent sample (D1-D4).
#   Rscript dev/realdata/beta/disp/disp_forecast.R TASK_FILE TASK_ID
# Task columns: id, variant (base | smooth | drift1e3 | drift1e4), F (forecast rate multiple).
# Output: DEV_RESULTS/realdata/beta/disp/forecasts/<variant>_F<F>/<id>.Rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
rd_load_karst()
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
fit_root <- file.path(dev_results, "realdata", "beta", "salehi", "fits")
out_file <- file.path(dev_results, "realdata", "beta", "disp", "forecasts", paste0(task$variant, "_F", task$F), paste0(task$id, ".Rds"))
if (file.exists(out_file)) quit(save = "no")
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "beta", "salehi", "bench_unused"))
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
max_depth <- 4L; minobs <- 10L; design <- "ends"; pm <- 5e-5 * task$F; nb <- 1e7
done <- function(status, rows = NULL, kary = NULL, dist_mass = NULL, extra = NULL) {
  dev_save_rds(list(task = task, status = status, rows = rows, kary = kary, dist_mass = dist_mass, extra = extra), out_file)
  cat(task$id, task$variant, task$F, status, "\n")
  quit(save = "no")
}
fit_dir <- file.path(fit_root, task$id, design, "x1", paste0("minobs_", minobs))
run <- file.path(fit_dir, "run.Rds")
if (!file.exists(run) || !identical(readRDS(run)$status, "ok") || !file.exists(file.path(fit_dir, "landscape.Rds"))) done("no fit")

## ---- as forecast_score.R ----
meta <- utils::read.csv(paths$salehi_meta, stringsAsFactors = FALSE)
lineages <- readRDS(paths$salehi_lineages)
inputs <- list.files(paths$salehi_inputs)
li <- lineages[[task$id]]
parent <- utils::tail(li$ids, 1)
tp <- function(uid) as.numeric(sub("X", "", meta$timepoint[match(uid, meta$uid)]))
days_per_p <- ifelse(meta$PDX_id[match(parent, meta$uid)] %in% c("SA906", "SA039"), 5, 15)
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
support_globs <- as.vector(outer(fit_root, c("ends", "all"), function(r, d)
  file.path(r, task$id, d, "*", paste0("minobs_", minobs), "landscape.Rds")))
support_files <- unique(Sys.glob(support_globs))
ks <- sort(unique(c(k, unlist(lapply(support_files, function(f) as.character(readRDS(f)$k))))))
settings <- bench_settings$forecast_salehi
horizon <- max(desc$days)

## ---- parent sample, earlier passages, distances ----
s2m <- function(kk) do.call(rbind, lapply(strsplit(kk, ".", fixed = TRUE), as.numeric))
parent_all <- yi$x[, ncol(yi$x)]
parent_all <- parent_all[parent_all > 0]
parent_k <- names(parent_all)
n_parent <- sum(parent_all)
gt_unseen <- min(0.5, sum(parent_all == 1) / n_parent)
earlier <- list()
for (uid in utils::head(li$ids, -1)) { v <- counts_of(uid); if (!is.null(v)) earlier[[length(earlier) + 1L]] <- v }
earlier <- if (length(earlier)) tapply(unlist(earlier), names(unlist(earlier)), sum) else numeric(0)
parent_m <- s2m(parent_k)
min_dist <- function(kk) {
  if (!length(kk)) return(numeric(0))
  m <- s2m(kk)
  vapply(seq_len(nrow(m)), function(i) min(rowSums(abs(sweep(parent_m, 2, m[i, ])))), numeric(1))
}
k_dist <- stats::setNames(min_dist(k), k)
dist_class <- function(d) cut(d, c(-Inf, 0, 1, 2, Inf), labels = c("0", "1", "2", "3+"))

## ---- variant ----
abm_pop_size <- settings$abm_pop_size; abm_max_pop <- nb; abm_culling <- settings$abm_culling_survival
smooth_target <- NA_character_
if (task$variant == "smooth") {
  unseen <- k[!k %in% parent_k]
  w <- stats::setNames(as.numeric(earlier[unseen]), unseen); w[is.na(w)] <- 0
  if (sum(w) > 0) { g <- w / sum(w); smooth_target <- "earlier passages" } else {
    nb1 <- unseen[k_dist[unseen] == 1]
    g <- stats::setNames(numeric(length(unseen)), unseen)
    if (length(nb1)) { g[nb1] <- 1 / length(nb1); smooth_target <- "one-step neighbours" } else smooth_target <- "none"
  }
  if (sum(g) > 0) { x0 <- (1 - gt_unseen) * x0; x0[names(g)] <- x0[names(g)] + gt_unseen * g }
} else if (task$variant %in% c("drift1e3", "drift1e4")) {
  b <- if (task$variant == "drift1e3") 1e3 else 1e4
  abm_pop_size <- b; abm_max_pop <- 10 * b; abm_culling <- 0.1
} else if (task$variant != "base") stop("unknown variant ", task$variant)

## ---- forecast (as forecast_score.R) ----
started <- proc.time()[["elapsed"]]
fc <- lapply(seq_len(settings$reps), function(rep) {
  y <- alfakR::predict_evo(lscape = lscape, p = pm, times = seq(0, horizon, by = 1), x0 = x0,
                           prediction_type = "ABM", abm_pop_size = abm_pop_size,
                           abm_delta_t = settings$abm_delta_t, abm_max_pop = abm_max_pop,
                           abm_record_interval = settings$abm_record_interval,
                           abm_culling_survival = abm_culling, abm_seed = rep)
  m <- as.matrix(y[, -1, drop = FALSE])
  m <- m[, match(k, colnames(m)), drop = FALSE]
  m[is.na(m)] <- 0
  list(time = y$time, m = m)
})
seconds <- proc.time()[["elapsed"]] - started
loglik <- function(y, q, eps) sum(y * log((1 - eps) * q + eps / length(q))) / sum(y)

rows <- list(); kary <- list(); dmass <- list()
for (i in seq_len(nrow(desc))) {
  y <- counts_of(desc$uid[i])
  if (is.null(y)) next
  if (sum(y[names(y) %in% k]) <= 0) next
  t_index <- vapply(fc, function(f) which.min(abs(f$time - desc$days[i])), integer(1))
  f <- Reduce(`+`, Map(function(f, j) f$m[j, ], fc, t_index)) / length(fc)
  f <- f / sum(f)
  names(f) <- k   # positional, as in forecast_score.R (fs[k] <- f)
  ys <- stats::setNames(numeric(length(ks)), ks)
  ins <- names(y)[names(y) %in% ks]
  ys[ins] <- y[ins]
  fs <- stats::setNames(numeric(length(ks)), ks)
  fs[k] <- f
  eps <- 1e-2
  rows[[length(rows) + 1L]] <- data.frame(id = task$id, uid = desc$uid[i], depth = desc$depth[i], days = desc$days[i],
                                          cells_desc = sum(y), coverage_s = sum(ys) / sum(y),
                                          ll_fc_s_e2 = loglik(ys, fs, eps), n_parent = n_parent, gt_unseen = gt_unseen,
                                          smooth_target = smooth_target, forecast_seconds = seconds, stringsAsFactors = FALSE)
  # every observed karyotype of the passage (inside and outside the support)
  obs <- names(y)
  cls <- ifelse(obs %in% parent_k, "A", ifelse(obs %in% names(earlier), "B", NA))
  d <- min_dist(obs)
  cls[is.na(cls)] <- ifelse(d[is.na(cls)] == 1, "C", "D")
  q <- fs[obs]
  kary[[length(kary) + 1L]] <- data.frame(uid = desc$uid[i], days = desc$days[i], karyotype = obs, y = as.numeric(y),
                                          in_support = obs %in% ks, class = cls, dist = d, q = as.numeric(q),
                                          contrib = ifelse(obs %in% ks, as.numeric(y) * log((1 - eps) * as.numeric(q) + eps / length(ks)), NA),
                                          stringsAsFactors = FALSE)
  dm <- tapply(f, dist_class(k_dist[k]), sum)
  dmass[[length(dmass) + 1L]] <- data.frame(uid = desc$uid[i], days = desc$days[i], dist = names(dm), mass = as.numeric(dm),
                                            stringsAsFactors = FALSE)
}
if (!length(rows)) done("no scored passage")
done("ok", do.call(rbind, rows), do.call(rbind, kary), do.call(rbind, dmass),
     list(parent_k = parent_k, n_parent = n_parent, earlier_k = names(earlier), landscape_size = length(k)))
