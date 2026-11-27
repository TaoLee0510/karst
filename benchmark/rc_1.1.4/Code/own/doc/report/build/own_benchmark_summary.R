# Summary of KARST's own simulation benchmark (benchmark/sparse_two_timepoint),
# held-out split, core scenarios: error ratios against the null predictor,
# Spearman, coverage, direction accuracy of one-copy changes, the sampling
# designs, and each held-out landscape on its own (the replicates).
# The scoring functions are the benchmark's (sourced from the repository).
#
#   Rscript build/own_benchmark_summary.R RESULTS_DIR OUT.json    (from doc/report)
#
# RESULTS_DIR is the copy made by build/fetch_own.sh (data/own_<round>), or the round's
# folder on the HPC: with 500 held-out landscapes the direction accuracy (every fit
# scored again) is best computed there, with OWN_CORES workers (default 2).
# ROUND (or KARST_OWN_ROUND; default: the folder's name) is written into the output as `round`,
# which build/extract_data.R checks. Pooled values come with 95% intervals from 2000
# bootstrap draws of the landscapes of a wavelength (seed 1).

args <- commandArgs(TRUE)
if (!length(args) %in% 2:3) stop("usage: Rscript build/own_benchmark_summary.R RESULTS_DIR OUT.json [ROUND]")
res <- normalizePath(args[1])
# rc freeze: the inputs (datasets, nodes, index, far nodes) are in STT_DATA_DIR
dat <- normalizePath(Sys.getenv("STT_DATA_DIR", res))
cores <- as.integer(Sys.getenv("OWN_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "2")))
own_round_label <- if (length(args) >= 3) args[3] else Sys.getenv("KARST_OWN_ROUND", basename(res))
out_file <- args[2]
repo <- normalizePath("../..")
old_wd <- setwd(repo)
for (part in c("landscape.R", "config.R", "score.R")) source(file.path("benchmark/sparse_two_timepoint", part))
config <- stt_config()
setwd(old_wd)

methods <- c(default = "karst_sparse_joint_balance", alfak5 = "alfak_minobs5", alfak10 = "alfak_minobs10",
             alfak20 = "alfak_minobs20", null = "zero")
fitted <- setdiff(names(methods), "null")   # the methods drawn, KARST and alfak at each minobs
# every wavelength class of the configuration, smallest (most rugged) first
waves <- config$landscapes$wavelength_class[order(config$landscapes$wavelength)]
waves <- unique(waves)
designs <- c("very_sparse_50_50_180d", "sparse_200_200_180d", "sparse_200_200_360d", "typical_500_700_360d",
             "long_500_700_1080d", "bottleneck_500_700_360d", "unbalanced_100_3000_360d", "unbalanced_3000_100_360d")
cols <- c("rmse", "rmse_zero", "spearman", "coverage95")
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

near <- read.csv(file.path(res, "metrics_heldout.csv"), stringsAsFactors = FALSE)
far <- read.csv(file.path(res, "sweeps", "far", "metrics_far_heldout.csv"), stringsAsFactors = FALSE)
keep <- c("dataset_id", "method", "landscape_id", "wavelength_class", "scenario", "role", "scoreboard", "status",
          "stratum", "seconds", cols)
m <- rbind(near[near$stratum %in% c("d0", "d1", "d2"), keep], far[far$stratum %in% c("d3", "d4"), keep])
m <- m[m$role == "core" & m$scoreboard & m$method %in% methods, ]

# mean over designs, then over the landscapes of a wavelength (as the benchmark reports it)
cell <- aggregate(m[cols], m[c("method", "wavelength_class", "landscape_id", "scenario", "stratum")], mean_na)
land <- aggregate(cell[cols], cell[c("method", "wavelength_class", "landscape_id", "stratum")], mean_na)
wl <- aggregate(land[cols], land[c("method", "wavelength_class", "stratum")], mean_na)
wl$ratio <- wl$rmse / wl$rmse_zero
ds <- aggregate(m[m$stratum %in% c("d0", "d1", "d2"), cols], m[m$stratum %in% c("d0", "d1", "d2"), c("method", "wavelength_class", "scenario", "stratum")], mean_na)
ds$ratio <- ds$rmse / ds$rmse_zero

pick <- function(tab, method, wave, stratum, column) {
  v <- tab[tab$method == method & tab$wavelength_class == wave & tab$stratum == stratum, column]
  if (length(v) && is.finite(v)) round(v, 3) else NA_real_
}
series <- function(tab, column) lapply(methods, function(mm) setNames(lapply(waves, function(w)
  vapply(paste0("d", 0:4), function(s) pick(tab, mm, w, s, column), 0)), waves))
count_below <- function(method) {
  x <- ds[ds$method == method, ]
  list(below = sum(x$ratio < 1, na.rm = TRUE), total = sum(is.finite(x$ratio)),
       min = round(min(x$ratio, na.rm = TRUE), 2), max = round(max(x$ratio, na.rm = TRUE), 2))
}
design_tab <- function(method) do.call(rbind, lapply(designs, function(s) data.frame(scenario = s,
  t(vapply(waves, function(w) pick(ds[ds$scenario == s, ], method, w, "d1", "ratio"), 0)), check.names = FALSE)))
fits <- unique(near[near$role == "core" & near$scoreboard & near$method %in% methods, c("dataset_id", "method", "status", "seconds", "wavelength_class")])
status <- do.call(rbind, lapply(split(fits, list(fits$method, fits$wavelength_class), drop = TRUE), function(d) data.frame(
  method = d$method[1], wave = d$wavelength_class[1], fits = nrow(d), failed = sum(d$status != "ok"),
  seconds_median = stats::median(d$seconds, na.rm = TRUE), seconds_total = sum(d$seconds, na.rm = TRUE))))

# cells per observed karyotype (the data, not a fit): median over the landscapes of a wavelength
sp <- near[near$method == "zero" & near$stratum == "d0" & near$role == "core" & near$scoreboard, ]
sp$cells_per_karyotype <- (sp$depth_t0 + sp$depth_t1) / sp$n_karyotypes
sparsity <- aggregate(cells_per_karyotype ~ scenario + wavelength_class, sp, stats::median)
cells <- do.call(rbind, lapply(designs, function(s) data.frame(scenario = s, t(vapply(waves, function(w) {
  v <- sparsity$cells_per_karyotype[sparsity$scenario == s & sparsity$wavelength_class == w]
  if (length(v) == 1) round(v, 1) else NA_real_
}, 0)), check.names = FALSE)))

# direction accuracy of one-copy changes: the share of neighbouring pairs (one copy of
# one chromosome apart) whose predicted fitness change has the true sign; ties count 1/2
direction <- function(nodes, truth, prediction, keep_distances) {
  table <- stt_score_predictions(nodes, truth, prediction, config$scoring$reference_min_count)$nodes
  cn <- stt_parse_karyotypes(table$karyotype)
  index_of <- stats::setNames(seq_len(nrow(table)), table$karyotype)
  edges <- do.call(rbind, lapply(seq_len(ncol(cn)), function(chr) {
    moved <- cn
    moved[, chr] <- moved[, chr] + 1
    child <- index_of[stt_format_karyotypes(moved)]
    keep <- which(!is.na(child))
    cbind(keep, unname(child[keep]))
  }))
  d <- pmax(table$distance[edges[, 1]], table$distance[edges[, 2]])
  dt <- table$truth_relative[edges[, 2]] - table$truth_relative[edges[, 1]]
  dp <- table$pred_relative[edges[, 2]] - table$pred_relative[edges[, 1]]
  dp[abs(dp) < 1e-12] <- 0
  score <- ifelse(dp == 0, 0.5, as.numeric(sign(dp) == sign(dt)))
  ok <- dt != 0
  do.call(rbind, lapply(keep_distances, function(k) {
    use <- ok & d == k
    data.frame(stratum = paste0("d", k), n_pairs = sum(use), direction = if (any(use)) mean(score[use]) else NA_real_)
  }))
}
get_pred <- function(path) {
  if (!file.exists(path)) return(data.frame(karyotype = character(), fitness = numeric()))
  p <- readRDS(path)$prediction
  if (is.null(p)) data.frame(karyotype = character(), fitness = numeric()) else p
}
index <- read.csv(file.path(dat, "dataset_index_heldout.csv"), stringsAsFactors = FALSE)
index <- index[index$scoreboard & index$role == "core", ]
dir_rows <- parallel::mclapply(index$dataset_id, function(id) {
  dataset <- readRDS(file.path(dat, "datasets", paste0(id, ".rds")))
  nn <- readRDS(file.path(dat, "nodes", paste0(id, ".rds")))
  ff <- readRDS(file.path(dat, "far", "nodes", paste0(id, ".rds")))
  tn <- stt_true_growth_rate(nn$karyotype, dataset$truth)
  tf <- stt_true_growth_rate(ff$karyotype, dataset$truth)
  do.call(rbind, lapply(fitted, function(mm) {
    a <- direction(nn, tn, get_pred(file.path(res, "fits", sprintf("%s__%s.rds", id, methods[[mm]]))), 0:2)
    b <- direction(ff, tf, get_pred(file.path(res, "sweeps", "far", "fits", sprintf("%s__%s.rds", id, methods[[mm]]))), 3:4)
    cbind(dataset_id = id, method = mm, rbind(a, b))
  }))
}, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(dir_rows, inherits, logical(1), "try-error")
if (any(bad)) stop("direction failed for: ", paste(index$dataset_id[bad], collapse = ", "))
x <- merge(do.call(rbind, dir_rows), index[c("dataset_id", "landscape_id", "wavelength_class", "scenario")])
dl <- aggregate(direction ~ method + wavelength_class + landscape_id + stratum, x, mean_na)
dw <- aggregate(direction ~ method + wavelength_class + stratum, dl, mean_na)
dir_series <- lapply(stats::setNames(fitted, fitted), function(mm) setNames(lapply(waves, function(w)
  vapply(paste0("d", 0:4), function(s) { v <- dw$direction[dw$method == mm & dw$wavelength_class == w & dw$stratum == s]
    if (length(v)) round(v, 3) else NA_real_ }, 0)), waves))

# Replicates: each held-out landscape (one dataset per sampling design) on its own.
# Error ratio per landscape = its mean RMSE / its mean null RMSE over the designs; per
# landscape and design = that dataset's ratio. The pooled ratios above are the ratio of
# the means over landscapes, a weighted mean of these, so they lie within their range.
land$ratio <- land$rmse / land$rmse_zero
ratio_land <- land[land$method %in% methods, c("method", "wavelength_class", "landscape_id", "stratum", "ratio")]
ratio_land$method <- names(methods)[match(ratio_land$method, methods)]
dset <- cell[cell$stratum == "d1" & cell$method %in% methods[fitted], ]
design_land <- data.frame(method = names(methods)[match(dset$method, methods)], wavelength_class = dset$wavelength_class,
                          landscape_id = dset$landscape_id, scenario = dset$scenario, ratio = dset$rmse / dset$rmse_zero)
cell$ratio <- cell$rmse / cell$rmse_zero
per_dataset <- lapply(stats::setNames(fitted, fitted), function(mm) {
  x <- cell[cell$method == methods[[mm]] & is.finite(cell$ratio), ]
  list(below = sum(x$ratio < 1), total = nrow(x))
})

# 95% intervals of the pooled values: the landscapes of a wavelength drawn with
# replacement (the same draws for every method, stratum and measure)
set.seed(1)
n_boot <- 2000L
boot_ix <- lapply(stats::setNames(waves, waves), function(w) {
  ids <- sort(unique(land$landscape_id[land$wavelength_class == w]))
  lapply(seq_len(n_boot), function(b) sample(ids, length(ids), replace = TRUE))
})
boot_ci <- function(tab, value) {   # tab: method, wavelength_class, landscape_id, stratum, value columns
  lapply(methods_ci, function(mm) setNames(lapply(waves, function(w) lapply(paste0("d", 0:4), function(st) {
    x <- tab[tab$method == mm & tab$wavelength_class == w & tab$stratum == st, ]
    if (!nrow(x)) return(c(NA_real_, NA_real_))
    v <- vapply(boot_ix[[w]], function(ids) value(x[match(ids, x$landscape_id), , drop = FALSE]), 0)
    round(unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)), 3)
  })), waves))
}
methods_ci <- methods[c(fitted, "null")]
ratio_ci <- boot_ci(land, function(x) mean(x$rmse, na.rm = TRUE) / mean(x$rmse_zero, na.rm = TRUE))
names(ratio_ci) <- c(fitted, "null")
methods_ci <- stats::setNames(fitted, fitted)
direction_ci <- boot_ci(dl, function(x) mean(x$direction, na.rm = TRUE))

out <- list(results = res, round = own_round_label, waves = waves,
            ratio_ci = ratio_ci, direction_ci = direction_ci, n_boot = n_boot,
            ratio_land = ratio_land, direction_land = dl, design_land = design_land, per_dataset = per_dataset,
            ratio = series(wl, "ratio")[c(fitted, "null")],
            spearman = series(wl, "spearman")[fitted],
            coverage = series(wl, "coverage95")[c("default")],
            direction = dir_series,
            designs = lapply(stats::setNames(fitted, fitted), function(mm) design_tab(methods[[mm]])),
            cells = cells,
            below_null = lapply(stats::setNames(fitted, fitted), function(mm) count_below(methods[[mm]])),
            status = status)
jsonlite::write_json(out, out_file, auto_unbox = TRUE, pretty = TRUE, digits = NA, dataframe = "rows")
cat("wrote", out_file, "\n")
