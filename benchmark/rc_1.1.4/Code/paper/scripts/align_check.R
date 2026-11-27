# Checks of version B's building blocks (R/align.R), before its round:
#
#   Rscript scripts/align_check.R [N]
#
# 1. align_abm_forecast() against predict_evo(prediction_type = "ABM") on
#    alfak's own lists: the round's forecast of N ABM and N Salehi fits
#    (alfak, minobs 10, seed 1) made both ways must give the same frequencies
#    (largest difference 0).
# 2. The stable Kriging model rebuilt from bootstrap_res.Rds
#    (align_alfak_landscape()) against the one alfak saved: the alfak refits of
#    dev/realdata/extend (landscape_output = TRUE, landscape_data.Rds; read
#    only, DEV_EXTEND_ALFAK) of N Salehi trajectories, predicted on their
#    landscape's karyotypes and on 2000 karyotypes one step beyond it.
# Prints the results and writes align_checks.Rds into the version B folder
# (ALIGN_RESULTS, else BENCH_RESULTS/align_karst).

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
n <- if (length(args)) as.integer(args[1]) else 3L
base <- paths$results_base
set.seed(1)
out <- list()

## 1. sparse forecasts --------------------------------------------------------------------
abm_dirs <- Sys.glob(file.path(base, "fits", "abm", "*", "w4_ends", "alfak", "minobs_10"))
sal_dirs <- Sys.glob(file.path(base, "fits", "salehi", "*", "ends", "alfak", "minobs_10"))
pick <- function(d) d[sample.int(length(d), min(n, length(d)))]
ok_dir <- function(d) file.exists(file.path(d, "run.Rds")) && identical(readRDS(file.path(d, "run.Rds"))$status, "ok")
rows <- list()
lineages <- readRDS(paths$salehi_lineages)
meta <- utils::read.csv(paths$salehi_meta)
for (d in c(pick(Filter(ok_dir, abm_dirs)), pick(Filter(ok_dir, sal_dirs)))) {
  p <- rev(strsplit(d, "/", fixed = TRUE)[[1]])
  src <- p[5]; dataset <- p[4]; design <- p[3]
  lscape <- readRDS(file.path(d, "landscape.Rds"))
  if (src == "abm") {
    sim <- read_abm_sim(dataset, paths$abm_raw)
    s <- bench_settings$forecast_abm
    dense <- forecast_abm(lscape, sim$abm_output$x, bench_settings$abm$pm, 1)
    sx <- sim$abm_output$x
    times <- as.numeric(colnames(sx))
    start <- max(times[times < bench_settings$abm$last_time])
    at <- c(0, utils::head(times[times > start], s$horizon) - start)
    x0 <- sx[, times == start]
    names(x0) <- rownames(sx)
    x0 <- x0[names(x0) %in% lscape$k & x0 > 0]
    x0 <- x0 / sum(x0)
    sparse <- align_abm_forecast(lscape$k, lscape$mean, x0, at, bench_settings$abm$pm, s$abm_pop_size, s$abm_delta_t,
                                 s$abm_max_pop, s$abm_culling_survival, s$abm_record_interval, 1)
  } else {
    yi <- salehi_input(dataset, ends = TRUE, paths)
    sim_time <- salehi_forecast_time(dataset, lineages, meta)
    if (is.null(sim_time)) next
    s <- bench_settings$forecast_salehi
    x0 <- rep(0, nrow(lscape))
    names(x0) <- lscape$k
    last <- yi$x[rownames(yi$x) %in% lscape$k, ncol(yi$x), drop = FALSE]
    x0[rownames(last)] <- last[, 1]
    x0 <- x0 / sum(x0)
    at <- seq(0, sim_time, by = 1)
    dense <- alfakR::predict_evo(lscape = lscape, p = bench_settings$salehi$pm, times = at, x0 = x0,
                                 prediction_type = "ABM", abm_pop_size = s$abm_pop_size, abm_delta_t = s$abm_delta_t,
                                 abm_max_pop = bench_settings$salehi$nb, abm_record_interval = s$abm_record_interval,
                                 abm_culling_survival = s$abm_culling_survival, abm_seed = 1)
    sparse <- align_abm_forecast(lscape$k, lscape$mean, x0[x0 > 0], at, bench_settings$salehi$pm, s$abm_pop_size,
                                 s$abm_delta_t, bench_settings$salehi$nb, s$abm_culling_survival, s$abm_record_interval, 1)
  }
  wide <- align_densify(list(sparse), lscape$k)[[1]]
  dm <- as.matrix(dense[, -1])
  sm <- as.matrix(wide[, colnames(dm)])
  rows[[length(rows) + 1L]] <- data.frame(source = src, dataset = dataset, karyotypes = nrow(lscape),
                                          times = nrow(dm), max_abs_diff = max(abs(dm - sm)),
                                          same_times = isTRUE(all.equal(dense$time, wide$time)), stringsAsFactors = FALSE)
}
out$forecast <- do.call(rbind, rows)
cat("== sparse forecasts against predict_evo()\n")
print(out$forecast, row.names = FALSE)

## 2. Kriging rebuilt from bootstrap_res.Rds ------------------------------------------------
# rc freeze: this check compares with the dev study's saved alfak landscapes
# (outside the freeze); it runs only when DEV_EXTEND_ALFAK names that folder.
dev <- Sys.getenv("DEV_EXTEND_ALFAK")
if (!nzchar(dev)) { cat("DEV_EXTEND_ALFAK is not set: section 2 skipped\n"); dev <- tempfile("no_dev_") }
kd <- Sys.glob(file.path(dev, "*", "minobs_10"))
kd <- kd[file.exists(file.path(kd, "landscape_data.Rds")) & file.exists(file.path(kd, "bootstrap_res.Rds"))]
rows <- list()
loadNamespace("fields")
for (d in kd[sample.int(length(kd), min(n, length(kd)))]) {
  lsc <- readRDS(file.path(d, "landscape.Rds"))
  saved <- readRDS(file.path(d, "landscape_data.Rds"))$mean
  # alfak's list holds the frequent karyotypes' neighbours within two steps;
  # the probe adds 2000 karyotypes one step beyond it
  beyond <- setdiff(karst:::alfak_neighbours(lsc$k), lsc$k)
  probe <- c(lsc$k, beyond[sample.int(length(beyond), min(2000L, length(beyond)))])
  A <- align_alfak_landscape(d, unique(probe))
  km <- karst:::parse_karyotypes(unique(probe))
  storage.mode(km) <- "double"
  a <- as.numeric(stats::predict(saved, km))
  b <- as.numeric(stats::predict(A$krig, km))
  rows[[length(rows) + 1L]] <- data.frame(trajectory = basename(dirname(d)), points = length(a),
                                          max_abs_diff = max(abs(a - b)), sd = stats::sd(a),
                                          agreement_with_mean = A$agreement, stringsAsFactors = FALSE)
}
out$kriging <- do.call(rbind, rows)
cat("== stable Kriging rebuilt from bootstrap_res.Rds against alfak's saved model\n")
print(out$kriging, row.names = FALSE)
out_dir <- if (nzchar(Sys.getenv("ALIGN_RESULTS"))) Sys.getenv("ALIGN_RESULTS") else paths$results
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(out, file.path(out_dir, "align_checks.Rds"))
