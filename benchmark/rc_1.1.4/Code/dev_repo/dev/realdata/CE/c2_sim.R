# C2 of dev/realdata/PLAN.md (phase C): simulations matched to Salehi.
#
# One trajectory per landscape and setting: 22 chromosomes, a diploid founder,
# the GRF of the development benchmark (wavelengths 0.8 and 1.2, 6 landscapes
# each), the C2 ABM (dev/abm/abm_c2.cpp) with a passage every 15 days that
# thins the population to `bottleneck` cells, a cap of 2e5 cells between
# passages, and a burn-in of `burn_days` before the first sample. Samples of
# 300 and 450 cells are drawn just before the passages at t0 and t0 + 15, 30,
# 45 or 60 days (four datasets per trajectory); the true populations at every
# passage up to t0 + 120 are kept for checks further ahead.
#
# Truth: the ABM's expected per-day growth rate of a karyotype, net of the
# divisions that leave it by missegregation at the simulated rate:
#   log(1 + f (1 - 2 p_any) dt) / dt,  p_any = 1 - (1 - beta_k)^copies,
# beta_k = beta_sim (x beta_factor at ploidy 4 with WGD). Fits use beta 5e-5,
# as the Salehi benchmark, whatever the simulated rate.
#
# Sourced after dev/R/common.R. Results under DEV_RESULTS/realdata/CE/c2/:
#   grid/<setting>/{truth,datasets,nodes}/   tuning rounds
#   regimes/<regime>/{datasets,nodes}/        accepted datasets (H, M, W)
#   fits/<method>/<regime or setting>/<dataset_id>.rds

c2_root <- file.path(dev_results, "realdata", "CE", "c2")
c2_abm <- list(delta_t = 0.1, passage_days = 15, max_pop = 2e5, min_cn = 1L)
c2_depth <- c(300, 450)
c2_intervals <- c(15, 30, 45, 60)
c2_ahead_days <- seq(0, 120, by = 15)
c2_fit_beta <- 5e-5

c2_landscapes <- data.frame(landscape_id = sprintf("c2_L%02d", 1:12), wavelength = rep(c(0.8, 1.2), each = 6L),
                            landscape_seed = 64000L + 1:12, stringsAsFactors = FALSE)

c2_bench_source <- function() {
  if (exists("stt_score_predictions", mode = "function")) return(invisible())
  for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R")) {
    source(file.path(dev_repo, "benchmark", "sparse_two_timepoint", part))
  }
}

c2_load_abm <- function() {
  if (exists("c2_abm_cpp", mode = "function")) return(invisible())
  cache <- file.path(tempdir(), "c2_cpp")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(dev_repo, "dev", "abm", "abm_c2.cpp"), cacheDir = cache, env = globalenv())
  invisible()
}

# Settings of the tuning rounds: one row each. W rows take the accepted M
# setting's rate, bottleneck and burn-in (filled in when round W is written).
c2_settings_path <- function() file.path(dev_repo, "dev", "realdata", "CE", "tasks", "c2_settings.tsv")
c2_settings <- function() utils::read.delim(c2_settings_path(), stringsAsFactors = FALSE)
c2_setting <- function(name) {
  s <- c2_settings()
  row <- s[s$setting == name, , drop = FALSE]
  if (nrow(row) != 1L) stop("No C2 setting ", name, call. = FALSE)
  as.list(row)
}

c2_dir <- function(key) {
  regime <- file.path(c2_root, "regimes", key)
  if (dir.exists(regime)) regime else file.path(c2_root, "grid", key)
}

c2_truth_path <- function(setting, landscape_id) file.path(c2_root, "grid", setting, "truth", paste0(landscape_id, ".rds"))

c2_parse <- function(labels) {
  k <- do.call(rbind, lapply(strsplit(labels, ".", fixed = TRUE), as.integer))
  storage.mode(k) <- "integer"
  k
}

# One trajectory: the population just before every passage from t0 to t0 + 120.
c2_simulate <- function(setting, landscape_id) {
  path <- c2_truth_path(setting$setting, landscape_id)
  if (file.exists(path)) return(readRDS(path))
  c2_bench_source()
  c2_load_abm()
  land <- c2_landscapes[c2_landscapes$landscape_id == landscape_id, ]
  founder <- rep(2L, 22L)
  set.seed(land$landscape_seed)
  centroids <- stt_generate_centroids(founder, 10L, land$wavelength)
  steps_per_passage <- as.integer(round(c2_abm$passage_days / c2_abm$delta_t))
  t0 <- setting$burn_days
  record_days <- t0 + c2_ahead_days
  record_steps <- as.integer(round(record_days / c2_abm$delta_t))
  started <- proc.time()[["elapsed"]]
  raw <- c2_abm_cpp(
    initial_karyotypes = matrix(founder, nrow = 1L), initial_counts = setting$bottleneck,
    centroids = centroids, wavelength = land$wavelength, beta = setting$beta_sim, delta_t = c2_abm$delta_t,
    n_steps = max(record_steps), record_steps = record_steps, max_population = c2_abm$max_pop,
    culling_survival = 1, min_cn = c2_abm$min_cn, max_cn = as.integer(setting$max_cn),
    seed = land$landscape_seed * 100 + setting$seed_offset, constant_fitness = 0, landscape_type = 0L,
    omega = setting$omega, beta_factor = setting$beta_factor, wgd_gain = setting$gain,
    passage_steps = steps_per_passage, bottleneck = setting$bottleneck, ploidy_aware = as.integer(setting$ploidy_aware))
  snapshots <- lapply(record_steps, function(step) {
    use <- raw$step == step
    stats::setNames(raw$count[use], raw$karyotype[use])
  })
  names(snapshots) <- as.character(c2_ahead_days)
  passages <- as.data.frame(raw$passages)
  passages$day <- passages$step * c2_abm$delta_t
  truth <- list(landscape_id = landscape_id, setting = setting$setting, wavelength = land$wavelength,
                centroids = centroids, abm_delta_t = c2_abm$delta_t, t0 = t0, beta_sim = setting$beta_sim,
                bottleneck = setting$bottleneck, omega = setting$omega, gain = setting$gain,
                beta_factor = setting$beta_factor, ploidy_aware = setting$ploidy_aware, max_cn = setting$max_cn,
                snapshots = snapshots, passages = passages[passages$day >= t0 - 1e-9, , drop = FALSE],
                seconds = proc.time()[["elapsed"]] - started)
  dev_save_rds(truth, path)
  truth
}

c2_dataset_id <- function(landscape_id, interval) sprintf("%s_i%03d", landscape_id, as.integer(interval))

# The four datasets of a trajectory, with their scored nodes.
c2_make_datasets <- function(setting, landscape_id) {
  c2_bench_source()
  truth <- c2_simulate(setting, landscape_id)
  out <- character()
  for (interval in c2_intervals) {
    id <- c2_dataset_id(landscape_id, interval)
    path <- file.path(c2_root, "grid", setting$setting, "datasets", paste0(id, ".rds"))
    if (!file.exists(path)) {
      p0 <- truth$snapshots[["0"]]
      p1 <- truth$snapshots[[as.character(interval)]]
      if (!length(p0) || !length(p1)) stop("Empty population in ", setting$setting, " ", id, call. = FALSE)
      seed <- as.integer((c2_landscapes$landscape_seed[c2_landscapes$landscape_id == landscape_id] %% 1000L) * 1000L +
                           interval + setting$seed_offset * 7L)
      counts <- stt_sample_counts(p0, p1, c2_depth[[1L]], c2_depth[[2L]], interval, seed)
      storage.mode(counts) <- "integer"
      d <- list(dataset_id = id, landscape_id = landscape_id, setting = setting$setting, wavelength = truth$wavelength,
                t0 = truth$t0, interval = interval, depth_t0 = c2_depth[[1L]], depth_t1 = c2_depth[[2L]], seed = seed,
                beta = c2_fit_beta, max_cn = as.integer(setting$max_cn), counts = counts,
                truth = list(centroids = truth$centroids, wavelength = truth$wavelength, abm_delta_t = truth$abm_delta_t,
                             beta_sim = truth$beta_sim, beta_factor = truth$beta_factor, omega = truth$omega,
                             gain = truth$gain, ploidy_aware = truth$ploidy_aware,
                             frequency_t0 = p0 / sum(p0), frequency_t1 = p1 / sum(p1),
                             truth_file = file.path("grid", setting$setting, "truth", paste0(landscape_id, ".rds")),
                             ahead_days = c2_ahead_days[c2_ahead_days > interval]))
      dev_save_rds(d, path)
      nodes <- stt_scoring_nodes(counts, 1L, as.integer(setting$max_cn), 2000L, 2000L, seed)
      dev_save_rds(nodes, file.path(c2_root, "grid", setting$setting, "nodes", paste0(id, ".rds")))
    }
    out <- c(out, id)
  }
  out
}

c2_dataset <- function(key, id) readRDS(file.path(c2_dir(key), "datasets", paste0(id, ".rds")))
c2_nodes <- function(key, id) readRDS(file.path(c2_dir(key), "nodes", paste0(id, ".rds")))

# Expected per-day growth rate of each karyotype in the ABM (see the header).
c2_true_rate <- function(labels, d) {
  c2_load_abm()
  k <- c2_parse(labels)
  t <- d$truth
  f <- c2_landscape_fitness(k, t$centroids, t$wavelength, t$gain, as.integer(t$ploidy_aware))
  copies <- rowSums(k)
  ploidy4 <- rowMeans(k) >= 3
  beta_k <- ifelse(ploidy4 & t$omega > 0, pmin(t$beta_sim * t$beta_factor, 0.999), t$beta_sim)
  p_any <- 1 - (1 - beta_k)^copies
  log1p(pmax(f, 0) * (1 - 2 * p_any) * t$abm_delta_t) / t$abm_delta_t
}

# Share of a column's cells in karyotypes with mean copy number >= 3.
c2_ploidy4_share <- function(counts, column) {
  k <- c2_parse(rownames(counts))
  sum(counts[rowMeans(k) >= 3, column]) / sum(counts[, column])
}

# Statistics of the counts alone (acceptance; the fit gives coverage and dispersion).
c2_count_stats <- function(d) {
  x <- d$counts
  data.frame(dataset_id = d$dataset_id, interval = d$interval, wavelength = d$wavelength,
             karyotypes = nrow(x), singleton_share_first = sum(x[, 1] == 1) / sum(x[, 1]),
             singleton_share_last = sum(x[, 2] == 1) / sum(x[, 2]),
             shared = sum(x[, 1] > 0 & x[, 2] > 0),
             ploidy4_first = c2_ploidy4_share(x, 1L), ploidy4_last = c2_ploidy4_share(x, 2L),
             mean_cn_last = sum(rowMeans(c2_parse(rownames(x))) * x[, 2]) / sum(x[, 2]),
             stringsAsFactors = FALSE)
}

# Pearson dispersion at the second time point (rd_dispersion of dev/realdata/R/realdata.R).
c2_dispersion <- function(s) {
  seen <- s$count_t1 > 0 | s$expected_t1 >= 1
  if (!any(seen)) return(NA_real_)
  sum(((s$count_t1 - s$expected_t1)^2 / pmax(s$expected_t1, 1e-8))[seen]) / sum(seen)
}
