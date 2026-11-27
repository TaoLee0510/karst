# R10, WGD part (dev/phase3/r10_spec.md): data from the WGD ABM
# (dev/abm/abm_wgd.cpp) with the development benchmark's design, and the
# truth used for scoring. Sourced after common.R.

r10_load_abm <- function() {
  if (exists("wgd_abm_cpp", mode = "function")) return(invisible())
  cache <- file.path(dev_results, "cpp_cache_wgd", Sys.info()[["nodename"]])
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(dev_repo, "dev", "abm", "abm_wgd.cpp"), cacheDir = cache, env = globalenv())
  invisible()
}

r10_source_bench <- function() {
  if (exists("stt_make_dataset", mode = "function")) return(invisible())
  for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R")) {
    source(file.path(dev_repo, "benchmark", "sparse_two_timepoint", part))
  }
}

r10_omega <- 1e-4
r10_gain_sets <- c(`1` = 0.02, `2` = 0.05)

r10_config <- function(gain_set = "1") {
  r10_source_bench()
  config <- stt_config(file.path(dev_results, "phase3", "r10", paste0("wgd_g", gain_set)))
  config$abm$max_cn <- 10L
  config$scoring$max_cn <- 10L
  index <- seq_len(6L)
  config$landscapes <- data.frame(landscape_id = sprintf("wgd_L%02d", index), split = "wgd", wavelength_class = "w0.8",
                                  wavelength = 0.8, landscape_seed = 63000L + index, stringsAsFactors = FALSE)
  config$omega <- r10_omega
  config$gain <- r10_gain_sets[[gain_set]]
  config
}

r10_run_abm <- function(initial, config, record_days, seed, beta, culling_survival, centroids, wavelength) {
  r10_load_abm()
  abm <- config$abm
  record_steps <- as.integer(round(sort(unique(record_days)) / abm$delta_t))
  raw <- wgd_abm_cpp(
    initial_karyotypes = matrix(as.integer(do.call(rbind, strsplit(names(initial), ".", fixed = TRUE))), nrow = length(initial)),
    initial_counts = as.numeric(initial), centroids = centroids, wavelength = wavelength, beta = beta,
    delta_t = abm$delta_t, n_steps = max(record_steps), record_steps = record_steps,
    max_population = abm$max_pop, culling_survival = culling_survival, min_cn = abm$min_cn, max_cn = abm$max_cn,
    seed = seed, constant_fitness = 0, landscape_type = 0L, omega = config$omega, beta_factor = 2.5,
    wgd_gain = config$gain)
  out <- lapply(stats::setNames(record_steps, as.character(record_steps * abm$delta_t)), function(step) {
    use <- raw$step == step
    stats::setNames(raw$count[use], raw$karyotype[use])
  })
  out
}

r10_truth_path <- function(config, landscape_id) file.path(config$output_dir, "truth", paste0(landscape_id, "_P01.rds"))

r10_truth <- function(config, i) {
  landscape <- config$landscapes[i, ]
  path <- r10_truth_path(config, landscape$landscape_id)
  if (file.exists(path)) return(readRDS(path))
  founder <- rep(config$founder_cn, config$n_chr)
  set.seed(landscape$landscape_seed)
  centroids <- stt_generate_centroids(founder, config$n_waves, landscape$wavelength)
  initial <- stats::setNames(config$abm$pop_size, paste(founder, collapse = "."))
  snapshots <- r10_run_abm(initial, config, seq(0, config$abm$end_day, by = config$abm$record_days),
                           landscape$landscape_seed * 100 + 1, config$beta, config$abm$culling_survival,
                           centroids, landscape$wavelength)
  truth <- list(landscape_id = landscape$landscape_id, patient = 1L, split = "wgd", wavelength_class = "w0.8",
                wavelength = landscape$wavelength, centroids = centroids, beta = config$beta,
                abm_delta_t = config$abm$delta_t, snapshots = snapshots, omega = config$omega, gain = config$gain)
  dev_save_rds(truth, path)
  truth
}

# Core scenarios of the development benchmark, made with stt_make_dataset()
# while the WGD ABM stands in for stt_run_abm() (restarts).
r10_make_datasets <- function(config, i) {
  truth <- r10_truth(config, i)
  core <- config$scenarios[config$scenarios$role == "core", ]
  original <- get("stt_run_abm", envir = globalenv())
  on.exit(assign("stt_run_abm", original, envir = globalenv()))
  assign("stt_run_abm", function(initial, centroids, wavelength, config_inner, record_days, seed, beta = config_inner$beta) {
    r10_run_abm(initial, config, record_days, seed, beta, config_inner$abm$culling_survival, centroids, wavelength)
  }, envir = globalenv())
  ids <- character()
  for (s in seq_len(nrow(core))) {
    seed <- as.integer(config$landscapes$landscape_seed[[i]] * 1000L + s)
    path <- stt_make_dataset(config, truth, core[s, ], config$t0, seed)
    d <- readRDS(path)
    if (is.null(d$omega)) {
      d$omega <- config$omega
      d$gain <- config$gain
      saveRDS(d, path)
    }
    ids <- c(ids, d$dataset_id)
  }
  ids
}

# True net growth rate per day at the true missegregation rate (2.5 x beta at ploidy 4).
r10_true_rate <- function(labels, dataset) {
  r10_load_abm()
  k <- do.call(rbind, lapply(strsplit(labels, ".", fixed = TRUE), as.integer))
  f <- wgd_landscape_fitness(k, dataset$truth$centroids, dataset$truth$wavelength, dataset$gain)
  self <- 2 * r10_self(k, dataset$beta) - 1
  log1p(pmax(f, 0) * self * dataset$truth$abm_delta_t) / dataset$truth$abm_delta_t
}

# Probability that a daughter keeps the parent's karyotype, 2.5 x beta at ploidy 4.
r10_self <- function(k, beta) {
  stay <- function(b) vapply(0:max(k), function(cn) karst:::karst_pij_cpp(cn, cn, b), numeric(1))
  s2 <- exp(rowSums(matrix(log(stay(beta))[k + 1L], nrow = nrow(k))))
  s4 <- exp(rowSums(matrix(log(stay(2.5 * beta))[k + 1L], nrow = nrow(k))))
  ifelse(rowMeans(k) >= 3, s4, s2)
}

r10_nodes <- function(config, dataset) {
  path <- file.path(config$output_dir, "nodes", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  nodes <- stt_scoring_nodes(dataset$counts, 1L, 10L, 2000L, 2000L, dataset$seed)
  dev_save_rds(nodes, path)
  nodes
}

r10_ploidy4_share <- function(counts, column = 1L) {
  k <- do.call(rbind, lapply(strsplit(rownames(counts), ".", fixed = TRUE), as.numeric))
  sum(counts[rowMeans(k) >= 3, column]) / sum(counts[, column])
}
