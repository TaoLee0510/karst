# Non-GRF landscape families of D1 (dev/PLAN.md): rough Mount Fuji
# (landscape_type 1) and NK on the copy-number lattice with K = 2
# (landscape_type 2) of dev/abm/abm_dev.cpp, simulated with the design of the
# development benchmark (benchmark/sparse_two_timepoint: its ABM settings,
# truth trajectories, core scenarios, sampling and scored nodes).
#
# Sourced after common.R. Other scripts load a dataset and its truth with
#   config  <- dev_family_config("rmf")             # or "nk"
#   dataset <- dev_family_dataset(config, dataset_id)
#   nodes   <- dev_family_nodes(config, dataset)
#   truth   <- dev_family_node_truth(config, dataset, nodes)   # cached; or
#              dev_family_true_growth_rate(labels, dataset$truth) for any labels
# exactly as the development benchmark uses stt_dataset_nodes() and
# stt_true_growth_rate() (truth = per-day growth rate, scored relative to the
# reference by stt_score_predictions()).

dev_family_seed_base <- c(rmf = 61000L, nk = 62000L)

# Parameter sets: set 1 is the first choice; set 2 the one allowed adjustment
# (dev/PLAN.md, D1). Set 1 of RMF gave a median of 1 major (GRF 10% quantile
# 2): its populations sat on the summit, four steps from the founder. Set 2
# weakens selection (slope 0.015, noise 0.02 per karyotype) and puts the summit
# eight steps away, so that adaptation is still going on at day 1440. NK passed
# with set 1 and keeps it.
dev_family_parameters <- list(
  rmf = list(`1` = list(top = 1.0, slope = 0.05, sd = 0.06, summit_steps = 4L),
             `2` = list(top = 1.0, slope = 0.015, sd = 0.02, summit_steps = 8L)),
  nk = list(`1` = list(mean = 0.55, sd = 0.25, K = 2L))
)
# The parameter set in use per family.
dev_family_param_set <- c(rmf = "2", nk = "1")

dev_family_config <- function(family, param_set = dev_family_param_set[[family]]) {
  dev_source_validation()
  config <- stt_config(file.path(dev_results, "phase0", "d1", family))
  index <- seq_len(6L)
  config$landscapes <- data.frame(
    landscape_id = sprintf("%s_L%02d", family, index), split = family, wavelength_class = family,
    wavelength = NA_real_, landscape_seed = dev_family_seed_base[[family]] + index, stringsAsFactors = FALSE)
  config$family <- family
  config$param_set <- as.character(param_set)
  config$family_par <- dev_family_parameters[[family]][[config$param_set]]
  if (is.null(config$family_par)) stop("No parameter set ", param_set, " for ", family, call. = FALSE)
  config
}

# The landscape of one family member: type, parameters and reference karyotype
# for dev_abm_cpp() / dev_landscape_fitness(). Like stt_generate_centroids(),
# a draw is repeated until the diploid founder divides faster than 0.4 per day.
dev_family_landscape <- function(config, landscape_seed) {
  dev_load_abm()
  founder <- rep(config$founder_cn, config$n_chr)
  p <- config$family_par
  for (draw in 0:999) {
    seed <- landscape_seed * 1000 + draw
    if (config$family == "rmf") {
      set.seed(seed)
      ref <- founder
      chosen <- sample.int(config$n_chr, p$summit_steps)
      ref[chosen] <- ref[chosen] + sample(c(-1L, 1L), p$summit_steps, replace = TRUE)
      spec <- list(type = 1L, par = c(p$top, p$slope, p$sd, seed), ref = as.integer(ref))
    } else {
      spec <- list(type = 2L, par = c(p$mean, p$sd, p$K, seed), ref = integer())
    }
    f0 <- dev_landscape_fitness(matrix(as.integer(founder), nrow = 1L), spec$type, spec$par, spec$ref)
    if (f0 > 0.4) return(c(spec, list(founder_fitness = f0, draw = draw)))
  }
  stop("No landscape draw with founder fitness above 0.4.", call. = FALSE)
}

dev_family_fitness <- function(labels, spec) {
  dev_load_abm()
  k <- if (is.character(labels)) stt_parse_karyotypes(labels) else as.matrix(labels)
  storage.mode(k) <- "integer"
  dev_landscape_fitness(k, spec$type, spec$par, spec$ref)
}

dev_family_true_growth_rate <- function(labels, truth) {
  stt_effective_growth_rate(dev_family_fitness(labels, truth$landscape), truth$abm_delta_t)
}

# The true growth rates of a dataset's scored nodes, cached on disk (so that
# fits on HPC need not compile the dev ABM for dev_landscape_fitness()).
dev_family_node_truth <- function(config, dataset, nodes = dev_family_nodes(config, dataset)) {
  path <- file.path(config$output_dir, "truth_rates", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) {
    cached <- readRDS(path)
    if (identical(names(cached), nodes$karyotype)) return(unname(cached))
  }
  rate <- dev_family_true_growth_rate(nodes$karyotype, dataset$truth)
  dev_save_rds(stats::setNames(rate, nodes$karyotype), path)
  rate
}

# stt_run_abm() replaced for the duration of `expr`: the restart of the
# bottleneck scenario (stt_restart_snapshot()) then runs on the family landscape.
dev_family_with_abm <- function(spec, expr) {
  original <- get("stt_run_abm", envir = globalenv())
  assign("stt_run_abm", function(initial, centroids, wavelength, config, record_days, seed, beta = config$beta) {
    dev_run_abm(initial, config$abm, record_days = record_days, seed = seed, beta = beta,
                culling_survival = config$abm$culling_survival, landscape_type = spec$type,
                landscape_par = spec$par, landscape_ref = spec$ref)$populations
  }, envir = globalenv())
  on.exit(assign("stt_run_abm", original, envir = globalenv()), add = TRUE)
  force(expr)
}

# One truth trajectory (patient 1), as stt_generate_truth() makes it.
dev_family_truth <- function(config, landscape) {
  path <- stt_truth_path(config, landscape$landscape_id, 1L)
  if (file.exists(path)) return(path)
  spec <- dev_family_landscape(config, landscape$landscape_seed)
  founder <- rep(config$founder_cn, config$n_chr)
  abm <- config$abm
  initial <- stats::setNames(abm$pop_size, paste(founder, collapse = "."))
  run <- dev_run_abm(initial, abm, record_days = seq(0, abm$end_day, by = abm$record_days),
                     seed = landscape$landscape_seed * 100 + 1, beta = config$beta,
                     culling_survival = abm$culling_survival, landscape_type = spec$type,
                     landscape_par = spec$par, landscape_ref = spec$ref)
  truth <- list(landscape_id = landscape$landscape_id, patient = 1L, split = landscape$split,
                wavelength_class = landscape$wavelength_class, wavelength = NA_real_,
                centroids = matrix(0, 1L, config$n_chr), beta = config$beta, abm_delta_t = abm$delta_t,
                snapshots = run$populations, landscape = spec, family = config$family,
                param_set = config$param_set)
  dev_save_rds(truth, path)
  path
}

# The core scenarios of the development benchmark, with their indices in the
# full scenario table (the seed of a dataset depends on it, as in
# stt_dataset_plan()).
dev_family_plan <- function(config) {
  core <- which(config$scenarios$role == "core")
  plan <- list()
  for (i in seq_len(nrow(config$landscapes))) for (s in core) {
    plan[[length(plan) + 1L]] <- list(
      landscape = config$landscapes[i, ], scenario = config$scenarios[s, ],
      seed = as.integer(config$landscapes$landscape_seed[i] %% 100000L * 1000L + 1L * 20L + s))
  }
  plan
}

dev_family_make_dataset <- function(config, item) {
  truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, 1L))
  path <- dev_family_with_abm(truth$landscape,
                              stt_make_dataset(config, truth, item$scenario, config$t0, item$seed))
  dataset <- readRDS(path)
  if (is.null(dataset$truth$landscape)) {
    dataset$truth$landscape <- truth$landscape
    dataset$family <- config$family
    dev_save_rds(dataset, path)
  }
  path
}

dev_family_dataset <- function(config, dataset_id) readRDS(stt_dataset_path(config, dataset_id))

dev_family_nodes <- function(config, dataset) {
  path <- file.path(config$output_dir, "nodes", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  s <- config$scoring
  nodes <- stt_scoring_nodes(dataset$counts, s$min_cn, s$max_cn, s$d1_cap, s$d2_cap, dataset$seed)
  dev_save_rds(nodes, path)
  nodes
}

# Population statistics of the D1 rule, from the counts alone: number of
# majors and share of cells inside the model domain as fit_sparse_joint()
# builds it (major frequency 0.01, 3 cells, 5000 states, copy numbers 1-8),
# and the number of distinct karyotypes.
dev_family_count_stats <- function(counts) {
  storage.mode(counts) <- "integer"
  domain <- karst:::sparse_joint_domain(counts, 0.01, 3, 5000, 1L, 8L)
  inside <- rownames(counts) %in% domain$labels
  data.frame(n_major = domain$n_major, domain_share = sum(counts[inside, ]) / sum(counts),
             n_karyotypes = nrow(counts))
}

# The base fit of one dataset, as the development benchmark's
# karst_sparse_joint_balance method makes it, scored against the family truth.
dev_family_fit <- function(config, dataset_id) {
  dataset <- dev_family_dataset(config, dataset_id)
  nodes <- dev_family_nodes(config, dataset)
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(karst::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = 1,
                                          min_cn = 1L, max_cn = 8L, initial_state = "balance"),
                  error = function(e) e)
  seconds <- proc.time()[["elapsed"]] - started
  key <- list(dataset_id = dataset_id, family = config$family, method = "karst_sparse_joint_balance",
              seconds = seconds, library = dirname(find.package("karst")))
  if (inherits(fit, "error")) return(c(key, list(status = "failed", detail = conditionMessage(fit))))
  p <- stats::predict(fit, nodes$karyotype)
  prediction <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd,
                           model_sd = p$relative_sd_model, shell = p$shell, stringsAsFactors = FALSE)
  truth_rate <- dev_family_node_truth(config, dataset, nodes)
  score <- stt_score_predictions(nodes, truth_rate, prediction, config$scoring$reference_min_count)
  d <- fit$diagnostics
  c(key, list(status = "ok", n_states = d$n_states, n_major = d$n_major, history = d$history,
              convergence = d$convergence, prediction = prediction, metrics = score$metrics,
              nodes = score$nodes))
}
