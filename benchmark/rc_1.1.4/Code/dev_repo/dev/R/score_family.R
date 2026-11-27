# R6 dosage prior (dev/phase3/r6_dosage_spec.md): the score landscape family,
# simulated with the development benchmark's design as dev/R/families.R does
# for RMF/NK, but on landscape_type 3 of dev/abm/abm_score.cpp:
#   f(k) = base + S * sum_c dd_c (k_c - 2) - q * sum_c (k_c - 2)^2 + sd * z(k),
# dd = Lynch 2022 Driver Density chromosome scores (p + q arm).
# Sourced after common.R and families.R (reuses dev_family_plan,
# dev_family_count_stats; families.R is not modified).

score_table <- function() {
  x <- utils::read.csv(file.path(dev_repo, "dev", "phase3", "lynch2022_table1.csv"), comment.char = "#",
                       stringsAsFactors = FALSE)
  x$chr <- sub("[pq]$", "", x$arm)
  x <- x[x$chr %in% as.character(1:22), ]
  agg <- stats::aggregate(cbind(gene_abundance, driver_density, hybrid) ~ chr, x, sum)
  agg[order(as.integer(agg$chr)), ]
}
score_dd <- function() score_table()$driver_density
# The prior of the KARST patch: Gene Abundance, centred and scaled to unit sd,
# times 0.01 (kappa = 1: one sd of score = 1% per copy).
score_ga_prior <- function() {
  ga <- score_table()$gene_abundance
  0.01 * (ga - mean(ga)) / stats::sd(ga)
}

score_seed_base <- 63000L
score_parameters <- list(`1` = list(base = 1.0, S = 0.25, q = 0.02, sd = 0.02))
score_param_set <- "1"

score_config <- function(param_set = score_param_set) {
  dev_source_validation()
  config <- stt_config(file.path(dev_results, "phase3", "r6d", "score"))
  index <- seq_len(6L)
  config$landscapes <- data.frame(
    landscape_id = sprintf("score_L%02d", index), split = "score", wavelength_class = "score",
    wavelength = NA_real_, landscape_seed = score_seed_base + index, stringsAsFactors = FALSE)
  config$family <- "score"
  config$param_set <- as.character(param_set)
  config$family_par <- score_parameters[[config$param_set]]
  config
}

score_load_abm <- function() {
  if (exists("score_abm_cpp", mode = "function")) return(invisible())
  cache <- file.path(dev_results, "cpp_cache_score", Sys.info()[["nodename"]])
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(dev_repo, "dev", "abm", "abm_score.cpp"), cacheDir = cache, env = globalenv())
  invisible()
}

score_run_abm <- function(initial, abm, record_days, seed, beta, culling_survival, spec) {
  score_load_abm()
  record_steps <- as.integer(round(sort(unique(record_days)) / abm$delta_t))
  raw <- score_abm_cpp(
    initial_karyotypes = matrix(as.integer(do.call(rbind, strsplit(names(initial), ".", fixed = TRUE))),
                                nrow = length(initial)),
    initial_counts = as.numeric(initial), centroids = matrix(0, 1, length(spec$dd)), wavelength = 1, beta = beta,
    delta_t = abm$delta_t, n_steps = max(record_steps), record_steps = record_steps,
    max_population = abm$max_pop, culling_survival = culling_survival,
    min_cn = abm$min_cn, max_cn = abm$max_cn, seed = seed, constant_fitness = 0,
    landscape_type = 3L, landscape_par = as.numeric(spec$par), landscape_ref = integer())
  lapply(stats::setNames(record_steps, as.character(record_steps * abm$delta_t)), function(step) {
    use <- raw$step == step
    stats::setNames(raw$count[use], raw$karyotype[use])
  })
}

score_landscape <- function(config, landscape_seed) {
  score_load_abm()
  founder <- rep(config$founder_cn, config$n_chr)
  p <- config$family_par
  dd <- score_dd()
  for (draw in 0:999) {
    seed <- landscape_seed * 1000 + draw
    spec <- list(type = 3L, par = c(p$base, p$S, p$q, p$sd, seed, dd), ref = integer(), dd = dd)
    f0 <- score_landscape_fitness(matrix(as.integer(founder), nrow = 1L), 3L, spec$par, integer())
    if (f0 > 0.4) return(c(spec, list(founder_fitness = f0, draw = draw)))
  }
  stop("No landscape draw with founder fitness above 0.4.", call. = FALSE)
}

score_fitness <- function(labels, spec) {
  score_load_abm()
  k <- if (is.character(labels)) stt_parse_karyotypes(labels) else as.matrix(labels)
  storage.mode(k) <- "integer"
  score_landscape_fitness(k, 3L, spec$par, integer())
}
score_true_growth_rate <- function(labels, truth) {
  stt_effective_growth_rate(score_fitness(labels, truth$landscape), truth$abm_delta_t)
}

score_with_abm <- function(spec, expr) {
  original <- get("stt_run_abm", envir = globalenv())
  assign("stt_run_abm", function(initial, centroids, wavelength, config, record_days, seed, beta = config$beta) {
    score_run_abm(initial, config$abm, record_days = record_days, seed = seed, beta = beta,
                  culling_survival = config$abm$culling_survival, spec = spec)
  }, envir = globalenv())
  on.exit(assign("stt_run_abm", original, envir = globalenv()), add = TRUE)
  force(expr)
}

score_truth <- function(config, landscape) {
  path <- stt_truth_path(config, landscape$landscape_id, 1L)
  if (file.exists(path)) return(path)
  spec <- score_landscape(config, landscape$landscape_seed)
  founder <- rep(config$founder_cn, config$n_chr)
  abm <- config$abm
  initial <- stats::setNames(abm$pop_size, paste(founder, collapse = "."))
  snapshots <- score_run_abm(initial, abm, record_days = seq(0, abm$end_day, by = abm$record_days),
                             seed = landscape$landscape_seed * 100 + 1, beta = config$beta,
                             culling_survival = abm$culling_survival, spec = spec)
  truth <- list(landscape_id = landscape$landscape_id, patient = 1L, split = landscape$split,
                wavelength_class = landscape$wavelength_class, wavelength = NA_real_,
                centroids = matrix(0, 1L, config$n_chr), beta = config$beta, abm_delta_t = abm$delta_t,
                snapshots = snapshots, landscape = spec, family = "score", param_set = config$param_set)
  dev_save_rds(truth, path)
  path
}

score_make_dataset <- function(config, item) {
  truth <- readRDS(stt_truth_path(config, item$landscape$landscape_id, 1L))
  path <- score_with_abm(truth$landscape, stt_make_dataset(config, truth, item$scenario, config$t0, item$seed))
  dataset <- readRDS(path)
  if (is.null(dataset$truth$landscape)) {
    dataset$truth$landscape <- truth$landscape
    dataset$family <- "score"
    dev_save_rds(dataset, path)
  }
  path
}

score_dataset <- function(config, id) readRDS(stt_dataset_path(config, id))

score_node_truth <- function(config, dataset, nodes) {
  path <- file.path(config$output_dir, "truth_rates", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) {
    cached <- readRDS(path)
    if (identical(names(cached), nodes$karyotype)) return(unname(cached))
  }
  rate <- score_true_growth_rate(nodes$karyotype, dataset$truth)
  dev_save_rds(stats::setNames(rate, nodes$karyotype), path)
  rate
}

# Far nodes (distance 3-4) and their truth, cached.
score_far <- function(config, dataset) {
  path <- file.path(config$output_dir, "far_nodes", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  nodes <- stt_far_nodes(counts, 1L, 8L, 2000L, dataset$seed + 3L)
  nodes$truth <- NA_real_
  nodes$truth[nodes$scored] <- score_true_growth_rate(nodes$karyotype[nodes$scored], dataset$truth)
  dev_save_rds(nodes, path)
  nodes
}
