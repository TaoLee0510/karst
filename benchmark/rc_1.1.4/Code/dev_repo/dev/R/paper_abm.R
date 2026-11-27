# Evaluation harness on the ALFA-K paper's ABM simulations (dev/PLAN.md: R4,
# R2(iii), R6, R7). Fits KARST on an input exactly as the paper benchmark's
# karst arm does (benchmark/paper/R/fit.R fit_karst(): the window's count
# matrix, dt = diff(times), beta = pm, then karst::as_alfak() per minobs with
# nboot 500, n0, nb, correct_efflux TRUE and first_time), and scores the
# landscape with the benchmark's own code (benchmark/paper/R/metrics_abm.R
# landscape_scores(): r, rho, R, rnn, Rnn, rd2, Rxv, nfq), plus the
# cross-validation's held-out predictions against the true fitness
# (r_heldout_true, as scripts/initial_state_vs_truth.R). No forecasts.
#
# Sourced after dev/R/common.R. Data are read from PAPER_DATA (read only).
#
# Task table columns (dev/phase1/make_paper_abm_tasks.R):
#   method   name of the arm (output folder)
#   lib      karst library (dev/results/lib/<lib> or DEV_LIBDIR)
#   source   "abm" (raw/ABM) or "misseg" (raw/ABM_varyMisseg)
#   dataset  simulation id, e.g. w_0p8_m_5e-05_rep_01
#   design   w2_ends, w4_ends or w8_ends (benchmark/paper/R/config.R)
#   minobs   comma-separated minobs, e.g. "5,10,20"
#   pm       missegregation rate given to the fit (abm 5e-5; misseg the true rate)
#   extra    optional R expression returning a list of further
#            fit_sparse_joint() arguments. It is evaluated with these objects
#            in scope: info (abm_sim_info row: w, pmis), yi (input), times,
#            pm, passages (survivors of each passage between the first and the
#            last time point: n0 = 2e5 per passage, one per growth phase),
#            dataset, design. Example: "list(drift = dev_drift_variance(passages, 1.777))".

dev_paper_source <- function() {
  if (exists("landscape_scores", mode = "function") && exists("abm_input", mode = "function")) return(invisible())
  for (f in c("config.R", "abm.R", "metrics_abm.R")) {
    source(file.path(dev_repo, "benchmark", "paper", "R", f))
  }
  invisible()
}

dev_paper_data <- function() Sys.getenv("PAPER_DATA", file.path(Sys.getenv("RC_HOME"), "data", "paper"))  # rc freeze

dev_paper_raw <- function(source) {
  file.path(dev_paper_data(), "raw", if (source == "abm") "ABM" else "ABM_varyMisseg")
}

# Settings of the source (benchmark/paper/R/config.R bench_settings).
dev_paper_settings <- function(source) {
  dev_paper_source()
  s <- if (source == "abm") bench_settings$abm else bench_settings$misseg
  list(n0 = s$n0, nb = s$nb, nboot = bench_settings$nboot, correct_efflux = bench_settings$correct_efflux)
}

# The input of a design and its passage record. Each sampling time of the
# paper's ABM is a passage (sampled, then diluted to n0); between the first
# and the last time point of a window of W sampling times there are W - 1
# growth phases, each started by a dilution to n0.
dev_paper_input <- function(source, dataset, design_name) {
  dev_paper_source()
  sim <- read_abm_sim(dataset, dev_paper_raw(source))
  design <- abm_designs[abm_designs$design == design_name, ]
  if (nrow(design) != 1L) stop("Unknown design ", design_name)
  yi <- abm_input(sim, design)
  x <- as.matrix(yi$x)
  storage.mode(x) <- "double"
  times <- as.numeric(colnames(yi$x)) * yi$dt
  all_times <- as.numeric(colnames(sim$abm_output$x)) * sim$abm_output$dt
  phases <- sum(all_times >= min(times) & all_times < max(times))
  settings <- dev_paper_settings(source)
  list(sim = sim, info = abm_sim_info(dataset), yi = yi, x = x, times = times,
       passages = rep(settings$n0, phases), settings = settings)
}

dev_paper_fit <- function(task) {
  input <- dev_paper_input(task$source, task$dataset, task$design)
  minobs <- as.numeric(strsplit(as.character(task$minobs), ",", fixed = TRUE)[[1L]])
  extra <- list()
  if (!is.null(task$extra) && !is.na(task$extra) && nzchar(task$extra)) {
    scope <- list2env(list(info = input$info, yi = input$yi, times = input$times, pm = task$pm,
                           passages = input$passages, dataset = task$dataset, design = task$design),
                      parent = globalenv())
    extra <- eval(parse(text = task$extra), envir = scope)
  }
  key <- data.frame(method = task$method, source = task$source, dataset = task$dataset, design = task$design,
                    w = input$info$w, pmis = input$info$pmis, pm = task$pm, n_passages = length(input$passages),
                    stringsAsFactors = FALSE)
  started <- proc.time()[["elapsed"]]
  # extra arguments replace the defaults (e.g. beta)
  fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint,
                                           utils::modifyList(list(input$x, dt = diff(input$times), beta = task$pm),
                                                             extra))),
                  error = function(e) e)
  fit_seconds <- proc.time()[["elapsed"]] - started
  if (inherits(fit, "error")) {
    return(cbind(key, minobs = minobs, status = "failed", message = conditionMessage(fit), fit_seconds = fit_seconds,
                 stringsAsFactors = FALSE))
  }
  truth_of <- function(k) grf_fitness(k, input$info$w, input$sim$true_landscape)
  # Domain nodes for coverage (R2): the fit's relative fitness, model sd and
  # shell against the true net rate at the true missegregation rate, both
  # relative to the count-weighted mean of the majors. Saved beside the scores.
  nodes_path <- sub("\\.rds$", ".nodes.rds", dev_paper_output(task))
  try({
    s <- fit$summary
    cn <- karst:::parse_karyotypes(s$karyotype)
    net <- truth_of(s$karyotype) * (2 * karst:::sparse_joint_self_probability(cn, input$info$pmis) - 1)
    total <- rowSums(s[grep("^count_t", names(s))])
    w <- ifelse(s$major, total, 0); w <- w / sum(w)
    observed <- total > 0
    obs_cn <- cn[observed, , drop = FALSE]
    distance <- vapply(seq_len(nrow(cn)), function(i) min(rowSums(abs(sweep(obs_cn, 2, cn[i, ])))), numeric(1))
    dev_save_rds(data.frame(karyotype = s$karyotype, distance = distance, shell = s$shell, model_sd = s$relative_sd_model,
                            relative = s$relative_fitness, truth_relative = net - sum(w * net),
                            error = s$relative_fitness - (net - sum(w * net)), count_t0 = s$count_t0,
                            count_t1 = s[[grep("^count_t", names(s), value = TRUE)[2]]], stringsAsFactors = FALSE),
                 nodes_path)
  }, silent = TRUE)
  time0 <- "first_time" %in% names(formals(karst::as_alfak))
  rows <- lapply(minobs, function(m) {
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch({
      args <- list(fit, minobs = m, nboot = input$settings$nboot, n0 = input$settings$n0, nb = input$settings$nb,
                   correct_efflux = input$settings$correct_efflux)
      if (time0) args$first_time <- input$times[1]
      do.call(karst::as_alfak, args)
    }, error = function(e) e)
    if (inherits(out, "error")) {
      return(cbind(key, minobs = m, status = "failed", message = conditionMessage(out), fit_seconds = fit_seconds,
                   stringsAsFactors = FALSE))
    }
    lscape <- out$landscape
    s <- landscape_scores(lscape, truth_of(lscape$k), out$xval)
    p <- unname(out$xval$predictions)
    k <- out$xval$xval_data$k
    ok <- is.finite(p)
    held <- if (sum(ok) >= 3) {
      tr <- truth_of(k[ok])
      if (stats::sd(tr) > 0 && stats::sd(p[ok]) > 0) stats::cor(tr, p[ok]) else NA_real_
    } else NA_real_
    cbind(key, minobs = m, status = "ok", message = NA_character_, s, r_heldout_true = held,
          history = fit$diagnostics$history, objective = fit$diagnostics$objective,
          convergence = fit$diagnostics$convergence, fit_seconds = fit_seconds,
          alfak_seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
  })
  dev_paper_bind(rows)
}

# rbind of data frames whose columns differ (failed rows lack the scores).
dev_paper_bind <- function(rows) {
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    for (c in setdiff(cols, names(r))) r[[c]] <- NA
    r[cols]
  }))
}

dev_paper_output <- function(task) {
  file.path(dev_results, "phase1", "paper_abm", task$method, task$source,
            sprintf("%s__%s.rds", task$dataset, task$design))
}

# All rows of a method, one per input x design x minobs.
dev_paper_collect <- function(method) {
  files <- list.files(file.path(dev_results, "phase1", "paper_abm", method), recursive = TRUE,
                      pattern = "\\.rds$", full.names = TRUE)
  files <- files[!grepl("\\.nodes\\.rds$", files)]  # node tables of R2, not score rows
  dev_paper_bind(lapply(files, readRDS))
}
