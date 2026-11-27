# Version B of the paper's tests (alfak aligned to KARST, R/align.R) for one
# landscape of one input (source, data set, design): alfak at one minobs, or
# one KARST method (one landscape for all its minobs), remade on the input's
# extended list E. A forecast on E can take hours (cells that missegregate
# stay on the list), so each landscape is a task of its own.
#
#   ALIGN=karst Rscript scripts/align_task.R --make-tasks [N_PER_SOURCE]
#   ALIGN=karst Rscript scripts/align_task.R TASKS.tsv ID
#
# --make-tasks writes BENCH_RESULTS/tasks/inputs.tsv, one row per landscape
# (columns source, dataset, design, method, minobs: "5", "10", "20" for
# alfak, "all" for a KARST method) of every input of the round (BENCH_RESULTS
# is <round>/align_karst; ALIGN_RESULTS names another folder, e.g. for a
# pilot); with N_PER_SOURCE, a pilot: the landscapes of that many inputs per
# source, spread over its wavelengths, designs and data sets.
# A task, for each fit folder of its landscape (KARST: every minobs):
#   run.Rds         the round's (the fit is the same); a failed fit gets
#                   nothing else
#   landscape.Rds   k = E, mean, fq (the majors of KARST's default fit of the
#                   input), nn (their one-step neighbours), dist (0, 1, 2 to
#                   the observed karyotypes); alfak's also `source`
#                   ("landscape" or "kriging"). One per alfak minobs; one for
#                   each KARST method, the same in its minobs folders.
#   xval.Rds        alfak: the round's. KARST: folds centred on the fit's own
#                   majors, held-out predictions against targets from each
#                   karyotype's own counts (as_alfak(xval_target = "own"),
#                   resampled once with seed 1), R2R on the karyotypes with
#                   cells; `target` "own", `folds` "majors"
#   abm_preds/      forecasts on E (ABM and Salehi; sparse,
#                   align_abm_forecast()), as the round's (R/forecast.R),
#                   and forecast.Rds
#   score.Rds       ABM and missegregation: abm_fit_summary() (accuracy on E
#                   by distance class; forecasts scored on E)
#   bootstrap_res.Rds, karst_fit.Rds: links to the round's files
# and task_<minobs>.Rds in the method's folder: sizes, times and checks (E
# holds the input's karyotypes as they are written, alfak's Kriging on its own
# list, KARST's landscape on its list against the round's). Salehi forecasts
# are scored afterwards by scripts/score_salehi.R (ALIGN=karst).

for (f in list.files("R", full.names = TRUE)) source(f)
if (align_mode() != "karst") stop("run with ALIGN=karst")
paths <- bench_paths()
if (nzchar(Sys.getenv("ALIGN_RESULTS"))) paths$results <- Sys.getenv("ALIGN_RESULTS")
args <- commandArgs(TRUE)

if (identical(args[1], "--make-tasks")) {
  dirs <- Sys.glob(file.path(paths$results_base, "fits", c("abm", "misseg", "salehi"), "*", "*"))
  p <- strsplit(dirs, "/", fixed = TRUE)
  tasks <- data.frame(source = vapply(p, function(x) x[length(x) - 2L], ""),
                      dataset = vapply(p, function(x) x[length(x) - 1L], ""),
                      design = vapply(p, function(x) x[length(x)], ""), stringsAsFactors = FALSE)
  tasks <- tasks[order(tasks$source, tasks$dataset, tasks$design), ]
  if (length(args) >= 2L) {
    n <- as.integer(args[2])
    set.seed(1)
    group <- ifelse(tasks$source == "salehi", paste(tasks$source, tasks$design),
                    paste(tasks$source, sub("_rep_.*$", "", tasks$dataset), tasks$design))
    pick <- unlist(lapply(split(seq_len(nrow(tasks)), tasks$source), function(i) {
      g <- split(i, group[i])
      per <- max(1L, ceiling(n / length(g)))
      chosen <- unlist(lapply(g, function(j) j[sample.int(length(j), min(per, length(j)))]))
      chosen[sample.int(length(chosen), min(n, length(chosen)))]
    }))
    tasks <- tasks[sort(pick), ]
  }
  # one row per landscape
  tasks <- do.call(rbind, lapply(seq_len(nrow(tasks)), function(i) {
    d <- file.path(paths$results_base, "fits", tasks$source[i], tasks$dataset[i], tasks$design[i])
    ms <- intersect(c("alfak", karst_xv_methods), basename(list.dirs(d, recursive = FALSE)))
    do.call(rbind, lapply(ms, function(m) {
      mo <- if (m == "alfak") sub("minobs_", "", basename(Sys.glob(file.path(d, m, "minobs_*")))) else "all"
      if (!length(mo)) return(NULL)
      data.frame(tasks[i, ], method = m, minobs = mo, stringsAsFactors = FALSE)
    }))
  }))
  dir.create(file.path(paths$results, "tasks"), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tasks, file.path(paths$results, "tasks", "inputs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(table(tasks$source, tasks$design))
  quit(save = "no")
}

tasks <- utils::read.delim(args[1], stringsAsFactors = FALSE, colClasses = c(minobs = "character"))
task <- tasks[as.integer(args[2]), ]
src <- task$source; dataset <- task$dataset; design <- task$design
base_in <- file.path(paths$results_base, "fits", src, dataset, design)
out_in <- file.path(paths$results, "fits", src, dataset, design)
dir.create(file.path(out_in, task$method), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_in, task$method, paste0("task_", task$minobs, ".Rds"))
if (file.exists(log_file)) {
  cat("done already:", out_in, "\n")
  quit(save = "no")
}
elapsed <- function() proc.time()[["elapsed"]]
started <- elapsed()
record <- list(source = src, dataset = dataset, design = design, method = task$method, minobs = task$minobs)
settings <- bench_settings[[src]]
run_ok <- function(d) {
  p <- file.path(d, "run.Rds")
  file.exists(p) && identical(readRDS(p)$status, "ok")
}
link <- function(from, to) if (file.exists(from) && !file.exists(to)) file.symlink(from, to)

## Input and its list ---------------------------------------------------------------
if (src == "salehi") {
  yi <- salehi_input(dataset, ends = design == "ends", paths)
} else {
  sim <- read_abm_sim(dataset, if (src == "abm") paths$abm_raw else paths$misseg_raw)
  yi <- abm_input(sim, abm_designs[abm_designs$design == design, ])
}
x <- as.matrix(yi$x)
t0 <- elapsed()
L <- align_list(x)
E <- L$E
observed_k <- rownames(x)[rowSums(x) > 0]
record$E <- length(E)
record$observed <- sum(L$dist == 0L)
record$list_seconds <- elapsed() - t0
record$observed_as_written <- all(observed_k %in% E)
if (!record$observed_as_written) stop("the input's karyotypes are not written as E writes them")
cat(src, dataset, design, task$method, task$minobs, ": E", length(E), "observed", record$observed, "\n")
f_true <- NULL
if (src != "salehi") {
  t0 <- elapsed()
  f_true <- grf_fitness(E, abm_sim_info(dataset)$w, sim$true_landscape)
  record$truth_seconds <- elapsed() - t0
}

## The fits of the round on this input ---------------------------------------------------
present <- basename(list.dirs(base_in, recursive = FALSE))
methods <- task$method
minobs_dirs <- function(m) {
  d <- Sys.glob(file.path(base_in, m, "minobs_*"))
  if (m == "alfak") d[basename(d) == paste0("minobs_", task$minobs)] else d
}
load_fit <- function(m) {
  ok <- Filter(function(d) run_ok(d) && file.exists(file.path(d, "karst_fit.Rds")), minobs_dirs(m))
  if (!length(ok)) return(NULL)
  fit <- readRDS(file.path(ok[1], "karst_fit.Rds"))
  if (inherits(fit, "karst_sparse_joint_rate_fit")) fit <- fit$fits[[fit$point]]
  fit
}
# The frequent karyotypes of every landscape on this input: the majors of
# KARST's default fit (else of another KARST method; else, without a KARST
# fit, the observed karyotypes).
ref_fit <- NULL
for (m in intersect(karst_xv_methods, present)) {
  ref_fit <- load_fit(m)
  if (!is.null(ref_fit)) break
}
majors <- if (is.null(ref_fit)) observed_k else ref_fit$summary$karyotype[ref_fit$summary$major]
record$majors_from <- if (is.null(ref_fit)) "observed" else "karst"
nn <- setdiff(karst:::alfak_neighbours(majors), majors)
lscape_frame <- function(mean, source = NULL) {
  d <- data.frame(k = E, mean = mean, fq = E %in% majors, nn = E %in% nn, dist = L$dist, stringsAsFactors = FALSE)
  if (!is.null(source)) d$source <- source
  d
}

## Forecasts ---------------------------------------------------------------------------
fc <- bench_settings[[if (src == "salehi") "forecast_salehi" else "forecast_abm"]]
sim_time <- NULL
if (src == "salehi") {
  lineages <- readRDS(paths$salehi_lineages)
  meta <- utils::read.csv(paths$salehi_meta)
  sim_time <- salehi_forecast_time(dataset, lineages, meta)
}
forecast_landscape <- function(mean) {
  if (src == "misseg" || (src == "salehi" && is.null(sim_time))) return(NULL)
  if (src == "abm") {
    sx <- sim$abm_output$x
    times <- as.numeric(colnames(sx))
    start <- max(times[times < settings$last_time])
    report <- utils::head(times[times > start], fc$horizon) - start
    x0 <- sx[, times == start]
    names(x0) <- rownames(sx)
    at <- c(0, report)
    max_pop <- fc$abm_max_pop
  } else {
    last <- x[, ncol(x)]
    x0 <- last
    at <- seq(0, sim_time, by = 1)
    max_pop <- settings$nb
  }
  x0 <- x0[names(x0) %in% E & x0 > 0]
  if (sum(x0) <= 0) stop("no cell of the forecast's start is on E")
  x0 <- x0 / sum(x0)
  lapply(seq_len(fc$reps), function(r) {
    align_abm_forecast(E, mean, x0, at, settings$pm %||% 5e-5, fc$abm_pop_size, fc$abm_delta_t, max_pop,
                       fc$abm_culling_survival, fc$abm_record_interval, r)
  })
}
`%||%` <- function(a, b) if (is.null(a)) b else a
write_forecast <- function(dir, preds, seconds) {
  if (is.null(preds)) return(invisible())
  dir.create(file.path(dir, "abm_preds"), showWarnings = FALSE)
  for (r in seq_along(preds)) {
    save_rds_atomic(preds[[r]], file.path(dir, "abm_preds", sprintf(if (src == "abm") "rep_%02d.Rds" else "rep_%d.Rds", r)))
  }
  save_rds_atomic(list(status = "ok", seconds = seconds, message = NA_character_, sim_time = sim_time),
                  file.path(dir, "forecast.Rds"))
}
score_abm <- function(dir, m, minobs) {
  if (src == "salehi") return(invisible())
  info <- cbind(abm_sim_info(dataset), design = design, method = m, minobs = minobs)
  score <- abm_fit_summary(dir, sim, info, with_forecast = src == "abm", f = f_true)
  save_rds_atomic(score, file.path(dir, "score.Rds"))
}

record$fits <- list()
for (m in methods) {
  dirs <- minobs_dirs(m)
  out_dirs <- file.path(out_in, m, basename(dirs))
  for (i in seq_along(dirs)) {
    dir.create(out_dirs[i], recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(dirs[i], "run.Rds"), file.path(out_dirs[i], "run.Rds"), overwrite = TRUE)
    link(file.path(dirs[i], "bootstrap_res.Rds"), file.path(out_dirs[i], "bootstrap_res.Rds"))
    link(file.path(dirs[i], "karst_fit.Rds"), file.path(out_dirs[i], "karst_fit.Rds"))
  }
  ok <- vapply(dirs, run_ok, logical(1))
  if (!any(ok)) next
  if (m == "alfak") {
    for (i in which(ok)) {
      t0 <- elapsed()
      A <- align_alfak_landscape(dirs[i], E)
      land_seconds <- elapsed() - t0
      saveRDS(lscape_frame(A$mean, A$source), file.path(out_dirs[i], "landscape.Rds"))
      file.copy(file.path(dirs[i], "xval.Rds"), file.path(out_dirs[i], "xval.Rds"), overwrite = TRUE)
      t0 <- elapsed()
      preds <- forecast_landscape(A$mean)
      fseconds <- elapsed() - t0
      write_forecast(out_dirs[i], preds, fseconds)
      score_abm(out_dirs[i], m, sub("minobs_", "", basename(dirs[i])))
      record$fits[[length(record$fits) + 1L]] <- data.frame(method = m, minobs = basename(dirs[i]),
        kriged = sum(A$source == "kriging"), kriging_agreement = A$agreement, landscape_seconds = land_seconds,
        xval_seconds = NA_real_, forecast_seconds = fseconds, check_original = NA_real_, stringsAsFactors = FALSE)
    }
  } else {
    fit <- load_fit(m)
    if (is.null(fit)) next
    dir10 <- dirs[ok][basename(dirs[ok]) == "minobs_10"]
    if (!length(dir10)) dir10 <- dirs[ok][1]
    t0 <- elapsed()
    K <- align_karst_landscape(fit, dir10, E, bench_settings$correct_efflux)
    land_seconds <- elapsed() - t0
    # Cross-validation with folds centred on the fit's own majors, against
    # targets from each karyotype's own counts.
    t0 <- elapsed()
    fit_majors <- fit$summary$karyotype[fit$summary$major]
    members <- unique(c(fit_majors, karst:::alfak_neighbours(fit_majors)))
    members <- members[members %in% E]
    dv <- if (isTRUE(bench_settings$correct_efflux)) {
      2 * karst:::sparse_joint_self_probability(karst:::parse_karyotypes(members), fit$data$beta) - 1
    } else rep(1, length(members))
    names(dv) <- members
    own <- karst:::sparse_joint_own_fitness(fit, members, resample = TRUE, seed = 1L)
    target <- stats::setNames((as.numeric(own) - fit$trend$reference$fitness + K$level) / dv, members)
    observation <- stats::setNames(K$mean[match(members, E)], members)
    xv <- karst:::sparse_joint_alfak_xval(fit, fit_majors, observation, K$level, dv, target)
    attr(xv, "folds") <- NULL
    xv$target <- "own"
    xv$folds <- "majors"
    xval_seconds <- elapsed() - t0
    t0 <- elapsed()
    preds <- forecast_landscape(K$mean)
    fseconds <- elapsed() - t0
    lsc <- lscape_frame(K$mean)
    for (i in which(ok)) {
      saveRDS(lsc, file.path(out_dirs[i], "landscape.Rds"))
      saveRDS(xv, file.path(out_dirs[i], "xval.Rds"))
      write_forecast(out_dirs[i], preds, fseconds)
      score_abm(out_dirs[i], m, sub("minobs_", "", basename(dirs[i])))
    }
    record$fits[[length(record$fits) + 1L]] <- data.frame(method = m, minobs = "all", kriged = NA_integer_,
      kriging_agreement = NA_real_, landscape_seconds = land_seconds, xval_seconds = xval_seconds,
      forecast_seconds = fseconds, check_original = K$check_original, stringsAsFactors = FALSE)
  }
}
record$fits <- do.call(rbind, record$fits)
record$seconds <- elapsed() - started
record$peak_mb <- peak_rss_mb()
saveRDS(record, log_file)
print(record$fits)
cat(sprintf("done: E %d, %.0f s, peak %.0f MB\n", length(E), record$seconds, record$peak_mb))
