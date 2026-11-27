# The paper's Fig. 2d on the ABM simulations, where the true landscape is
# known (docs/fig2d_truth_6ae2a24.md; plan approved 2026-11-23).
#
# Fig. 2d ranks, across forecasts, the observed change from a parent sample
# to a daughter sample (Wasserstein distance) by the predicted change (parent
# to forecast). On the Salehi passages alfak ranks higher than KARST; with
# the passages' sampling noise held fixed, no reliable difference is left
# (scripts/fig2d_noise.R). Here the same test runs on simulations with a known
# truth, with three forecasts per row: KARST's, alfak's and the true
# landscape's (the same karyotype list, the same forecast program,
# R/forecast.R). The parent is the last time point fitted (the last sampling
# time before day 120); the daughters are the next three sampling times.
#
# Step 1 takes the round's fits (1000 cells per time point) and adds sampling
# noise to the observed change only: parent and daughter subsampled to
# n = 100, 200 and 400 cells (20 draws each) besides the 1000 cells stored.
# Step 2 takes fits of inputs subsampled to the cell counts of the Salehi
# two-time-point inputs (ABM_CELLS, R/abm.R; the "tasks" stage writes their
# task tables, scripts/run_task.R runs them into F2T_ROOT/step2/); their
# observed change is the subsampled one, as on Salehi, besides the 1000 cells.
#
#   Rscript scripts/fig2d_truth.R --stage=cells    F2T_ROOT/salehi_cells.tsv (from the round's tables)
#   Rscript scripts/fig2d_truth.R --stage=tasks    F2T_ROOT/tasks/{sims,tasks}_{pilot,rest}.tsv
#   Rscript scripts/fig2d_truth.R --stage=engine   step 0: the stored forecasts of a few round fits redone with
#                                                  3 sampling times; F2T_ROOT/engine_check.tsv
#   Rscript scripts/fig2d_truth.R --stage=rows --step=1|2 --sims=pilot|rest --index=I
#                                                  one simulation: F2T_ROOT/step<k>/rows/<sim>.Rds
#                                                  (truth forecasts kept in F2T_ROOT/step<k>/truth/)
#   Rscript scripts/fig2d_truth.R --stage=analyse --step=1|2   (scripts/fig2d_truth_analyse.R)
#
# BENCH_RESULTS: the round for step 1, F2T_ROOT/step2 for step 2 (the fits
# the rows are made from). ABM_CELLS: F2T_ROOT/salehi_cells.tsv in step 2.
#
# --source=misseg (docs/fig2d_misseg_6ae2a24.md; plan approved 2026-11-23):
# step 2 on the paper's ABM_varyMisseg simulations (wavelength 0.8,
# missegregation rates 5e-5 to 8e-4, 25 each), whose populations are more
# spread out, as the Salehi passages are. Both methods are fitted twice,
# with the rate fixed at 5e-5 as on Salehi ("beta_fixed", for every rate) and
# with the true rate ("beta_true", rates above 5e-5 only; at 5e-5 the two
# coincide); each fit is forecast at the rate it was fitted with
# (MISSEG_FORECAST=TRUE, scripts/run_task.R). The true landscape is
# forecast at the true rate ("truth") and at 5e-5 ("truth_fixed").
#   Rscript scripts/fig2d_truth.R --source=misseg --stage=cells|tasks|dispersion
#   Rscript scripts/fig2d_truth.R --source=misseg --stage=rows --arm=beta_fixed|beta_true --sims=pilot|rest --index=I
#                         F2T_ROOT/<arm>/rows/<sim>.Rds; fits in BENCH_RESULTS = F2T_ROOT/<arm>; truth
#                         forecasts shared by the arms in F2T_ROOT/truth/<sim>/<design>/minobs_<m>/pm_<rate>/
# dispersion: F2T_ROOT/dispersion.tsv, how spread out the parent and daughter
# samples are (300 cells each) in these simulations, the main set's
# wavelength 0.8 and the Salehi passages, and the share of the Salehi
# parent cells on alfak's karyotype list.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
stage <- arg("stage")
source_set <- arg("source", "abm")
stopifnot(source_set %in% c("abm", "misseg"))
raw_dir <- if (source_set == "abm") paths$abm_raw else paths$misseg_raw
misseg_rates <- c(5e-5, 1e-4, 2e-4, 4e-4, 8e-4)
beta_fixed <- bench_settings$abm$pm
rate_label <- function(pm) format(pm, scientific = TRUE)
root <- Sys.getenv("F2T_ROOT", "")
if (!nzchar(root)) stop("F2T_ROOT is not set")
designs <- c("w2_ends", "w4_ends", "w8_ends")
methods <- c("karst", "alfak")
n_draws <- 20L
n_sub <- c(100L, 200L, 400L)
horizon <- 3L

s2m <- function(k) do.call(rbind, lapply(strsplit(k, ".", fixed = TRUE), as.numeric))
# Wasserstein distance (Euclidean ground distance, p = 1, as the paper)
# between two distributions given as masses on karyotype labels.
# transport's default solver (networkflow) can stall: on ABM_varyMisseg, 111
# parent karyotypes against a KARST forecast on 13,692 did not finish in an
# hour, while revsimplex found the same optimum (0.30297) in 3 s. Large
# problems therefore run networkflow in a child process and switch to
# revsimplex after wdist_timeout seconds; every other call is unchanged.
wdist_timeout <- as.numeric(Sys.getenv("WDIST_TIMEOUT", "60"))
wdist_fallbacks <- 0L
wdist <- function(ka, ma, kb, mb) {
  a <- is.finite(ma) & ma > 0
  b <- is.finite(mb) & mb > 0
  if (!any(a) || !any(b)) return(NaN)
  pa <- transport::wpp(s2m(ka[a]), mass = ma[a] / sum(ma[a]))
  pb <- transport::wpp(s2m(kb[b]), mass = mb[b] / sum(mb[b]))
  solve <- function(method) tryCatch(transport::wasserstein(pa, pb, method = method), error = function(e) NaN)
  if (sum(a) * sum(b) < 1e5) return(solve("networkflow"))
  job <- parallel::mcparallel(solve("networkflow"), silent = TRUE)
  res <- parallel::mccollect(job, wait = FALSE, timeout = wdist_timeout)
  if (!is.null(res) && is.numeric(res[[1]])) return(res[[1]])
  if (is.null(res)) {
    tools::pskill(job$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(job))
  }
  wdist_fallbacks <<- wdist_fallbacks + 1L
  solve("revsimplex")
}
# counts named by karyotype, restricted to the list `k` (zeros dropped)
on_list <- function(v, k) {
  v <- v[names(v) %in% k]
  v[v > 0]
}
w_counts <- function(a, b) wdist(names(a), a, names(b), b)
# Sampling noise of a sample, as scripts/fig2d_noise.R: the mean distance
# between two random halves of its cells (20 splits).
halves <- function(p, reps = 20L) {
  if (sum(p) < 4) return(NA_real_)
  cells <- rep(names(p), p)
  mean(vapply(seq_len(reps), function(r) {
    set.seed(r)
    i <- sample.int(length(cells), length(cells) %/% 2L)
    w_counts(table(cells[i]), table(cells[-i]))
  }, numeric(1)))
}
# n cells drawn without replacement from the counts `v`.
subsample <- function(v, n) {
  idx <- rep.int(seq_along(v), v)
  out <- stats::setNames(tabulate(idx[sample.int(length(idx), min(n, length(idx)))], length(v)), names(v))
  out[out > 0]
}
column <- function(x, time) {
  times <- as.numeric(colnames(x))
  v <- stats::setNames(as.numeric(x[, times == time]), rownames(x))
  v[v > 0]
}
fit_ok <- function(dir) {
  run <- file.path(dir, "run.Rds")
  file.exists(run) && identical(readRDS(run)$status, "ok") && file.exists(file.path(dir, "landscape.Rds"))
}
# mean of a fit's forecast replicates (rows: start and the next sampling
# times); NULL unless complete
mean_forecast <- function(files) {
  if (!length(files)) return(NULL)
  preds <- lapply(files, readRDS)
  mat <- Reduce(`+`, lapply(preds, function(x) as.matrix(x[, -1]))) / length(preds)
  list(time = preds[[1]]$time, mat = mat)
}
sim_list <- function(which, arm = NULL) {
  utils::read.delim(file.path(root, "tasks", paste0("sims_", if (!is.null(arm)) paste0(arm, "_"), which, ".tsv")),
                    stringsAsFactors = FALSE)$sim_id
}

## cells: the Salehi inputs' cells per time point --------------------------
if (stage == "cells") {
  si <- readRDS(file.path(paths$results_base, "tables", "salehi_inputs.Rds"))
  e <- si[si$design == "ends", ]
  cells <- data.frame(fi = rep(e$fi, 2), time_point = rep(c(1L, 2L), each = nrow(e)), cells = c(e$cells_1, e$cells_2))
  cells <- cells[is.finite(cells$cells) & cells$cells > 0, ]
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(cells, file.path(root, "salehi_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(summary(cells$cells))
}

## tasks: the simulations and step 2's fit tasks ----------------------------
# pilot: replicates 01-05 of each wavelength (20 simulations); rest: the
# other 380. KARST: one task per input (every minobs); alfak: one task per
# minobs (one fit per task).
if (stage == "tasks" && source_set == "misseg") {
  # pilot: replicates 01-02 of each rate (10 simulations); rest: the other
  # 115. beta_true leaves out the rate 5e-5 (its fits are beta_fixed's).
  sims <- sort(sub("\\.Rds$", "", list.files(raw_dir, pattern = "^w_.*_rep_.*\\.Rds$")))
  info <- do.call(rbind, lapply(sims, abm_sim_info))
  stopifnot(length(sims) == 125L, all(info$pmis %in% misseg_rates))
  dir.create(file.path(root, "tasks"), recursive = TRUE, showWarnings = FALSE)
  s <- bench_settings$misseg
  for (arm in c("beta_fixed", "beta_true")) for (which in c("pilot", "rest")) {
    keep <- (as.integer(info$abmrep) <= 2L) == (which == "pilot") & (arm == "beta_fixed" | info$pmis != beta_fixed)
    ids <- sims[keep]
    utils::write.table(data.frame(sim_id = ids), file.path(root, "tasks", paste0("sims_", arm, "_", which, ".tsv")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    rows <- list()
    for (id in ids) for (d in designs) {
      pm <- if (arm == "beta_fixed") beta_fixed else abm_sim_info(id)$pmis
      rows[[length(rows) + 1L]] <- data.frame(source = "misseg", dataset = id, design = d, method = "karst",
                                              minobs = "5,10,20", pm = pm)
      for (m in bench_settings$minobs) {
        rows[[length(rows) + 1L]] <- data.frame(source = "misseg", dataset = id, design = d, method = "alfak",
                                                minobs = as.character(m), pm = pm)
      }
    }
    t <- do.call(rbind, rows)
    t <- cbind(task = seq_len(nrow(t)), t[c("source", "dataset", "design", "method", "minobs")], n0 = s$n0,
               nb = s$nb, pm = t$pm)
    utils::write.table(t, file.path(root, "tasks", paste0("tasks_", arm, "_", which, ".tsv")), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    cat(arm, which, ":", length(ids), "simulations,", nrow(t), "tasks\n")
  }
}
if (stage == "tasks" && source_set == "abm") {
  sims <- sort(sub("\\.Rds$", "", list.files(paths$abm_raw, pattern = "^w_.*_m_5e-05_rep_.*\\.Rds$")))
  stopifnot(length(sims) == 400L)
  info <- do.call(rbind, lapply(sims, abm_sim_info))
  pilot <- as.integer(info$abmrep) <= 5L
  dir.create(file.path(root, "tasks"), recursive = TRUE, showWarnings = FALSE)
  s <- bench_settings$abm
  for (which in c("pilot", "rest")) {
    ids <- sims[if (which == "pilot") pilot else !pilot]
    utils::write.table(data.frame(sim_id = ids), file.path(root, "tasks", paste0("sims_", which, ".tsv")),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    rows <- list()
    for (id in ids) for (d in designs) {
      rows[[length(rows) + 1L]] <- data.frame(source = "abm", dataset = id, design = d, method = "karst",
                                              minobs = "5,10,20")
      for (m in bench_settings$minobs) {
        rows[[length(rows) + 1L]] <- data.frame(source = "abm", dataset = id, design = d, method = "alfak",
                                                minobs = as.character(m))
      }
    }
    t <- do.call(rbind, rows)
    t <- cbind(task = seq_len(nrow(t)), t, n0 = s$n0, nb = s$nb, pm = s$pm)
    utils::write.table(t, file.path(root, "tasks", paste0("tasks_", which, ".tsv")), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    cat(which, ":", length(ids), "simulations,", nrow(t), "tasks\n")
  }
}

## engine: step 0 ----------------------------------------------------------
# The forecast program with 3 sampling times gives the first rows of the
# round's stored forecasts (10 sampling times), replicate by replicate.
if (stage == "engine") {
  ids <- sim_list("pilot")
  ids <- ids[!duplicated(vapply(ids, function(i) abm_sim_info(i)$w, numeric(1)))]  # one per wavelength
  st <- bench_settings$forecast_abm
  st$horizon <- horizon
  out <- list()
  for (id in ids) {
    sim <- readRDS(file.path(paths$abm_raw, paste0(id, ".Rds")))
    for (m in methods) {
      dir <- file.path(paths$results, "fits", "abm", id, "w4_ends", m, "minobs_10")
      if (!fit_ok(dir) || !forecast_complete(dir)) next
      lscape <- readRDS(file.path(dir, "landscape.Rds"))
      for (r in seq_len(st$reps)) {
        stored <- readRDS(file.path(dir, "abm_preds", sprintf("rep_%02d.Rds", r)))
        t0 <- proc.time()[["elapsed"]]
        new <- forecast_abm(lscape, sim$abm_output$x, bench_settings$abm$pm, r, settings = st)
        sec <- proc.time()[["elapsed"]] - t0
        a <- as.matrix(new)
        b <- as.matrix(stored[seq_len(nrow(new)), ])
        out[[length(out) + 1L]] <- data.frame(sim_id = id, method = m, rep = r, rows = nrow(new),
                                              same_shape = identical(dim(a), dim(b)),
                                              max_abs_diff = if (identical(dim(a), dim(b))) max(abs(a - b)) else NA,
                                              seconds = sec)
      }
    }
  }
  out <- do.call(rbind, out)
  print(out, row.names = FALSE)
  utils::write.table(out, file.path(root, "engine_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

## dispersion: how spread out the samples are -----------------------------
# Parent (last sampling time before day 120) and daughter (two sampling times
# later) of every simulation of ABM_varyMisseg and of the main set's
# wavelength 0.8, and the last two passages of every Salehi input with at
# least 300 cells in both: each subsampled to 300 cells (10 draws), the
# Simpson diversity, the share of the most frequent karyotype, the share of
# cells whose karyotype occurs once, and the share of daughter cells whose
# karyotype the parent lacks. Salehi: also the share of the last input
# passage's cells on alfak's karyotype list (round fits, design "ends").
if (stage == "dispersion") {
  n_cells <- 300L
  measure <- function(p, d, reps = 10L) {
    p <- p[p > 0]
    d <- d[d > 0]
    r <- vapply(seq_len(reps), function(s) {
      set.seed(s)
      a <- subsample(p, n_cells)
      b <- subsample(d, n_cells)
      pa <- a / sum(a)
      c(simpson = 1 - sum(pa^2), top = max(pa), once = sum(a == 1) / sum(a),
        new_in_daughter = sum(b[!names(b) %in% names(a)]) / sum(b))
    }, numeric(4))
    as.list(rowMeans(r))
  }
  sim_rows <- function(dir, set, pattern) {
    ids <- sub("\\.Rds$", "", list.files(dir, pattern = pattern))
    do.call(rbind, lapply(ids, function(id) {
      x <- as.matrix(readRDS(file.path(dir, paste0(id, ".Rds")))$abm_output$x)
      tt <- as.numeric(colnames(x))
      i0 <- max(which(tt < bench_settings$abm$last_time))
      data.frame(set = set, id = id, line = NA_character_, pmis = abm_sim_info(id)$pmis, dt = tt[i0 + 2L] - tt[i0],
                 measure(column(x, tt[i0]), column(x, tt[i0 + 2L])), on_list = NA_real_)
    }))
  }
  out <- rbind(sim_rows(paths$misseg_raw, "ABM_varyMisseg", "^w_.*[.]Rds$"),
               sim_rows(paths$abm_raw, "ABM main set", "^w_0p8_m_5e-05_rep_.*[.]Rds$"))
  sal <- list()
  for (fi in salehi_trajectories(paths)$traj_id) {
    y <- readRDS(file.path(paths$salehi_inputs, paste0(fi, ".Rds")))
    x <- y$x
    p <- stats::setNames(as.numeric(x[, ncol(x) - 1L]), rownames(x))
    d <- stats::setNames(as.numeric(x[, ncol(x)]), rownames(x))
    if (sum(p) < n_cells || sum(d) < n_cells) next
    # alfak's list (every minobs) against the last passage of the "ends" input
    shares <- vapply(bench_settings$minobs, function(m) {
      f <- file.path(paths$results_base, "fits", "salehi", fi, "ends", "alfak", paste0("minobs_", m), "landscape.Rds")
      if (!file.exists(f)) return(NA_real_)
      sum(d[names(d) %in% readRDS(f)$k]) / sum(d)
    }, numeric(1))
    tt <- as.numeric(colnames(x))
    sal[[length(sal) + 1L]] <- data.frame(set = "Salehi", id = fi, line = sub("_.*$", "", fi), pmis = NA_real_,
                                          dt = (tt[ncol(x)] - tt[ncol(x) - 1L]) * y$dt, measure(p, d),
                                          on_list = mean(shares, na.rm = TRUE))
  }
  out <- rbind(out, do.call(rbind, sal))
  # one row per distinct Salehi pair of passages (trajectories share them)
  out <- out[!(out$set == "Salehi" & duplicated(out[c("set", "simpson", "top", "once", "new_in_daughter")])), ]
  utils::write.table(out, file.path(root, "dispersion.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  out$group <- ifelse(is.na(out$pmis), out$set, paste(out$set, rate_label(out$pmis)))
  print(stats::aggregate(cbind(simpson, top, once, new_in_daughter, on_list) ~ group, out, stats::median,
                         na.action = stats::na.pass))
  # the Salehi lines differ widely in how spread out they are
  print(stats::aggregate(cbind(simpson, top, once, new_in_daughter, on_list) ~ line, out[out$set == "Salehi", ],
                         stats::median, na.action = stats::na.pass))
}

## rows: one simulation ----------------------------------------------------
if (stage == "rows") {
  arm <- NULL
  if (source_set == "misseg") {
    arm <- arg("arm")
    stopifnot(arm %in% c("beta_fixed", "beta_true"))
    step <- 2L
  } else {
    step <- as.integer(arg("step"))
    stopifnot(step %in% 1:2)
  }
  ids <- sim_list(arg("sims", "pilot"), arm)
  index <- as.integer(arg("index", Sys.getenv("SLURM_ARRAY_TASK_ID")))
  id <- ids[[index]]
  info <- abm_sim_info(id)
  out_dir <- if (is.null(arm)) file.path(root, paste0("step", step)) else file.path(root, arm)
  dir.create(file.path(out_dir, "rows"), recursive = TRUE, showWarnings = FALSE)
  sim0 <- readRDS(file.path(raw_dir, paste0(id, ".Rds")))
  # the rate the fits were made and forecast with; the true landscape is
  # forecast at the true rate and, for the missegregation simulations, at 5e-5
  pm_fit <- if (identical(arm, "beta_true")) info$pmis else beta_fixed
  truth_pm <- c(truth = info$pmis)
  if (source_set == "misseg") truth_pm <- c(truth_pm, truth_fixed = beta_fixed)
  x_full <- sim0$abm_output$x
  if (step == 2L) {
    cells_file <- Sys.getenv("ABM_CELLS", "")
    if (!nzchar(cells_file)) stop("step 2 needs ABM_CELLS (the cells table the fits were subsampled with)")
    x_fit <- abm_subsample_x(x_full, id, utils::read.delim(cells_file)$cells)
  } else {
    x_fit <- x_full
  }
  x_full <- as.matrix(x_full)
  x_fit <- as.matrix(x_fit)
  times <- as.numeric(colnames(x_full))
  t0 <- max(times[times < bench_settings$abm$last_time])
  td <- utils::head(times[times > t0], horizon)
  p_full <- column(x_full, t0)
  p_fit <- column(x_fit, t0)
  d_full <- lapply(td, function(t) column(x_full, t))
  d_fit <- lapply(td, function(t) column(x_fit, t))
  # step 1: subsamples of parent and daughter (draw d, size n), the same
  # for every karyotype list of the simulation
  draws <- NULL
  if (step == 1L) {
    base_seed <- (abm_sim_seed(id) - 9100000L) * 10000L
    draws <- lapply(seq_along(td), function(j) lapply(seq_along(n_sub), function(a) lapply(seq_len(n_draws), function(d) {
      set.seed(base_seed + j * 1000L + a * 100L + d)
      list(p = subsample(p_full, n_sub[a]), d = subsample(d_full[[j]], n_sub[a]))
    })))
  }
  st <- bench_settings$forecast_abm
  st$horizon <- horizon
  preds <- list()
  obs <- list()
  drw <- list()
  fits <- list()
  started <- proc.time()[["elapsed"]]
  truth_seconds <- 0
  for (design in designs) for (mo in bench_settings$minobs) {
    dirs <- stats::setNames(file.path(paths$results, "fits", source_set, id, design, methods, paste0("minobs_", mo)),
                            methods)
    ok <- vapply(dirs, fit_ok, logical(1))
    if (!any(ok)) next
    # the karyotype list: alfak's (KARST's is the same list in version A,
    # tables/ksets.Rds); KARST's where alfak's fit failed
    lscape_list <- readRDS(file.path(dirs[[if (ok[["alfak"]]) "alfak" else "karst"]], "landscape.Rds"))
    k <- lscape_list$k
    f_true <- grf_fitness(k, info$w, sim0$true_landscape)
    key <- data.frame(sim_id = id, w = info$w, pmis = info$pmis, beta = pm_fit, abmrep = info$abmrep, design = design,
                      min_obs = mo, stringsAsFactors = FALSE)
    # each method's landscape accuracy and cross-validation score
    for (m in methods[ok]) {
      l <- readRDS(file.path(dirs[[m]], "landscape.Rds"))
      fm <- if (identical(l$k, k)) f_true else grf_fitness(l$k, info$w, sim0$true_landscape)
      acc <- l
      if (!is.null(l$dist)) {
        acc$fq <- l$dist == 0L
        acc$nn <- l$dist == 1L
      }
      sc <- landscape_scores(acc, fm, readRDS(file.path(dirs[[m]], "xval.Rds")))
      fits[[length(fits) + 1L]] <- cbind(key, method = m, same_list = identical(l$k, k), n_k = nrow(l), sc)
    }
    # the true landscape's forecast on the list (missegregation simulations:
    # shared by the two arms, which fit the same subsampled inputs, unless
    # their lists differ)
    lt <- data.frame(k = k, mean = f_true, stringsAsFactors = FALSE)
    fc <- list()
    for (ta in names(truth_pm)) {
      tdir <- file.path(out_dir, "truth", id, design, paste0("minobs_", mo))
      if (source_set == "misseg") {
        tdir <- file.path(root, "truth", id, design, paste0("minobs_", mo), paste0("pm_", rate_label(truth_pm[[ta]])))
        kfile <- file.path(tdir, "k.Rds")
        if (file.exists(kfile) && !identical(readRDS(kfile), k)) {
          tdir <- file.path(out_dir, "truth", id, design, paste0("minobs_", mo), paste0("pm_", rate_label(truth_pm[[ta]])))
          kfile <- file.path(tdir, "k.Rds")
        }
        dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(kfile)) save_rds_atomic(k, kfile)
      }
      dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
      t_start <- proc.time()[["elapsed"]]
      # a forecast that cannot run (e.g. no parent cell on the list) is left
      # out, as a method's failed forecast is
      done <- tryCatch({
        for (r in seq_len(st$reps)) {
          path <- file.path(tdir, sprintf("rep_%02d.Rds", r))
          if (!file.exists(path)) save_rds_atomic(forecast_abm(lt, x_fit, truth_pm[[ta]], r, settings = st), path)
        }
        TRUE
      }, error = function(e) {
        cat("truth forecast", ta, design, "minobs", mo, "failed:", conditionMessage(e), "\n")
        FALSE
      })
      truth_seconds <- truth_seconds + proc.time()[["elapsed"]] - t_start
      fc[[ta]] <- if (done) mean_forecast(sort(list.files(tdir, pattern = "^rep_", full.names = TRUE))) else NULL
    }
    for (m in methods) {
      fc[[m]] <- if (ok[[m]] && forecast_complete(dirs[[m]])) {
        mean_forecast(sort(list.files(file.path(dirs[[m]], "abm_preds"), full.names = TRUE)))
      } else NULL
    }
    pk <- on_list(p_fit, k)
    for (j in seq_along(td)) {
      # predicted change: parent (the sample the forecast starts from) to the forecast
      for (a in names(fc)) {
        pr <- fc[[a]]
        val <- NA_real_
        if (!is.null(pr) && nrow(pr$mat) > j) {
          fj <- pr$mat[j + 1L, ]
          kk <- colnames(pr$mat)
          val <- wdist(names(pk), pk, kk, fj)
        }
        preds[[length(preds) + 1L]] <- cbind(key, j = j, dt = td[j] - t0, arm = a, pred = val)
      }
      pf <- on_list(p_full, k)
      df <- on_list(d_full[[j]], k)
      dk <- on_list(d_fit[[j]], k)
      obs[[length(obs) + 1L]] <- cbind(key, j = j, dt = td[j] - t0, n_k = length(k),
                                       cells_parent_fit = sum(p_fit), cells_daughter_fit = sum(d_fit[[j]]),
                                       in_list_parent_fit = sum(pk) / sum(p_fit),
                                       in_list_daughter_fit = sum(dk) / sum(d_fit[[j]]),
                                       obs_1000 = w_counts(pf, df), noise_1000 = (halves(pf) + halves(df)) / 2,
                                       obs_fit = w_counts(pk, dk), noise_fit = (halves(pk) + halves(dk)) / 2)
      if (step == 1L) {
        for (a in seq_along(n_sub)) {
          vals <- vapply(draws[[j]][[a]], function(s) w_counts(on_list(s$p, k), on_list(s$d, k)), numeric(1))
          s1 <- draws[[j]][[a]][[1]]
          nz <- (halves(on_list(s1$p, k)) + halves(on_list(s1$d, k))) / 2
          drw[[length(drw) + 1L]] <- cbind(key, j = j, n = n_sub[a], draw = seq_len(n_draws), obs = vals,
                                           noise = nz)
        }
      }
    }
  }
  res <- list(preds = do.call(rbind, preds), obs = do.call(rbind, obs), draws = do.call(rbind, drw),
              fits = do.call(rbind, fits), step = step, t0 = t0, td = td,
              seconds = proc.time()[["elapsed"]] - started, truth_seconds = truth_seconds,
              wdist_fallbacks = wdist_fallbacks)
  save_rds_atomic(res, file.path(out_dir, "rows", paste0(id, ".Rds")))
  cat(id, "done:", NROW(res$preds), "prediction rows,", NROW(res$obs), "observation rows in",
      round(res$seconds), "s (truth forecasts", round(truth_seconds), "s; revsimplex fallbacks",
      wdist_fallbacks, ")\n")
}
