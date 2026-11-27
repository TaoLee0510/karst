# Study B, ALFA-K arms (added at the user's request after the report): the
# same input, list E, engine and scores as scripts/sim_task.R, with ALFA-K's
# landscape (minobs 10) as the rate map: its mean on its own list and its
# stable Kriging prediction elsewhere on E (as benchmark/paper/R/align.R's
# align_alfak_landscape(); the map codex used on Salehi). Rates enter the
# engine as they are, as in ../outside_rank step 3.
#   Rscript scripts/alfak_task.R TASK_FILE ID OUT_DIR
# Arms: alfak_g1 (the beta-fixed fit, forecast x1); misseg only:
# alfak_gtrue (same fit, forecast at the true rate) and alfakT_gtrue (the fit
# at the true beta, forecast at the true rate). Writes OUT_DIR/alfak_<source>/.
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "engine.R")); source(file.path(here, "R", "fr.R"))
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
out_file <- file.path(args[[3]], paste0("alfak_", task$source), paste0(task$sim_id, "_", task$design, ".rds"))
if (file.exists(out_file)) quit(save = "no")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
save <- function(x) { tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, out_file) }
started <- proc.time()[["elapsed"]]
k <- od_karst(use_package = TRUE)
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
paper_data <- rc_paper_data
source(rc_abm_R)
cells <- utils::read.delim(file.path(rc_fig2d_truth, "salehi_cells.tsv"))$cells
info <- abm_sim_info(task$sim_id); base_rate <- 5e-5; true_mult <- info$pmis / base_rate
fit_dir <- function(arm, method) if (task$source == "abm")
  file.path(rc_fig2d_truth, "step2", "fits", "abm", task$sim_id, task$design, method, "minobs_10") else
  file.path(rc_fig2d_misseg, arm, "fits", "misseg", task$sim_id, task$design, method, "minobs_10")
meta <- list(task = task, info = info, true_mult = true_mult)
a_fixed <- fit_dir("beta_fixed", "alfak"); a_true <- fit_dir("beta_true", "alfak")
k1 <- file.path(fit_dir("beta_fixed", "karst"), "karst_fit.Rds")
if (!file.exists(file.path(a_fixed, "bootstrap_res.Rds")) || !file.exists(k1)) { save(c(meta, list(status = "no fit"))); quit(save = "no") }

## input, parent, daughter, list E: as scripts/sim_task.R -----------------------------
raw_dir <- file.path(paper_data, "raw", if (task$source == "abm") "ABM" else "ABM_varyMisseg")
sim <- readRDS(file.path(raw_dir, paste0(task$sim_id, ".Rds")))
x_full <- as.matrix(sim$abm_output$x)
xs <- as.matrix(abm_subsample_x(sim$abm_output$x, task$sim_id, cells))
times <- as.numeric(colnames(x_full))
window <- c(w2_ends = 2L, w4_ends = 4L, w8_ends = 8L)[[task$design]]
tt <- abm_window_times(times, window, TRUE, 120)
train <- xs[, match(tt, times), drop = FALSE]; train <- train[rowSums(train) > 0, , drop = FALSE]
t0 <- tt[2]; later <- utils::head(times[times > t0], 3); t2 <- later[2]; horizon <- t2 - t0
parent <- train[, 2]; parent <- parent[parent > 0]
daughter <- stats::setNames(as.numeric(x_full[, times == t2]), rownames(x_full)); daughter <- daughter[daughter > 0]
canon <- function(v) k$format_karyotypes(k$parse_karyotypes(v))
observed <- unique(canon(rownames(train)))
daughter <- tapply(daughter, canon(names(daughter)), sum); daughter <- stats::setNames(as.numeric(daughter), names(daughter))
E <- fr_list_e(observed); cnE <- od_parse(E)
fit1 <- readRDS(k1)  # only for the domain and majors that group the scores, as in sim_task.R
dom <- fit1$summary$karyotype
dmajor <- rep(0, length(E)); outE <- !E %in% dom
dmajor[outE] <- k$sparse_joint_major_distance(cnE[outE, , drop = FALSE], fit1$trend$major)
eng <- engine_build(E)
n0 <- numeric(length(E)); n0[match(names(parent), E)] <- parent
steps <- as.integer(round(horizon * 10))

## ALFA-K's landscape on E ---------------------------------------------------------------
alfak_on_e <- function(dir, chunk = 20000L) {
  lsc <- readRDS(file.path(dir, "landscape.Rds"))
  boot <- readRDS(file.path(dir, "bootstrap_res.Rds"))
  fboot <- cbind(boot$final_fitness, boot$nn_fitness)
  strs <- c(colnames(boot$final_fitness), colnames(boot$nn_fitness))
  fmean <- colMeans(fboot, na.rm = TRUE)[strs]; valid <- is.finite(fmean)
  ktrain <- k$parse_karyotypes(strs); storage.mode(ktrain) <- "double"
  krig <- suppressWarnings(fields::Krig(ktrain[valid, , drop = FALSE], fmean[valid], cov.function = "stationary.cov",
                                        cov.args = alfakR:::krig_covariance_args(), nstep.cv = alfakR:::ALFAK_KRIG_NSTEP_CV,
                                        give.warnings = FALSE))
  mean <- rep(NA_real_, length(E)); inside <- match(lsc$k, E); keep <- !is.na(inside)
  mean[inside[keep]] <- lsc$mean[keep]
  outside <- which(is.na(mean))
  if (length(outside)) {
    km <- cnE[outside, , drop = FALSE]; storage.mode(km) <- "double"
    idx <- split(seq_along(outside), ceiling(seq_along(outside) / chunk))
    mean[outside] <- unlist(lapply(idx, function(i) as.numeric(stats::predict(krig, km[i, , drop = FALSE]))), use.names = FALSE)
  }
  list(rate = mean, own_list = sum(keep), own_list_outside_E = sum(!keep), kriged = length(outside))
}
run <- function(r, g) engine_run(eng, n0, r, base_rate * g, 0.1, steps = steps)[, 1]
maps <- list(fixed = alfak_on_e(a_fixed))
fc <- list(alfak_g1 = run(maps$fixed$rate, 1))
if (task$source == "misseg" && true_mult > 1) {
  fc$alfak_gtrue <- run(maps$fixed$rate, true_mult)
  if (file.exists(file.path(a_true, "bootstrap_res.Rds"))) { maps$true <- alfak_on_e(a_true); fc$alfakT_gtrue <- run(maps$true$rate, true_mult) }
}
dist_of <- function(labels) { d <- rep(0, length(labels)); o <- !labels %in% dom
  if (any(o)) d[o] <- k$sparse_joint_major_distance(od_parse(labels[o]), fit1$trend$major); d }
rows <- fr_score(fc, parent, daughter, observed, dmajor, dist_of)
save(c(meta, list(status = "ok", n_E = length(E), maps = lapply(maps, function(m) m[c("own_list", "own_list_outside_E", "kriged")]),
                  rows = rows, seconds = proc.time()[["elapsed"]] - started)))
