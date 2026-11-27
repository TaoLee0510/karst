# Study B (PLAN.md), one simulated input: deterministic forecasts from KARST
# fits at fitted rate f and forecast rate g (multiples of 5e-5), the mixture
# over g, and the true landscape, scored against the simulation's sample at
# the second sampling time after the fitted parent (as ../outside_domain).
#   Rscript scripts/sim_task.R TASK_FILE ID OUT_DIR
# TASK_FILE columns: source (abm | misseg), sim_id, design, rep, half.
# Fits (read only): x1 = the round's fits (fig2d_truth step2 / fig2d_misseg
# beta_fixed); x2-x16 (misseg only) = the dispersion-gate study's fits.
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "engine.R")); source(file.path(here, "R", "fr.R"))
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
out_file <- file.path(args[[3]], paste0("sim_", task$source), paste0(task$sim_id, "_", task$design, ".rds"))
if (file.exists(out_file)) quit(save = "no")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
save <- function(x) { tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, out_file) }
started <- proc.time()[["elapsed"]]
k <- od_karst(use_package = TRUE)
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
VM <- rc_gate_vm
paper_data <- rc_paper_data
source(rc_abm_R)
cells <- utils::read.delim(file.path(rc_fig2d_truth, "salehi_cells.tsv"))$cells
info <- abm_sim_info(task$sim_id)
grid <- c(1, 2, 4, 8, 16); base_rate <- 5e-5
fit_path <- function(m) {
  if (m == 1) {
    if (task$source == "abm") file.path(rc_fig2d_truth, "step2", "fits", "abm", task$sim_id, task$design, "karst", "minobs_10", "karst_fit.Rds") else
      file.path(rc_fig2d_misseg, "beta_fixed", "fits", "misseg", task$sim_id, task$design, "karst", "minobs_10", "karst_fit.Rds")
  } else file.path(VM, sprintf("x%02d", m), "fits", "misseg", task$sim_id, task$design, "karst", "minobs_10", "karst_fit.Rds")
}
meta <- list(task = task, info = info, karst = as.character(utils::packageVersion("karst")))
fits <- list()
for (m in if (task$source == "abm") 1 else grid) {
  p <- fit_path(m)
  if (file.exists(p)) { f <- readRDS(p); if (inherits(f, "karst_sparse_joint_fit")) fits[[as.character(m)]] <- f }
}
if (is.null(fits[["1"]])) { save(c(meta, list(status = "no x1 fit"))); quit(save = "no") }

## input, parent, daughter, list E (as ../outside_domain/scripts/sim_task.R) -----------
raw_dir <- file.path(paper_data, "raw", if (task$source == "abm") "ABM" else "ABM_varyMisseg")
sim <- readRDS(file.path(raw_dir, paste0(task$sim_id, ".Rds")))
x_full <- as.matrix(sim$abm_output$x)
xs <- as.matrix(abm_subsample_x(sim$abm_output$x, task$sim_id, cells))
times <- as.numeric(colnames(x_full))
window <- c(w2_ends = 2L, w4_ends = 4L, w8_ends = 8L)[[task$design]]
tt <- abm_window_times(times, window, TRUE, 120)
train <- xs[, match(tt, times), drop = FALSE]; train <- train[rowSums(train) > 0, , drop = FALSE]
same_input <- vapply(fits, function(f) { fc <- as.matrix(f$data$counts)
  isTRUE(all.equal(unname(colSums(fc)), unname(colSums(train)))) && setequal(rownames(fc)[rowSums(fc) > 0], rownames(train)) }, logical(1))
t0 <- tt[2]; later <- utils::head(times[times > t0], 3); t2 <- later[2]; horizon <- t2 - t0
parent <- train[, 2]; parent <- parent[parent > 0]
daughter <- stats::setNames(as.numeric(x_full[, times == t2]), rownames(x_full)); daughter <- daughter[daughter > 0]
canon <- function(v) k$format_karyotypes(k$parse_karyotypes(v))
observed <- unique(canon(rownames(train)))
daughter <- tapply(daughter, canon(names(daughter)), sum); daughter <- stats::setNames(as.numeric(daughter), names(daughter))
E <- fr_list_e(observed)
cnE <- od_parse(E)
fit1 <- fits[["1"]]
dom <- fit1$summary$karyotype
dmajor <- rep(0, length(E)); outE <- !E %in% dom
dmajor[outE] <- k$sparse_joint_major_distance(cnE[outE, , drop = FALSE], fit1$trend$major)
eng <- engine_build(E)
n0 <- numeric(length(E)); n0[match(names(parent), E)] <- parent
stopifnot(!anyNA(match(names(parent), E)))
steps <- as.integer(round(horizon * 10))

## division rates ----------------------------------------------------------------------
div_rate <- function(f) {
  pr <- k$predict.karst_sparse_joint_fit(f, E)
  pr$fitness_mean / (2 * k$sparse_joint_self_probability(cnE, f$data$beta) - 1)
}
rates <- lapply(fits, div_rate)
grf <- function(cn, w, centroids) {
  s <- 1 / (pi * sqrt(nrow(centroids))); acc <- numeric(nrow(cn))
  for (c in seq_len(nrow(centroids))) acc <- acc + sin(sqrt(rowSums(sweep(cn, 2L, centroids[c, ], "-")^2)) / w)
  acc * s
}
truth_rate <- grf(cnE, info$w, sim$true_landscape)
w_par <- n0 / sum(n0)
align <- function(r) r + (sum(w_par * truth_rate) - sum(w_par * r))

## forecasts -----------------------------------------------------------------------------
run <- function(r, g) engine_run(eng, n0, r, base_rate * g, 0.1, steps = steps)[, 1]
fc <- list()
for (g in grid) fc[[sprintf("f1_g%d", g)]] <- run(rates[["1"]], g)
for (g in grid) fc[[sprintf("f1a_g%d", g)]] <- run(align(rates[["1"]]), g)
for (m in setdiff(names(rates), "1")) fc[[sprintf("f%s_g%s", m, m)]] <- run(rates[[m]], as.numeric(m))
fc[["truth_g1"]] <- run(truth_rate, 1)
true_mult <- info$pmis / base_rate
if (!isTRUE(all.equal(true_mult, 1))) fc[["truth_gtrue"]] <- run(truth_rate, true_mult)
for (top in c(4, 8, 16)) fc[[sprintf("mix_1_%d", top)]] <- fr_mix(fc[sprintf("f1_g%d", grid[grid <= top])])

## scores ---------------------------------------------------------------------------------
dist_of <- function(labels) { d <- rep(0, length(labels)); o <- !labels %in% dom
  if (any(o)) d[o] <- k$sparse_joint_major_distance(od_parse(labels[o]), fit1$trend$major); d }
rows <- fr_score(fc, parent, daughter, observed, dmajor, dist_of)
save(c(meta, list(status = "ok", same_input = same_input, fits_found = names(fits), beta = vapply(fits, function(f) f$data$beta, numeric(1)),
                  true_mult = true_mult, t0 = t0, t2 = t2, horizon = horizon, n_E = length(E),
                  parent_cells = sum(parent), daughter_cells = sum(daughter), rows = rows,
                  seconds = proc.time()[["elapsed"]] - started)))
