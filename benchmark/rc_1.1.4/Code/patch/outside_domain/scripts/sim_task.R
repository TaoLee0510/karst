# Steps 1-2 of PLAN.md for one simulated input (one simulation and design):
# fitness of the karyotypes outside KARST's domain against the truth (A1), and
# deterministic forecasts on list E from the complete fitted parent to the
# second sampling time after it (A2), for the true landscape, KARST as it is
# and the candidate rules (B1, B2, B1B2), native and level-aligned.
#   Rscript scripts/sim_task.R TASK_FILE ID OUT_DIR
# TASK_FILE columns: source (abm | misseg), arm (step2 | beta_fixed |
# beta_true), sim_id, design (w2_ends | w4_ends | w8_ends).
# Runs with karst (the round's library) on the library path; reads the
# round's fits and the paper's simulations read only. Writes
# OUT_DIR/<source>_<arm>/<sim>_<design>.rds.
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(here, "R", "common.R")); source(file.path(here, "R", "engine.R")); source(file.path(here, "R", "variants.R"))
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
out_file <- file.path(args[[3]], paste0(task$source, "_", task$arm), paste0(task$sim_id, "_", task$design, ".rds"))
if (file.exists(out_file)) quit(save = "no")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
started <- proc.time()[["elapsed"]]
k <- od_karst(use_package = TRUE)
source(file.path(here, "..", "rc_paths.R"))  # rc freeze: roots from RC_HOME / RC_UPSTREAM
paper_data <- rc_paper_data
# the staged benchmark helpers of the Fig. 2d studies (Salehi-sized subsample, seeds, windows), read only
source(rc_abm_R)
cells <- utils::read.delim(file.path(rc_fig2d_truth, "salehi_cells.tsv"))$cells
raw_dir <- file.path(paper_data, "raw", if (task$source == "abm") "ABM" else "ABM_varyMisseg")
fit_path <- if (task$source == "abm") file.path(rc_fig2d_truth, "step2", "fits", "abm", task$sim_id, task$design, "karst", "minobs_10", "karst_fit.Rds") else
  file.path(rc_fig2d_misseg, task$arm, "fits", "misseg", task$sim_id, task$design, "karst", "minobs_10", "karst_fit.Rds")
save <- function(x) { tmp <- paste0(out_file, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, out_file) }
meta <- list(task = task, fit_path = fit_path, karst = as.character(utils::packageVersion("karst")))
if (!file.exists(fit_path)) { save(c(meta, list(status = "no fit"))); quit(save = "no") }
fit <- readRDS(fit_path)
if (!inherits(fit, "karst_sparse_joint_fit")) { save(c(meta, list(status = "fit failed"))); quit(save = "no") }

info <- abm_sim_info(task$sim_id)
sim <- readRDS(file.path(raw_dir, paste0(task$sim_id, ".Rds")))
x_full <- as.matrix(sim$abm_output$x)
xs <- as.matrix(abm_subsample_x(sim$abm_output$x, task$sim_id, cells))
times <- as.numeric(colnames(x_full))
window <- c(w2_ends = 2L, w4_ends = 4L, w8_ends = 8L)[[task$design]]
tt <- abm_window_times(times, window, TRUE, 120)
train <- xs[, match(tt, times), drop = FALSE]; train <- train[rowSums(train) > 0, , drop = FALSE]
# the fit's counts must be this input
fc <- as.matrix(fit$data$counts)
same_input <- isTRUE(all.equal(unname(colSums(fc)), unname(colSums(train)))) && setequal(rownames(fc)[rowSums(fc) > 0], rownames(train))
t0 <- tt[2]; later <- utils::head(times[times > t0], 3); t2 <- later[2]; horizon <- t2 - t0
parent <- train[, 2]; parent <- parent[parent > 0]
daughter <- stats::setNames(as.numeric(x_full[, times == t2]), rownames(x_full)); daughter <- daughter[daughter > 0]
beta_fit <- fit$data$beta

## list E, domain, distances ----------------------------------------------------
canon <- function(v) k$format_karyotypes(k$parse_karyotypes(v))
observed <- unique(canon(rownames(train)))
nb <- function(labels) {
  cn <- od_parse(labels); out <- character(0)
  for (i in seq_len(ncol(cn))) for (s in c(-1L, 1L)) {
    m <- cn; m[, i] <- m[, i] + s; ok <- m[, i] >= 1L
    out <- c(out, do.call(paste, c(as.data.frame(m[ok, , drop = FALSE]), sep = ".")))
  }
  unique(out)
}
one <- setdiff(nb(observed), observed); two <- setdiff(nb(one), c(observed, one))
E <- c(observed, one, two)
cnE <- od_parse(E)
dom <- fit$summary$karyotype
dmajor <- rep(0, length(E)); outE <- !E %in% dom
dmajor[outE] <- k$sparse_joint_major_distance(cnE[outE, , drop = FALSE], fit$trend$major)
eng <- engine_build(E)
n0 <- numeric(length(E)); n0[match(names(parent), E)] <- parent
stopifnot(!anyNA(match(names(parent), E)))

## rates ----------------------------------------------------------------------------
grf <- function(cn, w, centroids) {
  s <- 1 / (pi * sqrt(nrow(centroids))); acc <- numeric(nrow(cn))
  for (c in seq_len(nrow(centroids))) acc <- acc + sin(sqrt(rowSums(sweep(cn, 2L, centroids[c, ], "-")^2)) / w)
  acc * s
}
truth_rate <- grf(cnE, info$w, sim$true_landscape)
self <- 2 * k$sparse_joint_self_probability(cnE, beta_fit) - 1
truth_net <- truth_rate * self
s <- fit$summary
total <- rowSums(as.matrix(s[grep("^count_", names(s))]))
wm <- ifelse(s$major, total, 0); im <- match(s$karyotype, E)
stopifnot(!anyNA(im[s$major]))
truth_ref <- sum(wm[s$major] * truth_net[im[s$major]]) / sum(wm[s$major])
base <- k$predict.karst_sparse_joint_fit(fit, E)
obs_out <- E %in% fit$dropped
anchors <- od_anchor_index(eng, obs_out, dmajor)
variants <- c("current", "B1", "B2", "B1B2")
pv <- lapply(stats::setNames(variants, variants), function(v) od_predict_variant(k, fit, E, v, base = base, anchors = anchors))
w_par <- n0 / sum(n0)
align <- function(r) r + (sum(w_par * truth_rate) - sum(w_par * r))

## A1: fitness of out-of-domain karyotypes -------------------------------------------
nbr_idx <- unique(c(eng$plus[obs_out, ], eng$minus[obs_out, ])); nbr_idx <- nbr_idx[!is.na(nbr_idx)]
a1_idx <- c(which(obs_out), setdiff(nbr_idx[outE[nbr_idx] & !obs_out[nbr_idx] & !E[nbr_idx] %in% observed], which(obs_out)))
a1 <- data.frame(karyotype = E[a1_idx], kind = ifelse(obs_out[a1_idx], "observed", "neighbour"), dmajor = dmajor[a1_idx],
                 anchor_distance = anchors$distance[a1_idx], truth_rel = truth_net[a1_idx] - truth_ref, stringsAsFactors = FALSE)
for (v in variants) a1[[paste0("rel_", v)]] <- pv[[v]]$relative_fitness[a1_idx]
# inside the domain, for reference
in_idx <- which(!outE)
a1_inside <- data.frame(dmajor = 0, kind = ifelse(E[in_idx] %in% observed, "inside_observed", "inside_unseen"),
                        truth_rel = truth_net[in_idx] - truth_ref, rel = base$relative_fitness[in_idx])

## A2: forecasts --------------------------------------------------------------------------
steps <- as.integer(round(horizon * 10))
arms <- list(truth = list(rate = truth_rate, p = beta_fit))
if (task$source == "misseg" && task$arm == "beta_fixed" && info$pmis != beta_fit) arms$truth_true_beta <- list(rate = truth_rate, p = info$pmis)
for (v in variants) {
  arms[[paste0(v, "_native")]] <- list(rate = pv[[v]]$division_rate, p = beta_fit)
  arms[[paste0(v, "_aligned")]] <- list(rate = align(pv[[v]]$division_rate), p = beta_fit)
}
c_par <- od_centre(parent); obs_change <- od_centre(daughter) - c_par
dd <- function(labels) { d <- rep(0, length(labels)); o <- !labels %in% dom
  if (any(o)) d[o] <- k$sparse_joint_major_distance(od_parse(labels[o]), fit$trend$major); d }
d_par <- dd(names(parent)); d_dau <- dd(names(daughter))
share_row <- function(freq, dist) unlist(lapply(c(3, 4, 6), function(r) { sh <- od_shares(freq, od_group(dist, r)); stats::setNames(sh, paste0(names(sh), "_r", r)) }))
observed_shares <- rbind(parent = share_row(parent, d_par), daughter = share_row(daughter, d_dau))
fc_rows <- list()
for (a in names(arms)) {
  f <- engine_run(eng, n0, arms[[a]]$rate, arms[[a]]$p, 0.1, steps = steps)[, 1]
  keep <- f > 0
  fr <- f[keep]
  fc_rows[[a]] <- data.frame(arm = a, t(share_row(fr, dmajor[keep])), angle = od_angle(od_centre(fr) - c_par, obs_change),
                             overlap = od_overlap(fr, daughter),
                             w1 = tryCatch(od_w1(fr, daughter, floor = 1e-5), error = function(e) NA_real_),
                             parent_rate_mean = sum(w_par * arms[[a]]$rate), check.names = FALSE)
}
fc <- do.call(rbind, fc_rows)
save(c(meta, list(status = "ok", same_input = same_input, info = info, beta_fit = beta_fit, t0 = t0, t2 = t2, horizon = horizon,
                  n_E = length(E), n_domain = length(dom), n_observed_out = sum(obs_out),
                  parent_cells = sum(parent), daughter_cells = sum(daughter), observed_shares = observed_shares,
                  forecasts = fc, a1 = a1, a1_inside_summary = c(n = nrow(a1_inside), bias_median = stats::median(a1_inside$rel - a1_inside$truth_rel),
                                                                 spearman = suppressWarnings(stats::cor(a1_inside$rel, a1_inside$truth_rel, method = "spearman"))),
                  truth_ref = truth_ref, karst_ref = fit$trend$reference$fitness, seconds = proc.time()[["elapsed"]] - started)))
