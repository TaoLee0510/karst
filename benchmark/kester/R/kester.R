# Kester organoid forecasts (benchmark/rc_1.1.4/PLAN.md, part 3).
# Needs common.R (od_*), engine.R and support.R (fr_list_e, fr_mix).

kester_grid <- c(1, 2, 4, 8, 16, 32, 64)   # multiples of 5e-5 (D1a)
kester_days_per_passage <- 7                # weekly passaging (D4a)
kester_seed <- 20261126L

kester_inputs <- function(root = Sys.getenv("KESTER", "/share/lab_crd/taoli/Project/Kester")) {
  b <- file.path(root, "results")
  data.frame(culture = c("APKS1", "APKS2", "APKS3", "A3", "AP1"),
             kind = c("salehi", "salehi", "salehi", "gbm", "gbm"),
             path = c(file.path(b, "cohort_backbone", paste0(c("APKS1", "APKS2", "APKS3"), "_karst_counts.tsv")),
                      file.path(b, "cohort_backbone_a3_ap1", c("A3_karst_counts.tsv", "AP1_karst_counts.tsv"))),
             stringsAsFactors = FALSE)
}

# Integer count matrix: karyotype rows, numeric passage columns (as the prompt reads it).
kester_read <- function(path) {
  x <- as.matrix(utils::read.delim(path, row.names = 1, check.names = FALSE))
  storage.mode(x) <- "double"
  x[rowSums(x) > 0, order(as.numeric(colnames(x))), drop = FALSE]
}

# Salehi-like folds (training passages -> held-out passage; prompt section 1).
kester_salehi_folds <- data.frame(
  culture = c("APKS1", "APKS2", "APKS2", "APKS3", "APKS3", "APKS3", "APKS2", "APKS3", "APKS3"),
  train = c("10,22", "5,8", "5,8,14", "3,7", "3,7,15", "3,7,15,19", "5,14", "3,15", "3,19"),
  target = c(25, 14, 18, 15, 19, 23, 18, 19, 23),
  design = c(rep("all", 6), rep("ends", 3)), stringsAsFactors = FALSE)

# The cell partitions of kflInfer (PANcanKFLs_karst code/packages/kflInfer/
# Fitting_and_validation.R, karst_validation_folds(), two time points).
kester_cell_folds <- function(x, n_folds = 3L, seed = kester_seed) {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (is.null(old_seed)) { if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }
          else assign(".Random.seed", old_seed, envir = .GlobalEnv))
  set.seed(seed)
  partitions <- lapply(seq_len(ncol(x)), function(j) {
    group <- sample(rep(seq_len(n_folds), length.out = sum(x[, j])))
    as.matrix(table(factor(rep(seq_len(nrow(x)), x[, j]), levels = seq_len(nrow(x))), factor(group, levels = seq_len(n_folds))))
  })
  lapply(seq_len(n_folds), function(i) {
    test <- do.call(cbind, lapply(partitions, function(m) as.numeric(m[, i])))
    dimnames(test) <- dimnames(x)
    list(id = paste0("cells_", i), train = x - test, test = test)
  })
}

# The scoring of kflInfer (karst_score_counts()): per-cell negative log
# likelihood on a fixed support with a half-count pseudocount, against the
# training baseline. `prediction`: named expected composition.
kester_score_counts <- function(observed, prediction, baseline, common) {
  labels <- union(names(observed)[observed > 0], names(prediction))
  y <- as.numeric(observed[labels]); y[is.na(y)] <- 0
  used <- labels %in% common
  raw <- as.numeric(prediction[labels])
  p <- b <- rep(NA_real_, length(labels))
  scored_depth <- sum(y[used])
  if (any(used)) {
    b_counts <- as.numeric(baseline[labels[used]]); b_counts[is.na(b_counts)] <- 0
    n_train <- sum(baseline, na.rm = TRUE)
    pseudo <- 0.5 / sum(used)
    b[used] <- (b_counts + pseudo) / (sum(b_counts) + 0.5)
    p[used] <- (n_train * raw[used] / sum(raw[used]) + pseudo) / (n_train + 0.5)
  }
  loss <- function(q) if (scored_depth > 0) -sum(y[used] * log(q[used])) / scored_depth else NA_real_
  c(total_cells = sum(y), scored_cells = scored_depth, scored_coverage = if (sum(y) > 0) scored_depth / sum(y) else NA_real_,
    nll_per_cell = loss(p), baseline_nll_per_cell = loss(b))
}

# ALFA-K's landscape on list E: its mean on its own list, its stable Kriging
# prediction elsewhere (benchmark/paper/R/align.R, align_alfak_landscape()).
kester_alfak_on <- function(dir, E, k, chunk = 20000L) {
  lsc <- readRDS(file.path(dir, "landscape.Rds"))
  boot <- readRDS(file.path(dir, "bootstrap_res.Rds"))
  fboot <- cbind(boot$final_fitness, boot$nn_fitness)
  strs <- c(colnames(boot$final_fitness), colnames(boot$nn_fitness))
  fmean <- colMeans(fboot, na.rm = TRUE)[strs]; valid <- is.finite(fmean)
  mean <- rep(NA_real_, length(E)); inside <- match(lsc$k, E); keep <- !is.na(inside)
  mean[inside[keep]] <- lsc$mean[keep]
  outside <- which(is.na(mean))
  if (length(outside)) {
    if (sum(valid) < 3L) stop("too few karyotypes for Kriging")
    ktrain <- k$parse_karyotypes(strs[valid]); storage.mode(ktrain) <- "double"
    krig <- suppressWarnings(fields::Krig(ktrain, fmean[valid], cov.function = "stationary.cov",
                                          cov.args = alfakR:::krig_covariance_args(), nstep.cv = alfakR:::ALFAK_KRIG_NSTEP_CV,
                                          give.warnings = FALSE))
    km <- od_parse(E[outside]); storage.mode(km) <- "double"
    idx <- split(seq_along(outside), ceiling(seq_along(outside) / chunk))
    mean[outside] <- unlist(lapply(idx, function(i) as.numeric(stats::predict(krig, km[i, , drop = FALSE]))), use.names = FALSE)
  }
  mean
}

# KARST's native division rate on labels (fitness / (2 P_ii - 1)), as ../outside_rank step 3.
kester_karst_rates <- function(fit, labels, k) {
  k$predict.karst_sparse_joint_fit(fit, labels)$fitness_mean / (2 * k$sparse_joint_self_probability(od_parse(labels), fit$data$beta) - 1)
}

# Scores of a forecast `f` (shares on list E) against the held-out counts `d`
# from the parent counts `par`: full-centroid angle, W1, overlap, and the
# per-cell log score with weight `eps` on a uniform over E.
kester_engine_scores <- function(f, par, d, eps) {
  E <- names(f)
  fr <- f[f > 0]
  q <- (1 - eps) * f / sum(f) + eps / length(E)
  inE <- names(d) %in% E
  c(angle = od_angle(od_centre(fr) - od_centre(par), od_centre(d) - od_centre(par)),
    w1 = tryCatch(od_w1(fr, d, floor = 1e-5), error = function(e) NA_real_), overlap = od_overlap(fr, d),
    log_score = if (any(inE)) sum(d[inE] * log(q[names(d)[inE]])) / sum(d[inE]) else NA_real_,
    target_cells_outside_E = sum(d[!inE]), target_cells = sum(d))
}

# The paper's Fig. 2c-style scores (benchmark/paper/R/metrics_salehi.R,
# salehi_prediction_scores()) of the mean ABM forecast at its last time,
# on the landscape's karyotypes.
kester_paper_scores <- function(preds, par, d) {
  mat <- Reduce(`+`, lapply(preds, function(x) as.matrix(x[, -1]))) / length(preds)
  yp <- as.numeric(mat[nrow(mat), ]); kk <- colnames(mat)
  k <- od_parse(kk)
  yb <- yt <- numeric(length(kk)); names(yb) <- names(yt) <- kk
  pb <- par[names(par) %in% kk]; dt_ <- d[names(d) %in% kk]
  if (!length(pb) || !length(dt_)) return(c(angle = NA, wasserstein = NA, overlap = NA, coverage_target = NA))
  yb[names(pb)] <- pb / sum(pb); yt[names(dt_)] <- dt_ / sum(dt_)
  kt <- colSums(yt * k); kb <- colSums(yb * k); kp <- colSums(yp * k)
  ang <- 180 * acos(max(-1, min(1, sum((kt - kb) * (kp - kb)) / sqrt(sum((kt - kb)^2) * sum((kp - kb)^2))))) / pi
  c(angle = ang, wasserstein = tryCatch(transport::wasserstein(transport::wpp(k, mass = yt + 1e-12), transport::wpp(k, mass = yp + 1e-12), method = "revsimplex"), error = function(e) NA_real_),
    overlap = sum(pmin(yp, yt)) / min(sum(yp), sum(yt)), coverage_target = sum(dt_) / sum(d))
}
