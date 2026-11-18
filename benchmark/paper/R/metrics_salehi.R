# Scores of a fitted Salehi trajectory, as 03_process_alfak_fits.R of the
# paper computes them: for each daughter passage of the trajectory's last
# passage, every forecast day is compared with the daughter's observed
# karyotype distribution (prediction) and the parent's (baseline).

s2v <- function(x) as.numeric(unlist(strsplit(x, ".", fixed = TRUE)))

# One daughter passage `dj` (uid); NULL when it cannot be scored.
salehi_prediction_scores <- function(preds, dj, meta, inputs, inputs_dir, treat_lut) {
  mat <- Reduce(`+`, lapply(preds, function(x) as.matrix(x[, -1]))) / length(preds)
  preds <- as.data.frame(mat)
  fj <- paste(meta[meta$uid == dj, c("datasetname", "timepoint")], collapse = "_")
  if (sum(grepl(fj, inputs)) < 1) return(NULL)
  fj <- inputs[grepl(fj, inputs)][1]
  k <- do.call(rbind, lapply(colnames(preds), s2v))
  tru <- readRDS(file.path(inputs_dir, fj))$x
  preds <- preds[-1, ]
  yt <- 0 * preds[1, ]
  yb <- 0 * preds[1, ]
  tru <- tru[rownames(tru) %in% colnames(preds), (ncol(tru) - 1):ncol(tru)]
  if (ncol(tru) != 2) stop("expected 2 columns here")
  yb[rownames(tru)] <- tru[, 1] / sum(tru[, 1])
  yt[rownames(tru)] <- tru[, 2] / sum(tru[, 2])
  yt <- as.numeric(yt)
  yb <- as.numeric(yb)
  if (!is.finite(sum(yt)) || !is.finite(sum(yb))) return(NULL)
  mag <- function(v) sqrt(sum(v^2))
  getangle <- function(a, b) 180 * acos(sum(a * b) / (mag(a) * mag(b))) / pi
  cosine <- function(a, b) sum(a * b) / prod(sqrt(sum(a^2)), sqrt(sum(b^2)))
  overlap <- function(a, b) sum(pmin(a, b)) / min(sum(a), sum(b))
  entropy <- function(p) {
    p <- p[p > 0] / sum(p)
    -sum(p * log(p))
  }
  kt <- colSums(yt * k)
  kb <- colSums(yb * k)
  st <- transport::wpp(k, mass = yt)
  sb <- transport::wpp(k, mass = yb)
  # Beyond 03, for the coverage analysis: a "denoised" baseline, the parent
  # passage restricted to its frequent karyotypes (share >= 1% and >= 3
  # cells, karst's major rule applied to one time point) and renormalised, and
  # the entropies of the forecast, parent and daughter distributions.
  parent_cells <- stats::setNames(numeric(ncol(preds)), colnames(preds))
  parent_cells[rownames(tru)] <- tru[, 1]
  frequent <- yb >= 0.01 & parent_cells >= 3
  if (!any(frequent)) frequent <- yb > 0
  yd <- yb * frequent / sum(yb * frequent)
  kd <- colSums(yd * k)
  denoised <- c(euclidean = sqrt(sum((kd - kt)^2)), wasserstein = transport::wasserstein(st, transport::wpp(k, mass = yd)),
                cosine = cosine(yd, yt), overlap = overlap(yd, yt))
  df <- do.call(rbind, lapply(seq_len(nrow(preds)), function(tt) {
    yp <- as.numeric(preds[tt, ])
    if (!is.finite(sum(yp)) || sum(yt) <= 0 || min(yp) < 0) return(NULL)
    kp <- colSums(yp * k)
    sp <- transport::wpp(k, mass = yp)
    out <- rbind(
      data.frame(pred0 = sqrt(sum((kp - kb)^2)), pred = sqrt(sum((kp - kt)^2)), base = sqrt(sum((kb - kt)^2)),
                 metric = "euclidean"),
      data.frame(pred0 = transport::wasserstein(sb, sp), pred = transport::wasserstein(st, sp),
                 base = transport::wasserstein(st, sb), metric = "wasserstein"),
      data.frame(pred0 = cosine(yp, yb), pred = cosine(yp, yt), base = cosine(yb, yt), metric = "cosine"),
      data.frame(pred0 = overlap(yp, yb), pred = overlap(yp, yt), base = overlap(yb, yt), metric = "overlap"))
    out$tt <- tt
    out$angle <- getangle(kt - kb, kp - kb)
    out$denoised <- unname(denoised[out$metric])
    out$entropy_pred <- entropy(yp)
    out
  }))
  if (is.null(df)) return(NULL)
  df$test_treat <- as.character(treat_lut[paste0("u", dj)])
  df$entropy_base <- entropy(yb)
  df$entropy_true <- entropy(yt)
  df
}

## References of Fig. 2c, g-h that do not depend on the method (04) ----------

# Forked sister passages: two daughters i, j of one parent, each with its
# 2-passage input (parent, daughter). The daughters' distributions are
# compared with each other (value) and with the parent (baseline1, baseline2);
# "angle" is the angle between the two parent-to-daughter changes.
salehi_sister_metrics <- function(x1, x2) {
  clones <- union(rownames(x1), rownames(x2))
  get_counts <- function(mat, col) {
    v <- numeric(length(clones))
    names(v) <- clones
    present <- intersect(rownames(mat), clones)
    v[present] <- mat[present, col]
    v
  }
  v11 <- get_counts(x1, 1); v12 <- get_counts(x1, 2)
  v21 <- get_counts(x2, 1); v22 <- get_counts(x2, 2)
  k <- do.call(rbind, lapply(clones, s2v))
  vec <- function(y) colSums(k * y)
  mag <- function(v) sqrt(sum(v^2))
  getangle <- function(a, b) 180 * acos(sum(a * b) / (mag(a) * mag(b))) / pi
  overlap <- function(a, b) { ya <- a / sum(a); yb <- b / sum(b); sum(pmin(ya, yb)) / min(sum(ya), sum(yb)) }
  euclid <- function(a, b) { ya <- a / sum(a); yb <- b / sum(b); sqrt(sum((vec(ya) - vec(yb))^2)) }
  cosine <- function(a, b) { ya <- a / sum(a); yb <- b / sum(b); sum(ya * yb) / (sqrt(sum(ya^2)) * sqrt(sum(yb^2))) }
  wass <- function(a, b) {
    transport::wasserstein(transport::wpp(k, mass = a / sum(a)), transport::wpp(k, mass = b / sum(b)))
  }
  funs <- list(overlap = overlap, euclidean = euclid, cosine = cosine, wasserstein = wass)
  df <- do.call(rbind, lapply(names(funs), function(mn) {
    data.frame(metric = mn, value = funs[[mn]](v12, v22), baseline1 = funs[[mn]](v11, v12),
               baseline2 = funs[[mn]](v21, v22), stringsAsFactors = FALSE)
  }))
  delta1 <- vec(v12 / sum(v12)) - vec(v11 / sum(v11))
  delta2 <- vec(v22 / sum(v22)) - vec(v21 / sum(v21))
  rbind(df, data.frame(metric = "angle", value = getangle(delta1, delta2), baseline1 = NA, baseline2 = NA,
                       stringsAsFactors = FALSE))
}

# Angle between the first-to-second-passage changes of two unrelated inputs.
salehi_pair_angle <- function(a, b) {
  clones <- union(rownames(a), rownames(b))
  get_col <- function(mat, j) {
    v <- numeric(length(clones))
    names(v) <- clones
    v[match(rownames(mat), clones)] <- mat[, j]
    v
  }
  c11 <- get_col(a, 1); c12 <- get_col(a, 2)
  c21 <- get_col(b, 1); c22 <- get_col(b, 2)
  kmat <- do.call(rbind, lapply(clones, s2v))
  d1 <- colSums(kmat * (c12 / sum(c12) - c11 / sum(c11)))
  d2 <- colSums(kmat * (c22 / sum(c22) - c21 / sum(c21)))
  acos(sum(d1 * d2) / (sqrt(sum(d1^2)) * sqrt(sum(d2^2)))) * 180 / pi
}

# Lineage label of a time point or trajectory id (utils_lineage.R).
salehi_lineage_map <- c(SA1035T = "SA1035", SA1035U = "SA1035", SA532 = "SA532", SA609 = "SA609",
                        SA609UnBU = "SA609", SA906a = "p53 k.o", SA906b = "p53 k.o", SA039U = "p53 w.t",
                        SA535_CISPLATIN_CombinedT = "SA535", SA535_CISPLATIN_CombinedU = "SA535", SA000 = "SA609",
                        SA609R2 = "SA609", SA535_CISPLATIN_CombinedH = "SA535", SA001 = "SA609", SA609H2 = "SA609",
                        SA609H1 = "SA609", SA609R2b = "SA609", SA609TIIDT = "SA609")
salehi_assign_lineage <- function(fi, meta) {
  lut <- stats::setNames(meta$datasetname, paste(meta$datasetname, meta$timepoint, sep = "_"))
  id <- vapply(strsplit(fi, "_l_", fixed = TRUE), `[`, character(1), 1)
  unname(salehi_lineage_map[lut[id]])
}
# Fig. 2 (04) uses the first 12 entries only; trajectories of the other
# time points get no lineage there and drop out of the per-lineage panels.
salehi_lineage_fig2 <- function(fi, meta) {
  lut <- stats::setNames(meta$datasetname, paste(meta$datasetname, meta$timepoint, sep = "_"))
  pdx <- vapply(strsplit(fi, "_", fixed = TRUE), function(p) paste(p[1:(length(p) - 6)], collapse = "_"),
                character(1))
  unname(salehi_lineage_map[1:12][lut[pdx]])
}

# Is the median of the cluster medians of the angles below the 90 degrees of
# random directions in N dimensions? Monte Carlo, one-sided (utils_stats.R).
r_sphere_angle <- function(n, N = 22) acos(2 * stats::rbeta(n, (N - 1) / 2, (N - 1) / 2) - 1) * 180 / pi
sphere_null_test <- function(angles, cluster, N = 22, B = 1e4) {
  stopifnot(length(angles) == length(cluster))
  cl_sizes <- table(cluster)
  obs <- stats::median(tapply(angles, cluster, stats::median))
  null <- replicate(B, stats::median(vapply(cl_sizes, function(n) stats::median(r_sphere_angle(n, N)), numeric(1))))
  list(stat = obs, p = mean(null <= obs))
}

## Novel karyotypes, Fig. 3 (05) ------------------------------------------------

# 2-passage inputs (parent, daughter) of the daughters of a trajectory.
salehi_daughter_inputs <- function(traj_id, lineages, meta) {
  lin <- lineages[[traj_id]]
  if (length(lin$dec1) == 0) return(character(0))
  patterns <- paste(meta$datasetname[meta$uid %in% lin$dec1], meta$timepoint[meta$uid %in% lin$dec1], "l", 2,
                    sep = "_")
  as.character(unlist(lapply(patterns, function(p) names(lineages)[grepl(p, names(lineages))])))
}

# Rows for the regression of 05: for every landscape karyotype absent from the
# daughter's parent passage, the fraction of the parent population at
# distances 1-5, its fitness f, and whether it newly appears in the daughter.
salehi_novel_data <- function(lscape, x) {
  last_col <- ncol(x)
  prev_col <- last_col - 1
  present <- x[, last_col] > 0
  newk <- rownames(x)[present & (x[, prev_col] == 0)]
  oldk <- rownames(x)[present]
  oldv <- do.call(rbind, lapply(oldk, s2v))
  vn <- do.call(rbind, lapply(lscape$k, s2v))
  N <- nrow(vn)
  M <- nrow(oldv)
  # Manhattan distances between landscape and present karyotypes only (05
  # takes this block of the full distance matrix, which needs tens of GB for
  # large landscapes; the values are the same).
  tv <- t(vn)
  d_mat <- matrix(vapply(seq_len(M), function(j) colSums(abs(tv - oldv[j, ])), numeric(N)), nrow = N)
  prev_counts <- x[oldk, prev_col]
  result <- matrix(0, nrow = N, ncol = 6)
  for (d in 0:5) result[, d + 1] <- (d_mat == d) %*% prev_counts
  result <- result / sum(x[, prev_col])
  df <- as.data.frame(result)
  names(df) <- paste0("d", 0:5)
  df$f <- lscape$mean
  df$y <- lscape$k %in% newk
  df <- df[df$d0 == 0, -1, drop = FALSE]
  rownames(df) <- NULL
  df
}

# All rows of 03 for one fit directory (one trajectory, design, method and
# minobs): the fit's CV score and the forecast scores of every daughter.
salehi_fit_scores <- function(fit_dir, traj_id, lineages, meta, inputs_dir) {
  treat_lut <- stats::setNames(meta$on_treatment, paste0("u", meta$uid))
  li <- lineages[[traj_id]]
  # This benchmark's fits have run.Rds; for the paper's, a fit without
  # xval.Rds failed (as in 03).
  run_path <- file.path(fit_dir, "run.Rds")
  ok <- if (file.exists(run_path)) identical(readRDS(run_path)$status, "ok") else
    file.exists(file.path(fit_dir, "xval.Rds"))
  xv <- if (ok) xval_score(readRDS(file.path(fit_dir, "xval.Rds"))) else NaN
  head <- data.frame(fi = traj_id, min_obs = sub("minobs_", "", basename(fit_dir)), ntrain = length(li$ids),
                     t_train = meta$timepoint[meta$uid == utils::tail(li$ids, 1)],
                     train_treat = paste(treat_lut[paste0("u", li$ids)], collapse = ""), xv = xv,
                     stringsAsFactors = FALSE)
  empty <- data.frame(pred0 = NaN, pred = NaN, base = NaN, metric = NaN, tt = NaN, angle = NaN, denoised = NaN,
                      entropy_pred = NaN, test_treat = NaN, entropy_base = NaN, entropy_true = NaN)
  paths <- list.files(file.path(fit_dir, "abm_preds"), full.names = TRUE)
  if (!ok || !length(paths)) return(cbind(head, empty))
  preds <- lapply(paths, readRDS)
  inputs <- list.files(inputs_dir)
  xp <- lapply(li$dec1, salehi_prediction_scores, preds = preds, meta = meta, inputs = inputs,
               inputs_dir = inputs_dir, treat_lut = treat_lut)
  xp <- do.call(rbind, xp[!vapply(xp, is.null, logical(1))])
  if (is.null(xp)) cbind(head, empty) else cbind(head, xp)
}
