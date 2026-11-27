#!/usr/bin/env Rscript
# D4 of dev/realdata/PLAN.md (C-E, clarification C1): three interval forms
# computed after the fit, chosen on development data, judged once on held-out.
#   Rscript dev/realdata/CE/d4.R select   (development only: dev GRF dev split, Salehi dev)
#   Rscript dev/realdata/CE/d4.R judge    (held-out: GRF heldout, RMF, NK, e1-e3, Salehi confirm, GBM, Ovarian)
# Reads cd342db base fits: DEV_RESULTS/phase2/fits/ce_base/<kind>, phase1/eset/ce_base_x1/<set>,
# realdata/B (Salehi default, NUMBAT rate2 and halves). Writes realdata/CE/d4_<mode>.rds.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
source(file.path(dev_repo, "dev", "R", "calibration.R"))
options(width = 220)
mode <- commandArgs(TRUE)[1]
# Development coverage threshold: 0.95 as pre-registered; 0.945 by the user's decision 1 (PLAN.md, written 2026-11-20 01:15).
threshold <- as.numeric(Sys.getenv("D4_THRESHOLD", "0.945"))
stopifnot(mode %in% c("select", "judge"))
out_dir <- file.path(rd_results, "CE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
z975 <- stats::qnorm(0.975)
floor0 <- 0.009
bins <- c(0, 1, 2, 4, 8, Inf)

## ---------- nodes with truth (simulations) ----------
sim_nodes <- function(kind) {
  files <- if (kind %in% c("e1", "e2", "e3")) {
    list.files(file.path(dev_results, "phase1", "eset", "ce_base_x1", kind), full.names = TRUE)
  } else list.files(file.path(dev_results, "phase2", "fits", "ce_base", kind), full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f)
    if (!identical(r$status, "ok")) return(NULL)
    if (kind %in% c("e1", "e2", "e3")) {
      n <- r$nodes
      s <- r$summary
      at <- match(n$karyotype, s$karyotype)
      n$count_t0 <- ifelse(is.na(at), 0, s$count_t0[at]); n$count_t1 <- ifelse(is.na(at), 0, s$count_t1[at])
    } else n <- r$near_nodes
    n$dataset_id <- r$dataset_id
    n$kind <- kind
    n[c("dataset_id", "kind", "karyotype", "distance", "count_t0", "count_t1", "truth_relative", "pred_relative",
        "available", "model_sd", "shell", "relative")]
  }))
}
prep <- function(x, birth_rate = 1) {
  x <- x[x$available & is.finite(x$model_sd) & is.finite(x$truth_relative) & is.finite(x$pred_relative) &
           x$distance <= 2, ]
  x$est <- x$pred_relative
  x$err <- x$truth_relative - x$est       # truth - estimate
  x$u <- sqrt(x$model_sd^2 + (floor0 * birth_rate * (1 + x$shell))^2)
  x$class <- node_class(x$distance, x$count_t0 + x$count_t1)
  x$z <- x$est / x$model_sd
  x
}
node_class <- function(distance, total) {
  ifelse(distance == 0, ifelse(total >= 5, "obs5+", "obs<5"), paste0("d", distance))
}
# a model sd of 0 (the reference) puts a nonzero estimate in the top bin, a zero one in the first
zbin <- function(z) { z <- abs(z); z[is.nan(z)] <- 0; as.character(cut(z, bins, include.lowest = TRUE)) }

## ---------- weighted quantiles (each dataset the same weight) ----------
wq <- function(v, id, p) {
  w <- 1 / as.numeric(table(id)[id])
  o <- order(v); v <- v[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  v[which(cw >= p)[1]]
}

## ---------- the three forms ----------
fit_F1 <- function(x) {
  nodes <- data.frame(dataset_id = x$dataset_id, distance = x$distance, error = -x$err, model_sd = x$model_sd,
                      relative = x$est, shell = x$shell, count_t0 = x$count_t0, count_t1 = x$count_t1)
  list(form = "F1", constants = dev_calib_fit(nodes))
}
# bins with fewer than 200 development nodes merged into the next lower bin
bin_levels <- levels(cut(0, bins, include.lowest = TRUE))
merge_bins <- function(x) {
  map <- list()
  for (c in unique(x$class)) {
    counts <- table(factor(zbin(x$z[x$class == c]), levels = bin_levels))
    target <- bin_levels
    for (i in seq_along(bin_levels)) if (counts[[i]] < 200 && i > 1) target[i] <- target[i - 1]
    map[[c]] <- stats::setNames(target, bin_levels)
  }
  map
}
map_bins <- function(map, class, z) {
  zb <- zbin(z)
  out <- character(length(z))
  for (c in unique(class)) {
    i <- class == c
    m <- if (is.null(map[[c]])) stats::setNames(bin_levels, bin_levels) else map[[c]]
    out[i] <- m[zb[i]]
  }
  out
}
fit_F2 <- function(x, p) {
  map <- merge_bins(x)
  x$bin <- map_bins(map, x$class, x$z)
  s <- sign(x$est); s[s == 0] <- 1
  toward <- pmax(0, -x$err * s) / x$u
  away <- pmax(0, x$err * s) / x$u
  key <- paste(x$class, x$bin)
  q <- do.call(rbind, lapply(split(seq_len(nrow(x)), key), function(i) data.frame(
    key = key[i[1]], n = length(i), q_t = wq(toward[i], x$dataset_id[i], p), q_a = wq(away[i], x$dataset_id[i], p))))
  list(form = "F2", p = p, map = map, q = q)
}
fit_F3 <- function(x, p) {
  r <- abs(x$err) / x$u
  q <- vapply(split(seq_len(nrow(x)), x$class), function(i) wq(r[i], x$dataset_id[i], p), 0)
  list(form = "F3", p = p, q = q)
}
# interval [lo, hi] for rows with est, model_sd, u, class, z, relative, shell, distance
interval <- function(form, x, birth_rate = 1) {
  if (form$form == "F1") {
    k <- form$constants
    deficit <- ifelse(x$distance == 0, k[["deficit"]], k[["deficit"]] * 1.1 / 0.6)
    sd <- sqrt(x$model_sd^2 + (deficit * x$est)^2 + (k[["floor"]] * birth_rate * (1 + x$shell))^2)
    return(cbind(lo = x$est - z975 * sd, hi = x$est + z975 * sd))
  }
  if (form$form == "F3") {
    h <- form$q[x$class] * x$u
    return(cbind(lo = x$est - h, hi = x$est + h))
  }
  bin <- map_bins(form$map, x$class, x$z)
  at <- match(paste(x$class, bin), form$q$key)
  s <- sign(x$est); s[s == 0] <- 1
  toward <- form$q$q_t[at] * x$u; away <- form$q$q_a[at] * x$u
  lo <- ifelse(s > 0, x$est - toward, x$est - away)
  hi <- ifelse(s > 0, x$est + away, x$est + toward)
  cbind(lo = lo, hi = hi)
}
current_form <- list(form = "F1", constants = c(floor = 0.009, deficit = 0.6))

## ---------- simulation metrics ----------
sim_metrics <- function(form, x) {
  iv <- interval(form, x)
  x$covered <- x$truth_relative >= iv[, "lo"] & x$truth_relative <= iv[, "hi"]
  x$excl0 <- iv[, "lo"] > 0 | iv[, "hi"] < 0
  per <- stats::aggregate(covered ~ dataset_id + class, x, mean)
  cov <- stats::aggregate(covered ~ class, per, mean)
  ex <- x[x$excl0, ]
  sign_ok <- if (nrow(ex)) mean(sign(ex$truth_relative) == sign(ex$est)) else NA
  list(coverage = stats::setNames(cov$covered, cov$class), excl0 = mean(x$excl0), sign_ok = sign_ok,
       width = stats::median(iv[, "hi"] - iv[, "lo"]))
}

## ---------- real data ----------
B <- file.path(rd_results, "B")
salehi_a <- readRDS(file.path(rd_results, "A", "salehi.rds"))
salehi_a$group <- ifelse(salehi_a$line %in% c("SA609", "SA039"), "Salehi dev", "Salehi confirm")
numbat_sets <- readRDS(file.path(rd_results, "A", "numbat.rds"))$datasets[c("cancer", "sample")]
read_rec <- function(kind, cancer, sample, variant) {
  f <- file.path(B, kind, cancer, sprintf("%s__%s.rds", sample, variant))
  if (file.exists(f)) readRDS(f) else NULL
}
rec_birth_rate <- function(rec) if (!is.null(rec$birth_rate)) rec$birth_rate else 1
real_rows <- function(s, birth_rate) {
  x <- data.frame(karyotype = s$karyotype, est = s$relative_fitness, model_sd = s$relative_sd_model, shell = s$shell,
                  total = s$count_t0 + s$count_t1, stringsAsFactors = FALSE)
  x$distance <- ifelse(x$total > 0, 0L, 1L)
  x$u <- sqrt(x$model_sd^2 + (floor0 * birth_rate * (1 + x$shell))^2)
  x$class <- node_class(x$distance, x$total)
  x$z <- x$est / x$model_sd
  x[is.finite(x$est) & is.finite(x$model_sd) & x$model_sd > 0, ]
}
informativeness <- function(form, rec) {
  x <- real_rows(rec$summary, rec_birth_rate(rec))
  x <- x[x$total >= 10, ]
  if (nrow(x) < 3) return(NA)
  iv <- interval(form, x, rec_birth_rate(rec))
  mean(iv[, "lo"] > 0 | iv[, "hi"] < 0)
}
# half-width of each estimate's interval on the side facing `other`
facing_sd <- function(iv, est, other) ifelse(other >= est, iv[, "hi"] - est, est - iv[, "lo"]) / z975
lineage_share <- function(form, group) {
  ids <- salehi_a$id[salehi_a$group == group]
  recs <- lapply(stats::setNames(ids, ids), function(id) read_rec("salehi", "salehi", id, "default"))
  recs <- recs[vapply(recs, function(r) identical(r$status, "ok"), logical(1))]
  sal <- salehi_a[salehi_a$id %in% names(recs), ]
  keys <- split(sal, list(sal$line, sal$treat_first, sal$treat_last), drop = TRUE)
  out <- do.call(rbind, lapply(keys, function(g) {
    if (nrow(g) < 2) return(NULL)
    cmb <- utils::combn(nrow(g), 2)
    do.call(rbind, lapply(seq_len(ncol(cmb)), function(j) {
      a <- g[cmb[1, j], ]; b <- g[cmb[2, j], ]
      if (length(intersect(c(a$first_uid, a$last_uid), c(b$first_uid, b$last_uid)))) return(NULL)
      x <- real_rows(recs[[a$id]]$summary, rec_birth_rate(recs[[a$id]]))
      y <- real_rows(recs[[b$id]]$summary, rec_birth_rate(recs[[b$id]]))
      x <- x[x$total > 0, ]; y <- y[y$total > 0, ]
      at <- match(x$karyotype, y$karyotype)
      use <- !is.na(at)
      if (sum(use) < 3) return(NULL)
      x <- x[use, ]; y <- y[at[use], ]
      d <- x$est - y$est
      anchor <- x$total >= 10 & y$total >= 10
      shift <- if (any(anchor)) mean(d[anchor]) else mean(d)
      ix <- interval(form, x, rec_birth_rate(recs[[a$id]])); iy <- interval(form, y, rec_birth_rate(recs[[b$id]]))
      # each estimate shifted onto the other's reference before choosing the facing side
      sx <- facing_sd(ix, x$est, y$est + shift); sy <- facing_sd(iy, y$est, x$est - shift)
      z <- (d - shift) / sqrt(sx^2 + sy^2)
      data.frame(line = a$line, n = length(z), share_z = mean(abs(z) > z975))
    }))
  }))
  out
}
halves_within <- function(form, cancer_filter = NULL, kind = "numbat") {
  inputs <- numbat_sets
  if (!is.null(cancer_filter)) inputs <- inputs[inputs$cancer %in% cancer_filter, ]
  do.call(rbind, lapply(seq_len(nrow(inputs)), function(i) {
    h <- lapply(c("half1", "half2"), function(v) read_rec(kind, inputs$cancer[i], inputs$sample[i], v))
    if (any(vapply(h, function(r) !identical(r$status, "ok"), logical(1)))) return(NULL)
    full <- h[[1]]$observed
    a <- real_rows(h[[1]]$summary, rec_birth_rate(h[[1]])); b <- real_rows(h[[2]]$summary, rec_birth_rate(h[[2]]))
    at <- match(a$karyotype, b$karyotype)
    tot <- (full$count_t0 + full$count_t1)[match(a$karyotype, full$karyotype)]
    tot[is.na(tot)] <- 0
    use <- !is.na(at) & tot >= 10
    if (sum(use) < 3) return(NULL)
    a <- a[use, ]; b <- b[at[use], ]
    d <- a$est - b$est
    shift <- mean(d)
    ia <- interval(form, a, rec_birth_rate(h[[1]])); ib <- interval(form, b, rec_birth_rate(h[[2]]))
    sa <- facing_sd(ia, a$est, b$est + shift); sb <- facing_sd(ib, b$est, a$est - shift)
    data.frame(cancer = inputs$cancer[i], sample = inputs$sample[i], n = sum(use),
               within = mean(abs(d - shift) <= z975 * sqrt(sa^2 + sb^2)))
  }))
}
real_informative <- function(form, group = NULL, cancer = NULL) {
  if (!is.null(group)) {
    ids <- salehi_a$id[salehi_a$group == group]
    v <- vapply(ids, function(id) { r <- read_rec("salehi", "salehi", id, "default"); if (identical(r$status, "ok")) informativeness(form, r) else NA_real_ }, 0)
  } else {
    inputs <- numbat_sets; inputs <- inputs[inputs$cancer == cancer, ]
    v <- vapply(inputs$sample, function(s) { r <- read_rec("numbat", cancer, s, "rate2"); if (identical(r$status, "ok")) informativeness(form, r) else NA_real_ }, 0)
  }
  stats::median(v, na.rm = TRUE)
}

## ---------- select ----------
if (mode == "select") {
  dev <- prep(sim_nodes("dev"))
  cat("development nodes:", nrow(dev), "in", length(unique(dev$dataset_id)), "datasets\n")
  candidates <- list(current = current_form, F1 = fit_F1(dev))
  for (p in c(0.975, 0.976, 0.977, 0.978, 0.979, 0.98)) {
    f <- fit_F2(dev, p); m <- sim_metrics(f, dev)
    if (all(m$coverage >= threshold) || p == 0.98) { candidates$F2 <- f; break }
  }
  for (p in c(0.95, 0.952, 0.954, 0.956, 0.958, 0.96)) {
    f <- fit_F3(dev, p); m <- sim_metrics(f, dev)
    if (all(m$coverage >= threshold) || p == 0.96) { candidates$F3 <- f; break }
  }
  rows <- do.call(rbind, lapply(names(candidates), function(n) {
    f <- candidates[[n]]
    m <- sim_metrics(f, dev)
    lin <- lineage_share(f, "Salehi dev")
    data.frame(form = n, p = if (is.null(f$p)) NA else f$p,
               cov_obs_lt5 = m$coverage[["obs<5"]], cov_obs5 = m$coverage[["obs5+"]], cov_d1 = m$coverage[["d1"]],
               cov_d2 = m$coverage[["d2"]], dev_excl0 = m$excl0, dev_sign_ok = m$sign_ok, dev_width = m$width,
               salehi_dev_lineage = mean(lin$share_z), salehi_dev_informative = real_informative(f, "Salehi dev"),
               ovarian_informative = real_informative(f, cancer = "Ovarian"))
  }))
  rows$eligible <- with(rows, cov_obs_lt5 >= threshold & cov_obs5 >= threshold & cov_d1 >= threshold &
                          cov_d2 >= threshold & salehi_dev_lineage <= 0.07) & rows$form != "current"
  cat("development coverage threshold:", threshold, "\n")
  chosen <- if (any(rows$eligible)) rows$form[rows$eligible][which.max(rows$salehi_dev_informative[rows$eligible])] else NA
  print(rows, digits = 3, row.names = FALSE)
  cat("chosen:", chosen, "\n")
  if (!is.null(candidates$F2)) print(candidates$F2$q, digits = 3, row.names = FALSE)
  if (!is.null(candidates$F3)) print(candidates$F3$q, digits = 3)
  print(candidates$F1$constants)
  dev_save_rds(list(rows = rows, candidates = candidates, chosen = chosen, threshold = threshold),
               file.path(out_dir, sprintf("d4_select_%s.rds", threshold)))
  file.copy(file.path(out_dir, sprintf("d4_select_%s.rds", threshold)), file.path(out_dir, "d4_select.rds"), overwrite = TRUE)
}

## ---------- judge ----------
if (mode == "judge") {
  sel <- readRDS(file.path(out_dir, "d4_select.rds"))
  if (is.na(sel$chosen)) stop("No form met the selection threshold; D4 is not judged.")
  forms <- list(chosen = sel$candidates[[sel$chosen]], current = current_form)
  kinds <- c("heldout", "rmf", "nk", "e1", "e2", "e3")
  sims <- lapply(stats::setNames(kinds, kinds), function(k) prep(sim_nodes(k)))
  res <- list()
  for (fn in names(forms)) {
    f <- forms[[fn]]
    sm <- lapply(sims, function(x) sim_metrics(f, x))
    lin <- lineage_share(f, "Salehi confirm")
    hv <- halves_within(f, c("GBM", "Ovarian"))
    res[[fn]] <- list(
      sim = do.call(rbind, lapply(names(sm), function(k) data.frame(form = fn, kind = k, t(sm[[k]]$coverage),
                                                                   excl0 = sm[[k]]$excl0, sign_ok = sm[[k]]$sign_ok,
                                                                   width = sm[[k]]$width, check.names = FALSE))),
      lineage = mean(lin$share_z), lineage_rows = lin,
      halves = hv, halves_med = tapply(hv$within, hv$cancer, stats::median),
      informative = c(salehi_confirm = real_informative(f, "Salehi confirm"), GBM = real_informative(f, cancer = "GBM")))
  }
  ch <- res$chosen; cu <- res$current
  sim <- ch$sim
  obs_ok <- all(sim$`obs<5`[sim$kind %in% c("heldout", "rmf", "e1", "e2", "e3")] >= 0.93 &
                  sim$`obs5+`[sim$kind %in% c("heldout", "rmf", "e1", "e2", "e3")] >= 0.93, na.rm = TRUE)
  un_ok <- all(sim$d1[sim$kind %in% c("heldout", "rmf", "e1", "e2", "e3")] >= 0.90 &
                 sim$d2[sim$kind %in% c("heldout", "rmf", "e1", "e2", "e3")] >= 0.90, na.rm = TRUE)
  nk_ch <- unlist(sim[sim$kind == "nk", c("obs<5", "obs5+", "d1", "d2")])
  nk_cu <- unlist(cu$sim[cu$sim$kind == "nk", c("obs<5", "obs5+", "d1", "d2")])
  verdict <- c(
    D4_1 = obs_ok && un_ok && all(nk_ch >= nk_cu - 0.02, na.rm = TRUE),
    D4_2 = ch$lineage <= 0.10,
    D4_3 = all(ch$halves_med >= 0.90),
    D4_4 = stats::median(ch$informative) >= 0.20 &&
      all(sim$sign_ok[sim$kind %in% c("heldout", "rmf", "e1", "e2", "e3")] >= 0.95, na.rm = TRUE))
  print(rbind(ch$sim, cu$sim), digits = 3, row.names = FALSE)
  cat("lineage (Salehi confirm):", signif(ch$lineage, 3), "current", signif(cu$lineage, 3), "\n")
  cat("halves 10+ (median):\n"); print(rbind(chosen = ch$halves_med, current = cu$halves_med), digits = 3)
  cat("informative:\n"); print(rbind(chosen = ch$informative, current = cu$informative), digits = 3)
  print(verdict)
  cat("VERDICT D4:", if (all(verdict)) "PASSES" else "FAILS", "\n")
  dev_save_rds(list(selected = sel$chosen, results = res, verdict = verdict), file.path(out_dir, "d4_judge.rds"))
}
