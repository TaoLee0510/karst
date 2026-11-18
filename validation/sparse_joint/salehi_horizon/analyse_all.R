#!/usr/bin/env Rscript
# Salehi forecasts against later passages, for fits of every passage of a
# trajectory ("all" design, trajectories of three or more passages): karst
# anchored (the default), pairwise and state_space, alfak, and karst "free",
# and beside them the two-time-point ("ends") fits of the same trajectories.
# Score and intervals as in analyse.R: log-likelihood per cell of the later
# passage under the forecast minus under persistence (mixture weight 0.01),
# paired on the rows (trajectory, minobs, later passage) the compared fits
# share; bootstrap over lineage types and over samples. Horizon in days, as
# in analyse.R.
#
# The two designs choose their landscapes' karyotypes from different time
# points: the "all" landscape is larger, and each fit's own scores
# (ll_fc_e2, ll_p0_e2) mix the uniform distribution over its own landscape.
# Between designs (section 2) only scores on one common support are
# comparable: forecast_score.R's ll_*_s_e2, on the union of the landscapes of
# every fit of the trajectory and minobs, with one uniform floor and one
# persistence baseline. Rows written before those columns existed have none;
# section 2 is then left out (it needs forecast_score.R to be run again).
#   Rscript validation/sparse_joint/salehi_horizon/analyse_all.R [RESULTS]
# RESULTS: the folder with rows/ and rows_all/ (default results/ next to this
# script); the report goes to RESULTS/report_all.
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
res <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else file.path(script_dir, "results")
out_dir <- file.path(res, "report_all")
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
read_rows <- function(dir, design) {
  parts <- lapply(list.files(dir, pattern = "[.]Rds$", full.names = TRUE), readRDS)
  status <- data.frame(design = design, method = vapply(parts, function(p) p$task$method, character(1)),
                       status = vapply(parts, `[[`, character(1), "status"))
  rows <- do.call(rbind, lapply(parts, function(p) {
    r <- p$rows
    if (is.null(r)) return(NULL)
    r$design <- design
    # scores on the common support (forecast_score.R from 2026-11-18 on)
    for (v in c("ll_fc_s_e2", "ll_p0_s_e2", "support_size")) if (is.null(r[[v]])) r[[v]] <- rep(NA_real_, nrow(r))
    r[c("id", "design", "method", "minobs", "uid", "depth", "days", "treat_parent", "treat_desc", "cells_desc",
        "coverage", "ll_fc_e2", "ll_p0_e2", "ll_fc_e3", "ll_p0_e3", "w_fc", "w_p0", "noise", "ll_fc_s_e2",
        "ll_p0_s_e2", "support_size")]
  }))
  list(rows = rows, status = status)
}
a <- read_rows(file.path(res, "rows_all"), "all")
e <- read_rows(file.path(res, "rows"), "ends")
st <- rbind(a$status, e$status)
cat("Tasks:\n"); print(stats::aggregate(n ~ design + method + status, data = transform(st, n = 1), FUN = sum), row.names = FALSE)
x <- rbind(a$rows, e$rows)
label <- c(karst = "anchored", karst_pairwise = "pairwise", karst_state_space = "state_space", alfak = "alfak",
           karst_free = "free")
x$arm <- paste(x$design, label[x$method], sep = ":")
lineage_map <- c(SA1035T = "SA1035", SA1035U = "SA1035", SA532 = "SA532", SA609 = "SA609", SA609UnBU = "SA609",
                 SA906a = "p53 k.o", SA906b = "p53 k.o", SA039U = "p53 w.t", SA535_CISPLATIN_CombinedT = "SA535",
                 SA535_CISPLATIN_CombinedU = "SA535", SA000 = "SA609", SA609R2 = "SA609",
                 SA535_CISPLATIN_CombinedH = "SA535", SA001 = "SA609", SA609H2 = "SA609", SA609H1 = "SA609",
                 SA609R2b = "SA609", SA609TIIDT = "SA609")
x$sample <- sub("_X[0-9]+_l_.*$", "", x$id)
x$type <- unname(lineage_map[x$sample])
x$gain <- x$ll_fc_e2 - x$ll_p0_e2
x$gain3 <- x$ll_fc_e3 - x$ll_p0_e3
x$gain_common <- x$ll_fc_s_e2 - x$ll_p0_s_e2
x$closer <- x$w_fc < x$w_p0
# days ahead, not metadata levels ("depth"): see analyse.R
horizons <- c("<=15", "16-30", "31-60", ">60")
x$horizon <- as.character(cut(x$days, c(0, 15, 30, 60, Inf), labels = horizons))
x$key <- paste(x$id, x$minobs, x$uid)
x <- x[as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", x$id)) >= 3L, ]
cat("\nRows by arm and days ahead:\n"); print(table(x$arm, factor(x$horizon, horizons)))

boot <- function(d, statistic, cluster, draws = 4000L) {
  g <- split(seq_len(nrow(d)), d[[cluster]])
  set.seed(1)
  v <- vapply(seq_len(draws), function(i) statistic(d[unlist(g[sample.int(length(g), replace = TRUE)],
                                                        use.names = FALSE), , drop = FALSE]), numeric(1))
  stats::quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
# score: "gain" (on each fit's own landscape) or "gain_common" (on the common
# support).
wide <- function(arms, score = "gain") {
  keys <- Reduce(intersect, lapply(arms, function(m) x$key[x$arm == m & is.finite(x[[score]])]))
  w <- x[x$arm == arms[[1L]] & x$key %in% keys, c("key", "id", "minobs", "uid", "horizon", "days", "type", "sample")]
  for (m in arms) {
    r <- x[x$arm == m & x$key %in% keys, ]
    r <- r[match(w$key, r$key), ]
    w[[paste0("gain|", m)]] <- r[[score]]
    w[[paste0("ll|", m)]] <- r$ll_fc_e2
    w[[paste0("closer|", m)]] <- r$closer
  }
  w
}
summarise <- function(w, arms) {
  do.call(rbind, lapply(c("all", intersect(horizons, w$horizon)), function(h) {
    d <- if (h == "all") w else w[w$horizon == h, ]
    do.call(rbind, lapply(arms, function(m) {
      ci <- boot(d, function(b) mean(b[[paste0("gain|", m)]]), "type")
      data.frame(horizon = h, arm = m, rows = nrow(d), trajectories = length(unique(d$id)),
                 gain = mean(d[[paste0("gain|", m)]]), low = ci[1L], high = ci[2L],
                 closer_than_parent = mean(d[[paste0("closer|", m)]]), stringsAsFactors = FALSE)
    }))
  }))
}
paired <- function(w, pairs) {
  do.call(rbind, lapply(c("all", intersect(horizons, w$horizon)), function(h) {
    d <- if (h == "all") w else w[w$horizon == h, ]
    do.call(rbind, lapply(pairs, function(p) {
      # Gains over persistence. Within one design they differ as the
      # log-likelihoods do (one landscape, one persistence baseline); between
      # designs only on the common support (section 2).
      stat <- function(z) mean(z[[paste0("gain|", p[1])]] - z[[paste0("gain|", p[2])]])
      lt <- boot(d, stat, "type"); sa <- boot(d, stat, "sample")
      data.frame(horizon = h, comparison = paste(p[1], "-", p[2]), rows = nrow(d), difference = stat(d),
                 low_lineage = lt[1L], high_lineage = lt[2L], low_sample = sa[1L], high_sample = sa[2L],
                 share_first_better = mean(d[[paste0("gain|", p[1])]] > d[[paste0("gain|", p[2])]]),
                 stringsAsFactors = FALSE)
    }))
  }))
}
show <- function(t, title) {
  cat("\n==", title, "==\n")
  num <- vapply(t, is.numeric, logical(1))
  t[num] <- lapply(t[num], function(v) signif(v, 3))
  print(t, row.names = FALSE)
}

# 1. The fits of every passage: karst's three methods against alfak.
arms <- c("all:anchored", "all:pairwise", "all:state_space", "all:alfak")
w <- wide(arms)
s1 <- summarise(w, arms)
p1 <- paired(w, list(c("all:anchored", "all:alfak"), c("all:pairwise", "all:alfak"),
                     c("all:state_space", "all:alfak"), c("all:state_space", "all:anchored"),
                     c("all:pairwise", "all:anchored")))
show(s1, "All time points: log score over persistence per cell")
show(p1, "All time points: paired differences")
utils::write.csv(s1, file.path(out_dir, "gain_all.csv"), row.names = FALSE)
utils::write.csv(p1, file.path(out_dir, "paired_all.csv"), row.names = FALSE)

# 2. Does fitting every passage help? Each method against its two-time-point
# fit of the same trajectory, and the best of each design; on the common
# support only.
arms2 <- c("all:anchored", "all:state_space", "all:alfak", "ends:anchored", "ends:alfak")
w2 <- wide(arms2, "gain_common")
s2 <- NULL
if (nrow(w2)) {
  s2 <- summarise(w2, arms2)
  p2 <- paired(w2, list(c("all:anchored", "ends:anchored"), c("all:state_space", "ends:anchored"),
                        c("all:alfak", "ends:alfak"), c("ends:anchored", "ends:alfak")))
  show(s2, "All time points against two time points (same trajectories, common support)")
  show(p2, "Paired: every passage against two time points (common support)")
  utils::write.csv(s2, file.path(out_dir, "gain_all_vs_ends.csv"), row.names = FALSE)
  utils::write.csv(p2, file.path(out_dir, "paired_all_vs_ends.csv"), row.names = FALSE)
} else {
  cat("\nSection 2 left out: no rows with scores on the common support (forecast_score.R before",
      "2026-11-18); the gains of the two designs are on different landscapes and are not compared.\n")
}

# 3. karst "free" (every passage), where its fits exist.
arms3 <- c("all:anchored", "all:free", "all:alfak")
w3 <- wide(arms3)
p3 <- paired(w3, list(c("all:free", "all:anchored"), c("all:free", "all:alfak")))
show(p3, "karst free, every passage")
utils::write.csv(p3, file.path(out_dir, "paired_all_free.csv"), row.names = FALSE)
saveRDS(list(all = w, all_vs_ends = w2, free = w3), file.path(out_dir, "rows.Rds"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  theme_set(theme_bw(base_size = 11))
  # every passage against two (common support) where available, else the
  # fits of every passage alone
  both <- !is.null(s2)
  s <- if (both) s2[s2$horizon != "all", ] else s1[s1$horizon != "all", ]
  s$horizon <- factor(s$horizon, levels = horizons)
  s$arm <- factor(s$arm, levels = if (both) arms2 else arms)
  f1 <- ggplot(s, aes(horizon, gain, colour = arm, group = arm)) + geom_hline(yintercept = 0, colour = "grey50") +
    geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    labs(x = "days ahead of the parent", y = "log score over persistence (per cell)", colour = NULL,
         title = if (both) "Trajectories of 3+ passages: every passage (all) and two (ends), common support" else
           "Trajectories of 3+ passages: fits of every passage")
  ggsave(file.path(out_dir, "figures", "horizon_all.png"), f1, width = 8, height = 4.5, dpi = 150)
  ggsave(file.path(out_dir, "figures", "horizon_all.pdf"), f1, width = 8, height = 4.5)
}
