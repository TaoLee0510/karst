# Obsolete (2026-11-18): the "auto" arm failed its criterion in round
# 05a8855 and karst_auto is no method of the current code (R/fit.R,
# R/labels.R); kept for the record of the alfak2 era, with compose_auto.R.
#
# karst's "auto" arm: per input with more than two time points, which of the
# round's "all" fits (anchored, pairwise, state space) a rule picks. Nothing
# is refitted; scripts/compose_auto.R links the picked fit as method
# karst_auto.
#
#   Rscript scripts/auto_picks.R [SHARE_CUT [WEIGHT_CUT]]     (default 0.3 0.5)
#
# The rule (multi-timepoint session; thresholds 0.3 / 0.5 chosen on its e1-e3
# evaluation before any E4 scoring, replacing the provisional 0.3 / 0.8):
#   ends_major  the least share, over the time points in between, of a time point's
#               cells on the majors of the first and last time point. Majors:
#               karst:::sparse_joint_domain(counts, 0.01, 3, 5000, min_cn,
#               max_cn, c(1, K)), its first n_major labels; min_cn as the fits
#               (0 for Salehi, the default otherwise), max_cn the largest
#               observed copy number + 1 (the fit's default).
#   weight      the anchored "all" fit's diagnostics$intermediate_weight
#               (run.Rds); NA when it has none (no time point left in between,
#               or no fit).
#   pick        ends_major < SHARE_CUT: pairwise when weight is NA or
#               >= WEIGHT_CUT, otherwise state space; else anchored.
# An input left with two time points by the pair fallback counts as the picked
# method (its fit is a two-time-point fit).
#
# Writes BENCH_RESULTS/tables/auto_picks.tsv: source, dataset, design,
# time points, ends_major, weight, pick (a method name), share_cut, weight_cut.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
share_cut <- if (length(args) >= 1L) as.numeric(args[[1]]) else 0.3
weight_cut <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.5
stopifnot(is.finite(share_cut), is.finite(weight_cut))

auto_pick <- function(ends_major, weight) {
  ifelse(ends_major < share_cut, ifelse(is.na(weight) | weight >= weight_cut, "karst_pairwise", "karst_state_space"),
         "karst")
}

# Least share of a time point in between on the first and last time point's majors.
ends_major_share <- function(yi, pm, min_cn = NULL) {
  x <- as.matrix(yi$x)
  storage.mode(x) <- "double"
  counts <- karst_counts(x, input_times(yi), pm)
  k <- ncol(counts)
  if (k <= 2L) return(NA_real_)
  if (is.null(min_cn)) min_cn <- eval(formals(karst::fit_sparse_joint)$min_cn)
  max_cn <- max(karst:::parse_karyotypes(rownames(counts))) + 1L
  domain <- karst:::sparse_joint_domain(counts, 0.01, 3, 5000, as.integer(min_cn), as.integer(max_cn), c(1L, k))
  majors <- domain$labels[seq_len(domain$n_major)]
  on <- rownames(counts) %in% majors
  min(colSums(counts[on, 2:(k - 1L), drop = FALSE]) / colSums(counts[, 2:(k - 1L), drop = FALSE]))
}

anchored_weight <- function(source, dataset, design) {
  dirs <- Sys.glob(file.path(paths$results, "fits", source, dataset, design, "karst", "minobs_*"))
  for (d in dirs) {
    path <- file.path(d, "run.Rds")
    if (!file.exists(path)) next
    w <- readRDS(path)$intermediate_weight
    if (!is.null(w) && length(w)) return(w[[1L]])
  }
  NA_real_
}

rows <- list()
for (src in c("abm", "misseg")) {
  raw <- if (src == "abm") paths$abm_raw else paths$misseg_raw
  sims <- sub("\\.Rds$", "", list.files(raw, pattern = "\\.Rds$"))
  designs <- abm_designs[!abm_designs$ends, ]
  for (id in sims) {
    sim <- read_abm_sim(id, raw)
    pm <- if (src == "abm") bench_settings$abm$pm else abm_sim_info(id)$pmis
    for (j in seq_len(nrow(designs))) {
      yi <- abm_input(sim, designs[j, ])
      rows[[length(rows) + 1L]] <- data.frame(source = src, dataset = id, design = designs$design[j],
                                              samples = ncol(yi$x), ends_major = ends_major_share(yi, pm),
                                              weight = anchored_weight(src, id, designs$design[j]))
    }
  }
}
traj <- salehi_trajectories(paths)
s <- bench_settings$salehi
for (id in traj$traj_id[traj$length > 2]) {
  yi <- salehi_input(id, ends = FALSE, paths)
  rows[[length(rows) + 1L]] <- data.frame(source = "salehi", dataset = id, design = "all", samples = ncol(yi$x),
                                          ends_major = ends_major_share(yi, s$pm, s$karst_args$min_cn),
                                          weight = anchored_weight("salehi", id, "all"))
}
picks <- do.call(rbind, rows)
picks$pick <- auto_pick(picks$ends_major, picks$weight)
picks$share_cut <- share_cut
picks$weight_cut <- weight_cut
out <- file.path(paths$results, "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.table(picks, file.path(out, "auto_picks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("auto picks (share_cut", share_cut, " weight_cut", weight_cut, "):\n")
print(table(source = picks$source, pick = picks$pick))
cat("weight NA:", sum(is.na(picks$weight)), " ends_major NA:", sum(is.na(picks$ends_major)), "\n")
