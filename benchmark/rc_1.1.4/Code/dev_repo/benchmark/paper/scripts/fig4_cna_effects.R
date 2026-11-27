# The paper's Fig. 4 (06_cna_fitness_effects.R) on the Salehi fits of one arm
# (R/paper_fig456.R): bootstrap draws of non-overlapping lineages, the
# fitness effects (delta f) of every one-chromosome change of each frequent
# karyotype, and the paper's models:
#   4d  glmmTMB Gamma(log) of |delta f| (scaled) ~ context (PDX vs in vitro)
#       and ~ treatment (cisplatin vs control, PDX only), random effects
#       pdx/fi/k_id; Gaussian models of the scaled log variance of delta f
#       per frequent karyotype ~ context and ~ treatment
#   4f  lmer of the similarity (atanh of the correlation) of two frequent
#       karyotypes' delta f profiles ~ log1p(distance) * pair type +
#       treatment pair + (1|fi1) + (1|fi2)
# The variance of delta f per frequent karyotype is taken without na.rm, as
# 06 does (line 336): a frequent karyotype with a neighbour missing from the
# landscape (delta f NA) gets no variance and drops out of the variance
# models. Each draw records how many drop out (kvar_frequent, kvar_na) and
# the two variance models with na.rm = TRUE (PDX_kv_narm, cisplatin_kv_narm),
# as a sensitivity check beside the paper's.
# Fits as the paper selects them: xv > 0, per trajectory the minobs with the
# highest xv (of a tie the lowest minobs), lineages of at least 3 passages.
#
#   Rscript scripts/fig4_cna_effects.R --index N        (array task N)
#   Rscript scripts/fig4_cna_effects.R ARM_KEY FROM TO   (bootstrap draws FROM..TO)
#
# Task N runs arm ((N - 1) %/% CHUNKS) + 1 of fig456_arms, draws of chunk
# ((N - 1) %% CHUNKS) + 1 (BOOT = 200 draws in CHUNKS = 10 chunks of 20).
# Draw b uses seed 1000 + b in every arm; the first chunk also keeps one
# example draw (seed 42, as the paper's Fig. 4b-c) with its delta f table.
# FIG456_EXCLUDE=0 keeps the trajectories that are not analysable (the paper
# used all). Writes BENCH_RESULTS/analysis/fig456/fig4/<arm>_<from>-<to>.Rds.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages({
  library(glmmTMB)
  library(lme4)
  library(dplyr)
})
paths <- bench_paths()
BOOT <- 200L
CHUNKS <- 10L
args <- commandArgs(TRUE)
if (identical(args[1], "--index")) {
  n <- as.integer(args[2])
  arm_i <- (n - 1L) %/% CHUNKS + 1L
  size <- BOOT %/% CHUNKS
  from <- ((n - 1L) %% CHUNKS) * size + 1L
  to <- from + size - 1L
} else {
  arm_i <- match(args[1], fig456_arms$key)
  from <- as.integer(args[2])
  to <- as.integer(args[3])
}
stopifnot(!is.na(arm_i), arm_i <= nrow(fig456_arms))
arm <- fig456_arms[arm_i, ]
exclude <- Sys.getenv("FIG456_EXCLUDE", "1") != "0"
out <- file.path(paths$results, "analysis", "fig456", if (exclude) "fig4" else "fig4_all_inputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cat("Fig. 4, arm", arm$key, " draws", from, "-", to, " exclusion", exclude, "\n")

## Fits (load_fits of 06) ------------------------------------------------------------
rows <- fig456_fit_rows(arm, paths, exclude)
if (arm$method == "alfak_paper") rows <- rows[vapply(rows, ncol, integer(1)) == 13L]
x0 <- do.call(rbind, lapply(rows, function(r) utils::head(r[c("fi", "min_obs", "ntrain", "train_treat", "xv")], 1)))
x0 <- x0[!is.na(x0$xv) & x0$xv > 0, ]
# One fit per trajectory: karst derives every minobs from one fit, so xv is
# often tied between minobs (to rounding), and a trajectory kept twice breaks
# the profile correlations (06 keeps every tie; the paper's fits have none).
# The lowest minobs of a tie is kept (R/select.R).
x0 <- select_by_cv(x0, "fi", "min_obs", "xv")
x0 <- x0[!duplicated(x0$fi), ]
rownames(x0) <- NULL
lineages_raw <- readRDS(paths$salehi_lineages)
l <- lapply(lineages_raw, function(li) paste0(utils::head(li$ids, -1), "-", utils::tail(li$ids, -1)))
meta <- utils::read.csv(paths$salehi_meta, stringsAsFactors = FALSE)
cat("fits with xv > 0 (trajectories):", length(unique(x0$fi)), "\n")

## The paper's functions (06), the landscape path replaced ----------------------------
get_reps <- function(lineages, df) {
  lineages <- lineages[names(lineages) %in% df$fi & lengths(lineages) > 1]
  if (!length(lineages)) return(data.frame(fit_id = character(), n_pass = integer()))
  sizes <- lengths(lineages)
  covered <- chosen <- character()
  repeat {
    ok <- names(lineages)[!vapply(lineages, function(p) any(p %in% covered), logical(1))]
    if (!length(ok)) break
    fit <- sample(ok, 1, prob = sizes[ok])
    chosen <- c(chosen, fit)
    covered <- c(covered, lineages[[fit]])
  }
  data.frame(fit_id = chosen, n_pass = sizes[chosen], stringsAsFactors = FALSE)[sizes[chosen] > 1, ]
}
annotate_samples <- function(df) {
  key <- paste(meta$datasetname, meta$timepoint, sep = "_")
  lut_pdx <- stats::setNames(meta$PDX_id, key)
  lut_tp <- stats::setNames(meta$timepoint, key)
  roots <- sapply(strsplit(df$fi, "_"), function(v) {
    i <- which(v == "l")
    paste(v[1:(i - 1)], collapse = "_")
  })
  df$pdx <- lut_pdx[roots]
  df$tp <- lut_tp[roots]
  df
}
build_deltaf_df <- function(df) {
  lst <- lapply(seq_len(nrow(df)), function(i) {
    ls <- fig456_landscape(fig456_fit_dir(arm, df$fi[i], df$min_obs[i], paths))
    if (is.null(ls)) return(NULL)
    lut <- stats::setNames(ls$mean, ls$k)
    focal <- ls$k[ls$fq]
    dfi <- do.call(rbind, lapply(focal, function(fk) {
      pl <- mean(as.numeric(strsplit(fk, "\\.")[[1]]))
      nn <- apply(paper_gen_all_neighbours(fk, remove_nullisomes = FALSE), 1, paste, collapse = ".")
      data.frame(ploidy = pl, f = lut[fk], deltaf = lut[fk] - lut[nn],
                 mut_id = rep(paste0(rep(1:22, each = 2), c("-", "+")), 1), k_id = fk, stringsAsFactors = FALSE)
    }))
    trt <- strsplit(df$train_treat[i], "")[[1]]
    dfi$treat_switch <- length(unique(trt)) > 1
    dfi$treated <- grepl("y", df$train_treat[i], fixed = TRUE)
    dfi$ntrain <- df$ntrain[i]
    dfi$pdx <- df$pdx[i]
    dfi$fi <- df$fi[i]
    dfi
  })
  do.call(rbind, lst)
}
derive_cat_vars <- function(df) {
  df$wgd <- df$ploidy > 2.75
  df$context <- factor(ifelse(df$pdx %in% c("SA039", "SA906"), "invitro", "PDX"))
  df$treat <- factor(ifelse(df$treated, "cisplatin", "control"))
  df$feat <- interaction(df$mut_id)
  df$context <- stats::relevel(df$context, ref = "invitro")
  df$treat <- stats::relevel(df$treat, ref = "control")
  df
}
compute_pairwise <- function(df) {
  ylist <- split(df, interaction(df$k_id, df$fi))
  ylist <- ylist[sapply(ylist, nrow) > 0]
  kdat <- do.call(rbind, lapply(ylist, function(z) z[1, ]))
  cmb <- utils::combn(nrow(kdat), 2)
  pair <- data.frame(dk = NA, sim = NA, fi1 = kdat$fi[cmb[1, ]], fi2 = kdat$fi[cmb[2, ]],
                     pdx1 = kdat$pdx[cmb[1, ]], pdx2 = kdat$pdx[cmb[2, ]],
                     treat1 = kdat$treat[cmb[1, ]], treat2 = kdat$treat[cmb[2, ]],
                     k1 = kdat$k_id[cmb[1, ]], k2 = kdat$k_id[cmb[2, ]], stringsAsFactors = FALSE)
  for (j in seq_len(ncol(cmb))) {
    i1 <- cmb[1, j]
    i2 <- cmb[2, j]
    v1 <- as.integer(strsplit(pair$k1[j], "\\.")[[1]])
    v2 <- as.integer(strsplit(pair$k2[j], "\\.")[[1]])
    pair$dk[j] <- sum(abs(v1 - v2))
    pair$sim[j] <- stats::cor(ylist[[i1]]$deltaf, ylist[[i2]]$deltaf, use = "pairwise.complete.obs")
  }
  pair$pair_type <- with(pair, ifelse(fi1 == fi2, "same_traj", ifelse(pdx1 == pdx2, "parallel", "diff_line")))
  pair$treat_pair <- with(pair, ifelse(treat1 == treat2, treat1, "mixed"))
  pair$log_dk <- log1p(pair$dk)
  pair
}
deltaf_table <- function() {
  rep <- get_reps(l, x0)
  if (!nrow(rep)) return(NULL)
  x <- annotate_samples(x0[x0$fi %in% rep$fit_id, ])
  d <- build_deltaf_df(x)
  if (is.null(d) || !nrow(d)) NULL else derive_cat_vars(d)
}
# One draw (one_boot of 06); a model that cannot be fitted gives NA.
coef_or_na <- function(fit, term) {
  if (inherits(fit, "error")) return(NA_real_)
  cf <- summary(fit)$coefficients$cond
  if (term %in% rownames(cf)) cf[term, "Estimate"] else NA_real_
}
one_boot <- function(b) {
  set.seed(1000L + b)
  df <- deltaf_table()
  if (is.null(df)) {
    return(list(df1 = data.frame(rep = b, PDX_v = NA_real_, cisplatin_v = NA_real_, PDX_kv = NA_real_,
                                 cisplatin_kv = NA_real_, fits = 0L, frequent = 0L, pairs = 0L, kvar_frequent = 0L,
                                 kvar_na = 0L, PDX_kv_narm = NA_real_, cisplatin_kv_narm = NA_real_), df2 = NULL))
  }
  df$absdf <- abs(df$deltaf) + 1e-6
  df$absdf <- df$absdf / stats::sd(df$absdf, na.rm = TRUE)
  try_fit <- function(expr) tryCatch(expr, error = function(e) e)
  fit_var_context <- try_fit(glmmTMB(absdf ~ context + (1 | pdx / fi / k_id), family = Gamma("log"), data = df))
  fit_var_treat <- try_fit(glmmTMB(absdf ~ treat + (1 | pdx / fi / k_id), family = Gamma("log"),
                                   data = subset(df, context == "PDX")))
  # 06 line 336 as it is (no na.rm; see the header)
  kvar <- df %>% group_by(pdx, fi, k_id, context, treat) %>%
    summarise(v = stats::var(deltaf), .groups = "drop") %>%
    mutate(logv = scale(log(v + 1e-6)))
  fit_kvar_context <- try_fit(glmmTMB(logv ~ context + (1 | pdx / fi), family = gaussian(), data = kvar))
  fit_kvar_treat <- try_fit(glmmTMB(logv ~ treat + (1 | pdx / fi), family = gaussian(),
                                    data = subset(kvar, context == "PDX")))
  # sensitivity: the variance over the neighbours present
  kvar_narm <- df %>% group_by(pdx, fi, k_id, context, treat) %>%
    summarise(v = stats::var(deltaf, na.rm = TRUE), .groups = "drop") %>%
    mutate(logv = scale(log(v + 1e-6)))
  fit_kvar_context_narm <- try_fit(glmmTMB(logv ~ context + (1 | pdx / fi), family = gaussian(), data = kvar_narm))
  fit_kvar_treat_narm <- try_fit(glmmTMB(logv ~ treat + (1 | pdx / fi), family = gaussian(),
                                         data = subset(kvar_narm, context == "PDX")))
  pair <- compute_pairwise(df)
  pair$treat_pair[pair$treat_pair %in% c(1, 2)] <- "same"
  pair$z <- atanh(pair$sim)
  cross_fit <- try_fit(lmer(z ~ log1p(dk) * pair_type + treat_pair + (1 | fi1) + (1 | fi2), data = pair))
  df2 <- if (inherits(cross_fit, "error")) NULL else data.frame(stats::coef(summary(cross_fit)))
  if (!is.null(df2)) {
    df2$param <- rownames(df2)
    df2$rep <- b
  }
  df1 <- data.frame(rep = b, PDX_v = coef_or_na(fit_var_context, "contextPDX"),
                    cisplatin_v = coef_or_na(fit_var_treat, "treatcisplatin"),
                    PDX_kv = coef_or_na(fit_kvar_context, "contextPDX"),
                    cisplatin_kv = coef_or_na(fit_kvar_treat, "treatcisplatin"),
                    fits = length(unique(df$fi)), frequent = length(unique(paste(df$fi, df$k_id))),
                    pairs = nrow(pair), kvar_frequent = nrow(kvar), kvar_na = sum(is.na(kvar$v)),
                    PDX_kv_narm = coef_or_na(fit_kvar_context_narm, "contextPDX"),
                    cisplatin_kv_narm = coef_or_na(fit_kvar_treat_narm, "treatcisplatin"))
  list(df1 = df1, df2 = df2)
}

started <- proc.time()[["elapsed"]]
res <- lapply(from:to, function(b) {
  r <- one_boot(b)
  cat(sprintf("draw %d: %d fits, %d frequent karyotypes, %.0f s\n", b, r$df1$fits, r$df1$frequent,
              proc.time()[["elapsed"]] - started))
  r
})
example <- NULL
if (from == 1L) {
  set.seed(42)
  example <- deltaf_table()
}
saveRDS(list(arm = arm, exclude = exclude, df1 = do.call(rbind, lapply(res, `[[`, "df1")),
             df2 = do.call(rbind, lapply(res, `[[`, "df2")), example = example, x0 = x0),
        file.path(out, sprintf("%s_%03d-%03d.Rds", arm$key, from, to)))
cat("Done.\n")
