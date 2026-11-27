#!/usr/bin/env Rscript
# Step 1 of PLAN.md (rules fixed in "Step 1 rules"): the gate on the
# development data, the choice of k, the pick step and "is closer to the
# truth better?".
#   Rscript validation/sparse_joint/dispersion_gate/step1_analyse.R --v7=DIR --pick=DIR --vm=DIR --vm-x1=DIR \
#     --vm-raw=DIR --cells=FILE --out=DIR
# --v7: step1/v7 (x<m>/<dataset>.rds); --pick: step1/pick; --vm: step1/vm
# (x<m>/fits/misseg/...); --vm-x1: fig2d_misseg/beta_fixed/fits/misseg (the
# x1 KARST fits); --vm-raw: ABM_varyMisseg .Rds; --cells: salehi_cells.tsv.
# Reads results/step0/s0.tsv and results/step1/gate_v7_dev.tsv. Writes to
# --out: v7_fits.tsv, vm_fits.tsv, gate_cells.tsv, k_choice.tsv, closer.tsv,
# pick.tsv, step1.txt (verdict), step1.pdf.
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
# rc freeze: the paper benchmark's R/ as this analysis ran it (the newer
# working tree) is Code/paper; step 0's tables and gate_v7_dev.tsv are read from
# --gate-results (default RC_UPSTREAM/part4_rate_choice/dispersion_gate/results)
paper_dir <- file.path(Sys.getenv("RC_HOME"), "Code", "paper")
source(file.path(here, "common_gate.R"))
out <- arg("out"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_lines <- character()
say <- function(...) { line <- paste0(...); cat(line, "\n"); log_lines <<- c(log_lines, line) }
multiples <- c(1, 2, 4, 8, 16)
gate_results <- arg("gate-results", file.path(Sys.getenv("RC_UPSTREAM"), "part4_rate_choice", "dispersion_gate", "results"))
s0 <- utils::read.delim(file.path(gate_results, "step0", "s0.tsv"))
tau <- stats::setNames(s0$tau[s0$statistic == "once"], s0$n[s0$statistic == "once"])
gate_decide <- function(c0, c1) {
  n <- if (sum(c0) >= 300 && sum(c1) >= 300) 300L else if (sum(c0) >= 100 && sum(c1) >= 100) 100L else NA_integer_
  if (is.na(n)) return(data.frame(n = NA_integer_, once = NA_real_, raised = FALSE))
  once <- gate_stats(c0, c1, n)[["once"]]
  data.frame(n = n, once = once, raised = once >= tau[[as.character(n)]])
}

## v7 ----------------------------------------------------------------------
v7_files <- list.files(arg("v7"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
v7 <- do.call(rbind, lapply(v7_files, function(p) {
  r <- readRDS(p)
  ok <- identical(r$status, "ok")
  met <- if (ok) stats::setNames(r$metrics$rmse, r$metrics$stratum) else c(d0 = NA, d1 = NA, d2 = NA)
  data.frame(dataset_id = r$dataset_id, family = r$family, design = r$design, wavelength = r$wavelength,
             truth = r$truth_multiple, fitted = r$multiple, depth = paste0(r$depth_t0, "/", r$depth_t1), status = r$status,
             rmse_d0 = met[["d0"]], rmse_d1 = met[["d1"]], rmse_d2 = met[["d2"]],
             r_all = if (ok) r$r[["all"]] else NA, r_d1 = if (ok) r$r[["d1"]] else NA,
             objective = if (ok) r$objective else NA, posterior_objective = if (ok) r$posterior_objective else NA,
             stringsAsFactors = FALSE)
}))
gate_v7 <- utils::read.delim(file.path(gate_results, "step1", "gate_v7_dev.tsv"), stringsAsFactors = FALSE)
v7 <- merge(v7, gate_v7[, c("dataset_id", "source", "n", "once", "raised")], by = "dataset_id")
base <- v7[v7$fitted == 1, c("dataset_id", "rmse_d1")]; names(base)[2] <- "rmse_d1_x1"
v7 <- merge(v7, base, by = "dataset_id")
v7$ratio <- v7$rmse_d1 / v7$rmse_d1_x1
utils::write.table(v7, file.path(out, "v7_fits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("v7 fits: ", nrow(v7), " (", sum(v7$status == "ok"), " ok); datasets ", length(unique(v7$dataset_id)))

## ABM_varyMisseg -------------------------------------------------------------
old_wd <- setwd(paper_dir)
Sys.setenv(BENCH_RESULTS = file.path(tempdir(), "bench_unused"))
for (f in list.files("R", full.names = TRUE)) source(f)
setwd(old_wd)
cells <- utils::read.delim(arg("cells"))$cells
sims <- sub("[.]Rds$", "", list.files(arg("vm-raw"), pattern = "^w_.*[.]Rds$"))
designs <- c(w2_ends = 2L, w4_ends = 4L, w8_ends = 8L)
vm <- do.call(rbind, parallel::mclapply(sims, function(id) {
  sim <- readRDS(file.path(arg("vm-raw"), paste0(id, ".Rds")))
  info <- abm_sim_info(id)
  x <- as.matrix(abm_subsample_x(sim$abm_output$x, id, cells))
  times <- as.numeric(colnames(x))
  do.call(rbind, lapply(names(designs), function(dn) {
    tt <- utils::tail(times[times < 120], designs[[dn]])
    c0 <- stats::setNames(as.numeric(x[, times == tt[1]]), rownames(x))
    c1 <- stats::setNames(as.numeric(x[, times == tt[length(tt)]]), rownames(x))
    g <- gate_decide(c0, c1)
    do.call(rbind, lapply(multiples, function(m) do.call(rbind, lapply(c(5, 10, 20), function(mo) {
      root <- if (m == 1) arg("vm-x1") else file.path(arg("vm"), sprintf("x%02d", m), "fits", "misseg")
      dir <- file.path(root, id, dn, "karst", paste0("minobs_", mo))
      row <- data.frame(dataset_id = paste(id, dn, sep = "_"), sim = id, design = dn, minobs = mo,
                        truth = round(info$pmis / 5e-5), fitted = m, n = g$n, once = g$once, raised = g$raised,
                        status = NA_character_, r = NA_real_, rnn = NA_real_, rfq = NA_real_, Rxv = NA_real_,
                        objective = NA_real_, posterior_objective = NA_real_, stringsAsFactors = FALSE)
      if (!file.exists(file.path(dir, "run.Rds"))) { row$status <- "missing"; return(row) }
      run <- readRDS(file.path(dir, "run.Rds"))
      row$status <- run$status
      if (!identical(run$status, "ok")) return(row)
      l <- readRDS(file.path(dir, "landscape.Rds"))
      if (!is.null(l$dist)) { l$fq <- l$dist == 0L; l$nn <- l$dist == 1L }
      sc <- landscape_scores(l, grf_fitness(l$k, info$w, sim$true_landscape), readRDS(file.path(dir, "xval.Rds")))
      row$r <- sc$r; row$rnn <- sc$rnn; row$rfq <- sc$rfq; row$Rxv <- sc$Rxv
      row$objective <- if (is.null(run$objective)) NA_real_ else run$objective
      row$posterior_objective <- if (is.null(run$posterior_objective)) NA_real_ else run$posterior_objective
      row
    }))))
  }))
}, mc.cores = 3))
vm$err <- 1 - vm$r
vb <- vm[vm$fitted == 1, c("dataset_id", "minobs", "err")]; names(vb)[3] <- "err_x1"
vm <- merge(vm, vb, by = c("dataset_id", "minobs"))
vm$ratio <- vm$err / vm$err_x1
utils::write.table(vm, file.path(out, "vm_fits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("varyMisseg rows: ", nrow(vm), " (", sum(vm$status == "ok", na.rm = TRUE), " ok, ", sum(vm$status == "missing"), " missing)")

## gate on the development cells ---------------------------------------------
gv7 <- unique(v7[, c("dataset_id", "source", "wavelength", "truth", "raised")])
gv7$cell <- paste0(gv7$source, " w", gv7$wavelength)
gvm <- unique(vm[, c("dataset_id", "truth", "raised", "n")])
gvm$cell <- "paper_vm (Salehi cells) w0.8"
gcell <- rbind(stats::aggregate(raised ~ cell + truth, gv7, function(v) c(n = length(v), raised = mean(v))),
               stats::aggregate(raised ~ cell + truth, gvm, function(v) c(n = length(v), raised = mean(v))))
gcell <- data.frame(gcell[, c("cell", "truth")], gcell$raised)
utils::write.table(gcell, file.path(out, "gate_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\n== share raised by development cell and truth")
w <- reshape(gcell[, c("cell", "truth", "raised")], idvar = "cell", timevar = "truth", direction = "wide")
print(w, row.names = FALSE, digits = 3)
say("varyMisseg inputs by gate n: ", paste(names(table(gvm$n, useNA = "ifany")), table(gvm$n, useNA = "ifany"), collapse = ", "))
false_raise <- gcell[gcell$truth == 1, ]
fail_false <- any(false_raise$raised > 0.10)
say("max share raised at truth x1 over development cells: ", signif(max(false_raise$raised), 3),
    if (fail_false) "  -> FAIL (> 10%)" else "  (<= 10%)")

## choice of k ----------------------------------------------------------------
med <- function(v) if (length(v)) stats::median(v, na.rm = TRUE) else NA_real_
kc <- do.call(rbind, lapply(c(2, 4, 8), function(k) {
  a <- v7[v7$raised & v7$fitted == k & is.finite(v7$ratio), ]
  b <- vm[vm$raised & vm$fitted == k & is.finite(vm$ratio), ]
  data.frame(k = k, H = med(c(a$ratio[a$truth == 1], b$ratio[b$truth == 1])), H_rows = sum(a$truth == 1) + sum(b$truth == 1),
             G = med(a$ratio[a$truth >= 4]), G_rows = sum(a$truth >= 4),
             G_x4 = med(a$ratio[a$truth == 4]), G_x16 = med(a$ratio[a$truth == 16]),
             G_vm = med(b$ratio[b$truth >= 4]), G_vm_rows = sum(b$truth >= 4))
}))
kc$eligible <- is.finite(kc$H) & kc$H <= 1.25 & is.finite(kc$G) & kc$G <= 0.90
utils::write.table(kc, file.path(out, "k_choice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\n== choice of k (ratio = error at xk / error at x1, raised datasets)")
print(kc, row.names = FALSE, digits = 3)
k <- if (any(kc$eligible)) kc$k[kc$eligible][which.min(kc$G[kc$eligible])] else NA
fail_k <- is.na(k)
fail_vm <- !fail_k && !(kc$G_vm[kc$k == k] < 1)
say("chosen k: ", if (fail_k) "none -> FAIL" else k, if (!fail_k && fail_vm) "; varyMisseg G not below 1 -> FAIL" else "")

## is closer to the truth better? -----------------------------------------------
cl <- rbind(data.frame(set = "v7", stats::aggregate(ratio ~ truth + fitted, v7[is.finite(v7$ratio), ], stats::median)),
            data.frame(set = "varyMisseg", stats::aggregate(ratio ~ truth + fitted, vm[is.finite(vm$ratio), ], stats::median)))
utils::write.table(cl, file.path(out, "closer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\n== median error ratio (fitted xm / x1), all development datasets")
for (s in c("v7", "varyMisseg")) {
  say("-- ", s)
  print(reshape(cl[cl$set == s, c("truth", "fitted", "ratio")], idvar = "truth", timevar = "fitted", direction = "wide"),
        row.names = FALSE, digits = 3)
}

## pick step ------------------------------------------------------------------------
pick_files <- list.files(arg("pick"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
pk <- do.call(rbind, lapply(pick_files, function(p) {
  r <- readRDS(p)
  data.frame(dataset_id = r$dataset_id, fitted = r$multiple, status = r$status,
             nll_half = if (identical(r$status, "ok")) r$nll_half else NA_real_,
             b_inside = if (identical(r$status, "ok")) r$b_inside else NA_real_,
             domain = if (identical(r$status, "ok")) paste(length(r$domain), sum(nchar(r$domain))) else NA_character_,
             stringsAsFactors = FALSE)
}))
pick_rows <- list()
for (kk in c(2, 4, 8)) {
  a <- v7[v7$raised & v7$fitted == 1, c("dataset_id", "truth", "rmse_d1", "objective", "posterior_objective")]
  b <- v7[v7$raised & v7$fitted == kk, c("dataset_id", "rmse_d1", "objective", "posterior_objective")]
  ab <- merge(a, b, by = "dataset_id", suffixes = c("_1", "_k"))
  if (!is.null(pk)) {
    ab <- merge(ab, pk[pk$fitted == 1, c("dataset_id", "nll_half", "b_inside", "domain")], by = "dataset_id", all.x = TRUE)
    ab <- merge(ab, pk[pk$fitted == kk, c("dataset_id", "nll_half", "b_inside", "domain")], by = "dataset_id", all.x = TRUE,
                suffixes = c("_1", "_k"))
    ab$nll_half_1[is.na(ab$domain_1) | is.na(ab$domain_k) | ab$domain_1 != ab$domain_k] <- NA
  } else ab$nll_half_1 <- ab$nll_half_k <- NA_real_
  ab <- ab[is.finite(ab$rmse_d1_1) & is.finite(ab$rmse_d1_k), ]
  differ <- abs(ab$rmse_d1_k / ab$rmse_d1_1 - 1) > 0.02
  better_k <- ab$rmse_d1_k < ab$rmse_d1_1
  for (cr in c("objective", "posterior_objective", "nll_half")) {
    v1 <- ab[[paste0(cr, "_1")]]; vk <- ab[[paste0(cr, "_k")]]
    ok <- is.finite(v1) & is.finite(vk)
    choose_k <- vk < v1
    right <- (choose_k == better_k)[ok & differ]
    picked <- ifelse(choose_k, ab$rmse_d1_k, ab$rmse_d1_1)[ok]
    pick_rows[[length(pick_rows) + 1L]] <- data.frame(k = kk, criterion = cr, datasets = sum(ok), differ = sum(ok & differ),
      right = mean(right), picks_k = mean(choose_k[ok]), k_better = mean(better_k[ok & differ]),
      ratio_pick_gate = med(picked / ab$rmse_d1_k[ok]), pass = mean(right) >= 0.70 && med(picked / ab$rmse_d1_k[ok]) <= 1)
  }
}
pick <- do.call(rbind, pick_rows)
utils::write.table(pick, file.path(out, "pick.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\n== pick step (raised v7 datasets; default fit against xk fit)")
print(pick, row.names = FALSE, digits = 3)
if (!fail_k) {
  pk_k <- pick[pick$k == k, ]
  say("pick step for k = ", k, ": ", if (any(pk_k$pass)) paste("passes with", paste(pk_k$criterion[pk_k$pass], collapse = ", ")) else
    "no criterion passes -> the gate alone decides")
}
vmp <- vm[vm$raised & vm$fitted %in% c(1, if (!fail_k) k else 4) & is.finite(vm$err), ]
if (nrow(vmp)) {
  say("varyMisseg (reported): share of raised inputs where the higher Rxv (KARST CV) / lower objective picks the more accurate fit")
  kk <- if (!fail_k) k else 4
  a <- vmp[vmp$fitted == 1, ]; b <- vmp[vmp$fitted == kk, ]
  ab <- merge(a, b, by = c("dataset_id", "minobs"), suffixes = c("_1", "_k"))
  bk <- ab$err_k < ab$err_1
  say("  Rxv: ", signif(mean((ab$Rxv_k > ab$Rxv_1) == bk, na.rm = TRUE), 3), "  objective: ",
      signif(mean((ab$objective_k < ab$objective_1) == bk, na.rm = TRUE), 3), "  (inputs x minobs: ", nrow(ab), ")")
}

## verdict ------------------------------------------------------------------------------
fail <- fail_false || fail_k || fail_vm
say("\nStep 1 ", if (fail) "FAILS" else "PASSES", if (!fail) paste0(" with k = ", k) else "")
writeLines(log_lines, file.path(out, "step1.txt"))

## figure: closer to the truth --------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  cl$truth_f <- factor(paste0("truth ×", cl$truth), levels = paste0("truth ×", multiples))
  g <- ggplot(cl, aes(x = fitted, y = ratio, colour = truth_f)) + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_line() + geom_point(size = 1.5) + scale_x_log10(breaks = multiples) + facet_wrap(~set, scales = "free_y") +
    labs(x = "fitted rate (multiple of 5e-5)", y = "median error ratio against x1", colour = NULL,
         caption = "v7: RMSE at distance 1; varyMisseg: 1 - r over KARST's landscape.") + theme_bw(base_size = 9)
  ggsave(file.path(out, "step1.pdf"), g, width = 8, height = 3.5)
}
