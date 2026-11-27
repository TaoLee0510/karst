# Diagnostic: ALFA-K's rate map on list E against its own landscape, for a few
# Kester inputs at x1 (why its engine score is low). Descriptive only.
#   Rscript diag_alfak.R TASK_FILE ID OUT_DIR
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (f in c("common.R", "engine.R", "support.R", "kester.R")) source(file.path(here, "R", f))
out <- args[[3]]; k <- od_karst(use_package = TRUE)
canon <- function(v) k$format_karyotypes(k$parse_karyotypes(v))
for (id in c("APKS2_P5-14_to_P18_ends", "APKS3_P3-15_to_P19_ends", "APKS3_P3-7-15_to_P19_all", "AP1_cells_1")) {
  inp <- readRDS(file.path(out, "inputs", paste0(id, ".rds")))
  E <- fr_list_e(unique(canon(rownames(inp$yi$x))))
  first <- if (inp$kind == "salehi") inp$parent else { v <- inp$train_full[, 1]; v[v > 0] }
  for (mo in c(5, 10)) {
    dir <- file.path(out, "fits", id, "x01", "alfak", paste0("minobs_", mo))
    if (!file.exists(file.path(dir, "landscape.Rds"))) next
    lsc <- readRDS(file.path(dir, "landscape.Rds"))
    r <- kester_alfak_on(dir, E, k)
    own <- r[match(lsc$k, E)]
    obs <- r[match(names(first), E)]
    kf <- readRDS(file.path(out, "fits", id, "x01", "karst", "minobs_5", "karst_fit.Rds"))
    kr <- kester_karst_rates(kf, E, k)
    cat(sprintf("%s m%d | E %d, own list %d | ALFA-K rates: own list %.3f..%.3f (median %.3f), E %.3f..%.3f (median %.3f), parent karyotypes %.3f..%.3f | share of E above the parent max: %.3f | KARST rates E %.3f..%.3f, parent %.3f..%.3f\n",
                id, mo, length(E), nrow(lsc), min(own, na.rm = TRUE), max(own, na.rm = TRUE), median(own, na.rm = TRUE), min(r), max(r), median(r),
                min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), mean(r > max(obs, na.rm = TRUE)), min(kr), max(kr), min(kr[match(names(first), E)]), max(kr[match(names(first), E)])))
  }
}
