# Kester (benchmark/rc_1.1.4/PLAN.md, part 3): build every fit input, the
# fold manifest and the task table. Run once on the HPC (data stay on RED).
#   Rscript setup.R TASK_FILE ID OUT_DIR      (TASK_FILE, ID unused; array of one)
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (f in c("common.R", "engine.R", "support.R", "kester.R")) source(file.path(here, "R", f))
out <- args[[3]]; dir.create(file.path(out, "inputs"), recursive = TRUE, showWarnings = FALSE)
cult <- kester_inputs()
manifest <- list()
rebase <- function(x) { p <- as.numeric(colnames(x)); colnames(x) <- as.character(p - p[1]); x }

## Salehi-like folds ---------------------------------------------------------------------
for (i in seq_len(nrow(kester_salehi_folds))) {
  f <- kester_salehi_folds[i, ]
  x <- kester_read(cult$path[cult$culture == f$culture])
  tr <- as.character(strsplit(f$train, ",")[[1]]); tg <- as.character(f$target)
  stopifnot(all(c(tr, tg) %in% colnames(x)))
  xt <- x[, tr, drop = FALSE]; xt <- xt[rowSums(xt) > 0, , drop = FALSE]
  last <- xt[, ncol(xt)]; parent <- last[last > 0]
  earlier <- rownames(xt)[rowSums(xt[, -ncol(xt), drop = FALSE]) > 0]
  eps <- max(0.01, sum(parent[!names(parent) %in% earlier]) / sum(parent))
  target <- x[, tg]; target <- target[target > 0]
  id <- sprintf("%s_P%s_to_P%s_%s", f$culture, gsub(",", "-", f$train), tg, f$design)
  inp <- list(id = id, kind = "salehi", culture = f$culture, design = f$design, train_passages = as.numeric(tr), target_passage = as.numeric(tg),
              yi = list(x = rebase(xt), dt = kester_days_per_passage), parent = parent, target = target,
              horizon_days = (as.numeric(tg) - as.numeric(tr[length(tr)])) * kester_days_per_passage, eps = eps)
  saveRDS(inp, file.path(out, "inputs", paste0(id, ".rds")))
  manifest[[length(manifest) + 1L]] <- data.frame(input = id, kind = "salehi", culture = f$culture, design = f$design, train = f$train,
    target = tg, train_cells = paste(colSums(xt), collapse = ","), target_cells = sum(target), horizon_days = inp$horizon_days,
    eps = eps, karyotypes_train = nrow(xt), stringsAsFactors = FALSE)
}

## GBM-like: full fit and seeded three-fold cell holdout -----------------------------------
for (cu in cult$culture[cult$kind == "gbm"]) {
  x <- kester_read(cult$path[cult$culture == cu])
  stopifnot(ncol(x) == 2L)
  folds <- c(list(list(id = "full", train = x, test = NULL)), kester_cell_folds(x))
  for (fo in folds) {
    xt <- fo$train[rowSums(fo$train) > 0, , drop = FALSE]
    id <- sprintf("%s_%s", cu, fo$id)
    inp <- list(id = id, kind = "gbm", culture = cu, fold = fo$id, passages = as.numeric(colnames(x)),
                yi = list(x = rebase(xt), dt = kester_days_per_passage), train_full = fo$train, test = fo$test,
                horizon_days = diff(as.numeric(colnames(x))) * kester_days_per_passage)
    saveRDS(inp, file.path(out, "inputs", paste0(id, ".rds")))
    manifest[[length(manifest) + 1L]] <- data.frame(input = id, kind = "gbm", culture = cu, design = fo$id, train = paste(colnames(x), collapse = ","),
      target = if (is.null(fo$test)) "in-sample" else "held-out cells", train_cells = paste(colSums(fo$train), collapse = ","),
      target_cells = if (is.null(fo$test)) NA else paste(colSums(fo$test), collapse = ","), horizon_days = inp$horizon_days,
      eps = NA, karyotypes_train = nrow(xt), stringsAsFactors = FALSE)
  }
}
manifest <- do.call(rbind, manifest)
utils::write.table(manifest, file.path(out, "fold_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## task table: one fit per task -----------------------------------------------------------
tasks <- do.call(rbind, lapply(manifest$input, function(id) do.call(rbind, lapply(kester_grid, function(m)
  data.frame(input = id, multiple = m, method = c("karst", "alfak", "alfak", "alfak"), minobs = c("5,10,20", "5", "10", "20"),
             stringsAsFactors = FALSE)))))
tasks$task <- seq_len(nrow(tasks))
utils::write.table(tasks, file.path(out, "tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("karst", as.character(utils::packageVersion("karst")), find.package("karst")),
             paste("alfakR", as.character(utils::packageVersion("alfakR")), find.package("alfakR")),
             paste("seed", kester_seed), paste("days per passage", kester_days_per_passage), paste("multiples", paste(kester_grid, collapse = ","))),
           file.path(out, "versions.txt"))
print(manifest); cat("tasks", nrow(tasks), "\n")
