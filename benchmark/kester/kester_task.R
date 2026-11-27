# Kester (benchmark/rc_1.1.4/PLAN.md, part 3): one fit task. KARST (all three
# minobs in one task) or ALFA-K (one minobs), at missegregation rate
# 5e-5 x multiple, with the paper benchmark's own functions
# (R/paper/{config,fit,salehi}.R, copied from benchmark/paper/R). Salehi-like
# inputs also get the paper's ABM forecast (5 replicates) to the held-out
# passage. Writes OUT/fits/<input>/x<mm>/<method>/minobs_<m>/.
#   Rscript kester_task.R TASK_FILE ID OUT_DIR
args <- commandArgs(TRUE)
here <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (f in list.files(file.path(here, "R", "paper"), full.names = TRUE)) source(f)
tasks <- utils::read.delim(args[[1]], stringsAsFactors = FALSE, colClasses = c(minobs = "character"))
task <- tasks[tasks$task == as.integer(args[[2]]), ]; stopifnot(nrow(task) == 1L)
out <- args[[3]]
inp <- readRDS(file.path(out, "inputs", paste0(task$input, ".rds")))
pm <- 5e-5 * task$multiple; n0 <- 1e5; nb <- 1e7
outdir <- file.path(out, "fits", task$input, sprintf("x%02d", task$multiple), task$method)
minobs <- as.numeric(strsplit(task$minobs, ",")[[1]])
cat("Task", args[[2]], ":", paste(unlist(task), collapse = " | "), "\n")
for (p in c("alfakR", "karst")) cat(p, as.character(utils::packageVersion(p)), find.package(p), "\n")
cat("Input time points (days):", input_times(inp$yi), " pm:", pm, "\n")
todo <- minobs[!vapply(minobs, function(m) file.exists(file.path(outdir, paste0("minobs_", m), "run.Rds")), logical(1))]
if (length(todo)) {
  if (task$method == "alfak") for (m in todo) fit_alfak(inp$yi, outdir, m, n0, nb, pm) else
    fit_karst(inp$yi, outdir, todo, n0, nb, pm, min_cn = 0)
}
for (m in minobs) {
  dir <- file.path(outdir, paste0("minobs_", m)); r <- readRDS(file.path(dir, "run.Rds"))
  cat(sprintf("Fit minobs %s: %s (%.0f s) %s\n", m, r$status, r$seconds, ifelse(is.na(r$message), "", r$message)))
  if (inp$kind != "salehi" || !identical(r$status, "ok") || file.exists(file.path(dir, "forecast.Rds"))) next
  start <- proc.time()[["elapsed"]]
  msg <- tryCatch({ forecast_salehi_dir(dir, inp$yi, inp$horizon_days, pm, nb); NA_character_ }, error = function(e) conditionMessage(e))
  save_rds_atomic(list(status = if (is.na(msg)) "ok" else "failed", seconds = proc.time()[["elapsed"]] - start, message = msg,
                       horizon_days = inp$horizon_days, pm = pm), file.path(dir, "forecast.Rds"))
  cat(sprintf("Forecast minobs %s: %s %s\n", m, if (is.na(msg)) "ok" else "failed", ifelse(is.na(msg), "", msg)))
}
