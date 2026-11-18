# Benchmark results by wavelength class, as LaTeX tables (doc/karst_method.tex,
# Section 5).
#
#   Rscript validation/sparse_joint/wavelength_tables.R METRICS_CSV OUT_TEX
#
# Values are averaged over the scenarios of each landscape and then over the
# three landscapes of each wavelength; core scenarios and patient 1 only. The
# wavelength classes are those in the metrics, in order of wavelength. A
# failed fit is scored as the null predictor, with rank correlation zero
# (stt_null_rank() in benchmark/sparse_two_timepoint/report.R).

source("benchmark/sparse_two_timepoint/report.R")
args <- commandArgs(TRUE)
m <- stt_null_rank(read.csv(args[[1]], stringsAsFactors = FALSE))
tt <- function(x) paste0("\\texttt{", gsub("_", "\\\\_", x), "}")
methods <- c("Null predictor" = "zero",
             "alfak" = "alfak_minobs20",
             "Sparse joint, independent" = "karst_sparse_joint",
             "Sparse joint, balance" = "karst_sparse_joint_balance")
display <- function(label) {
  if (label == "alfak") return("\\texttt{alfak} (\\texttt{minobs = 20})")
  if (grepl("_", label)) return(tt(label))
  label
}
x <- m[m$method %in% methods & m$role == "core" & m$scoreboard & m$stratum %in% c("d0", "d1", "d2"), ]
cols <- c("rmse", "spearman", "coverage95")
f <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
cell <- aggregate(x[cols], x[c("method", "wavelength_class", "landscape_id", "scenario", "stratum")], f)

land <- aggregate(cell[cols], cell[c("method", "wavelength_class", "landscape_id", "stratum")], f)
wl <- aggregate(land[cols], land[c("method", "wavelength_class", "stratum")], f)
waves <- unique(wl$wavelength_class)
waves <- waves[order(as.numeric(sub("^w", "", waves)))]
cat("Wavelength classes:", waves, "\n")
value <- function(method, w, d, metric) {
  v <- wl[wl$method == method & wl$wavelength_class == w & wl$stratum == d, metric]
  if (!length(v) || is.na(v)) "--" else sprintf("%.3f", v)
}
table_tex <- function(metric, keep, header) {
  k <- seq_along(waves)
  out <- c("\\begin{center}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
           sprintf("\\begin{tabular}{l*{%d}{r}}", 3L * length(waves)), "\\toprule",
           paste0(header, paste(sprintf(" & \\multicolumn{3}{c}{$\\lambda=%s$}", sub("^w", "", waves)), collapse = ""),
                  " \\\\"),
           paste(sprintf("\\cmidrule(lr){%d-%d}", 3L * k - 1L, 3L * k + 1L), collapse = ""),
           paste0(paste(rep(" & $d=0$ & $d=1$ & $d=2$", length(waves)), collapse = ""), " \\\\"), "\\midrule")
  for (label in keep) {
    vals <- unlist(lapply(waves, function(w) vapply(c("d0", "d1", "d2"), function(d) value(methods[[label]], w, d, metric), "")))
    out <- c(out, paste0(display(label), " & ", paste(vals, collapse = " & "), " \\\\"))
  }
  c(out, "\\bottomrule", "\\end{tabular}", "\\end{center}")
}
has <- function(metric, label) any(is.finite(wl[wl$method == methods[[label]], metric]))
# The null predictor gives no ranking (its rank correlation is zero by
# definition) and no intervals.
keep_rank <- Filter(function(l) has("spearman", l), setdiff(names(methods), "Null predictor"))
keep_cov <- Filter(function(l) has("coverage95", l), names(methods))
tex <- c("\\paragraph{Root mean squared error.}", table_tex("rmse", names(methods), ""), "",
         "\\paragraph{Spearman correlation with the reference values.}", table_tex("spearman", keep_rank, ""), "",
         "\\paragraph{Coverage of nominal 95\\% intervals.}", table_tex("coverage95", keep_cov, ""))
writeLines(tex, args[[2]])
failed <- aggregate(status ~ method + wavelength_class, unique(x[c("dataset_id", "method", "wavelength_class", "status")]),
                    function(s) mean(s != "ok"))
print(failed[failed$status > 0, ], row.names = FALSE)
