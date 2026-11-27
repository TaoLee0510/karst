# rc freeze: GATE_RESULTS is the gate's results folder (inputs read from it,
# outputs written to it); the as-run copy read and wrote
# validation/sparse_joint/dispersion_gate/results from the repository root.
G <- Sys.getenv("GATE_RESULTS"); stopifnot(nzchar(G))
td <- commandArgs(TRUE)[1]
a <- rbind(read.delim(file.path(td, "tasks_beta_fixed_pilot.tsv"), stringsAsFactors = FALSE, colClasses = "character"),
           read.delim(file.path(td, "tasks_beta_fixed_rest.tsv"), stringsAsFactors = FALSE, colClasses = "character"))
k <- a[a$method == "karst", ]
cat("karst tasks:", nrow(k), " sims:", length(unique(k$dataset)), "\n")
for (m in c(2, 4, 8, 16)) {
  t <- k; t$pm <- format(5e-5 * m, scientific = TRUE); t$task <- seq_len(nrow(t))
  write.table(t, sprintf(file.path(G, "step1/tasks/fits_vm_x%02d.tsv"), m), sep = "\t", quote = FALSE, row.names = FALSE)
}
print(head(t, 3))
