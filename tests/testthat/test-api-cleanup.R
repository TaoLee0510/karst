test_that("legacy exports are absent", {
  ns <- getNamespaceExports("karst")
  expect_false(paste0("run_", "karst_", "bench", "mark") %in% ns)
  retired_exports <- c(
    paste0("al", "fak"),
    paste0("predict", "_", "evo"),
    paste0("al", "fak_", "cohort_", "transition"),
    paste0("simulate_", "to", "y_", "landscape"),
    # The two-layer and joint-count pipelines, removed in 1.0.0.
    "fit_local_posterior", "fit_extrapolation_layer", "fit_graph_posterior", "build_karyotype_graph",
    "summarize_karst", "plot_karst"
  )
  retired_patterns <- paste(
    paste0("nn", "_", "prior"),
    paste0("cohort", "_", "transition"),
    paste0("K", "rig"),
    sep = "|"
  )
  expect_false(any(retired_exports %in% ns))
  expect_false(any(grepl(retired_patterns, ns)))
})
