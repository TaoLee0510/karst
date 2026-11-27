# One interface to the datasets of every evaluation set, so that a method can be
# fitted and scored the same way on each:
#   kind "e1".."e4"        e sets (first and last time point), net-rate truth
#   kind "dev", "heldout"  development benchmark GRF datasets (scoreboard)
#   kind "rmf", "nk"       D1 families (dev/R/families.R)
# dev_adapter(kind, id) returns args (for karst::fit_sparse_joint, without
# method-specific arguments), nodes, truth (aligned with nodes), meta.
# Sourced after common.R and eset.R.

dev_adapter <- function(kind, id) {
  if (kind %in% c("e1", "e2", "e3", "e4")) {
    dev_source_validation()
    loaded <- dev_eset_load(kind, id)
    d <- loaded$dataset
    counts <- mt_counts(d, "ends")
    return(list(kind = kind, id = id, dataset = d, nodes = loaded$nodes,
                args = list(counts, times = d$times[c(1L, length(d$times))], beta = 5e-5, birth_rate = 1,
                            min_cn = 1L, max_cn = 8L),
                truth = sm_true_net_rate(loaded$nodes$karyotype, d$truth, d$beta),
                meta = list(bottleneck = d$bottleneck, spacing = d$spacing, depth = d$depth, samples = d$samples)))
  }
  if (kind %in% c("dev", "heldout")) {
    dev_devbench_source()
    input <- dev_devbench_input()
    d <- readRDS(file.path(input, "datasets", paste0(id, ".rds")))
    nodes <- readRDS(file.path(input, "nodes", paste0(id, ".rds")))
    return(list(kind = kind, id = id, dataset = d, nodes = nodes,
                args = list(d$counts, dt = d$interval, beta = d$beta, birth_rate = 1, min_cn = 1L, max_cn = 8L,
                            initial_state = "balance"),
                truth = stt_true_growth_rate(nodes$karyotype, d$truth),
                meta = list(scenario = d$scenario, wavelength_class = d$wavelength_class)))
  }
  if (kind %in% c("diverse", "diverse2", "paper3")) {
    # cloud_level sets (calling errors: dataset$call_error)
    dev_devbench_source()
    input <- {
      local <- file.path(dev_inputs, "cloud", kind)
      if (dir.exists(local)) local else file.path(dev_backup, "cloud_level", "results", "inputs", kind)
    }
    d <- readRDS(file.path(input, "datasets", paste0(id, ".rds")))
    nodes <- readRDS(file.path(input, "nodes", paste0(id, ".rds")))
    return(list(kind = kind, id = id, dataset = d, nodes = nodes,
                args = list(d$counts, dt = d$interval, beta = d$beta_assumed, birth_rate = 1, min_cn = 1L,
                            max_cn = 8L, initial_state = "balance"),
                truth = stt_true_growth_rate(nodes$karyotype, d$truth),
                meta = list(level = d$level, mechanism = d$mechanism, call_error = d$call_error,
                            scenario = d$scenario)))
  }
  if (kind == "regrowth") {
    # patient_regrowth (validation/sparse_joint/patient_regrowth): truth per
    # calendar day, the growing rate times the growing share
    dev_devbench_source()
    input <- {
      local <- file.path(dev_inputs, "regrowth")
      if (dir.exists(local)) local else file.path(dev_backup, "patient_regrowth", "results")
    }
    d <- readRDS(file.path(input, "datasets", paste0(id, ".rds")))
    nodes <- readRDS(file.path(input, "nodes", paste0(id, ".rds")))
    return(list(kind = kind, id = id, dataset = d, nodes = nodes,
                args = list(d$counts, dt = d$interval, beta = d$beta, birth_rate = 1, min_cn = 1L, max_cn = 8L,
                            initial_state = "balance"),
                truth = stt_true_growth_rate(nodes$karyotype, d$truth) * d$growing,
                meta = list(scenario = d$scenario, survivors = d$survivors, wavelength_class = d$wavelength_class)))
  }
  if (kind == "misseg") {
    # ALFA-K paper missegregation simulations; id = "<sim>__<design>" (ends
    # designs). Fitted at the default rate 5e-5 (the R2 grid starts there);
    # truth: the net rate at the simulation's true rate. Nodes as in the
    # development benchmark (stt_scoring_nodes on first and last counts).
    if (!exists("dev_paper_input", mode = "function")) source(file.path(dev_repo, "dev", "R", "paper_abm.R"))
    dev_devbench_source()
    parts <- strsplit(id, "__", fixed = TRUE)[[1L]]
    input <- dev_paper_input("misseg", parts[[1L]], parts[[2L]])
    x <- input$x
    collapsed <- cbind(x[, 1L], x[, ncol(x)])
    rownames(collapsed) <- rownames(x)
    seed <- sum(utf8ToInt(id)) %% 100000L
    nodes_path <- file.path(dev_results, "phase1", "misseg_nodes", paste0(id, ".rds"))
    nodes <- if (file.exists(nodes_path)) readRDS(nodes_path) else {
      nd <- stt_scoring_nodes(collapsed, 1L, 8L, 2000L, 2000L, seed)
      dev_save_rds(nd, nodes_path)
      nd
    }
    cn <- karst:::parse_karyotypes(nodes$karyotype)
    truth <- grf_fitness(nodes$karyotype, input$info$w, input$sim$true_landscape) *
      (2 * karst:::sparse_joint_self_probability(cn, input$info$pmis) - 1)
    return(list(kind = kind, id = id, dataset = list(seed = seed), nodes = nodes,
                args = list(x, dt = diff(input$times), beta = 5e-5),
                truth = truth, meta = list(pmis = input$info$pmis, w = input$info$w, design = parts[[2L]])))
  }
  if (kind %in% c("rmf", "nk")) {
    if (!exists("dev_family_config", mode = "function")) source(file.path(dev_repo, "dev", "R", "families.R"))
    config <- dev_family_config(kind)
    d <- dev_family_dataset(config, id)
    nodes <- dev_family_nodes(config, d)
    counts <- d$counts
    storage.mode(counts) <- "integer"
    return(list(kind = kind, id = id, dataset = d, nodes = nodes,
                args = list(counts, dt = d$interval, beta = d$beta, birth_rate = 1, min_cn = 1L, max_cn = 8L,
                            initial_state = "balance"),
                truth = dev_family_node_truth(config, d, nodes),
                meta = list(scenario = d$scenario)))
  }
  stop("Unknown kind ", kind)
}

# Score a prediction (data frame from predict()) of an adapter's nodes. Returns
# the metrics and the scored nodes with the fit's model sd, own relative
# fitness and shell (what the calibration formula needs).
dev_adapter_score <- function(a, prediction) {
  if (!exists("stt_score_predictions", mode = "function")) dev_devbench_source()
  pred <- data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                     fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                     shell = prediction$shell, stringsAsFactors = FALSE)
  score <- stt_score_predictions(a$nodes, a$truth, pred, 5)
  table <- score$nodes
  at <- match(table$karyotype, pred$karyotype)
  table$model_sd <- pred$model_sd[at]
  table$shell <- pred$shell[at]
  table$relative <- pred$fitness[at]
  list(metrics = score$metrics, nodes = table)
}
