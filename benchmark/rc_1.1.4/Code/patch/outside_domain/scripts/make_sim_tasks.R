# Task table of the simulated inputs (PLAN.md step 1): the main ABM set's 400
# simulations (fig2d_truth step 2 fits) and ABM_varyMisseg (fig2d_misseg,
# beta_fixed 125, beta_true 100), three designs each; `half` by replicate
# parity (odd: development, even: confirmation).
#   Rscript scripts/make_sim_tasks.R OUT_FILE
args <- commandArgs(TRUE)
w <- c(0.2, 0.4, 0.8, 1.6); designs <- c("w2_ends", "w4_ends", "w8_ends")
main <- expand.grid(rep = 1:100, w = w)
main$sim_id <- sprintf("w_%s_m_5e-05_rep_%02d", sub(".", "p", format(main$w, nsmall = 1), fixed = TRUE), main$rep)
rates <- c(5e-5, 1e-4, 2e-4, 4e-4, 8e-4)
vm <- expand.grid(rep = 1:25, pm = rates)
vm$sim_id <- sprintf("w_0p8_m_%s_rep_%02d", format(vm$pm, scientific = TRUE), vm$rep)
rows <- rbind(
  data.frame(source = "abm", arm = "step2", sim_id = rep(main$sim_id, each = 3), design = designs, rep = rep(main$rep, each = 3)),
  data.frame(source = "misseg", arm = "beta_fixed", sim_id = rep(vm$sim_id, each = 3), design = designs, rep = rep(vm$rep, each = 3)),
  data.frame(source = "misseg", arm = "beta_true", sim_id = rep(vm$sim_id[vm$pm != 5e-5], each = 3), design = designs,
             rep = rep(vm$rep[vm$pm != 5e-5], each = 3)))
rows$half <- ifelse(rows$rep %% 2 == 1, "dev", "confirm")
utils::write.table(rows, args[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
print(table(rows$source, rows$arm)); print(table(rows$half))
