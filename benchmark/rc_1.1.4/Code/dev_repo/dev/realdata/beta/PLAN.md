# Determining β on real data: plan and criteria written in advance

The user confirmed after 2026-11-20 01:28 EST: "Confirmed, run automatically according to the plan, stop after β2". This file was written at 2026-11-20 01:34 EST, before any β fit. It is not changed afterwards; anything that needs correcting goes into the "Corrections" section at the end, with the time and the reason.

β is the missegregation rate per division and per chromosome copy. The package default is 5e-5.

## Rules

- **Code and library:** everything goes in `dev/realdata/beta/`, no commit. The package is not changed; cd342db is used (`lib/karst_cd342db`).
- **HPC:** xxlarge, no cap, one dataset per task; the add-hours monitor uses small QOS, 3 days.
- **Simulations:** GRF only, no RMF (user decision 2).
- **Patient data:** only summary numbers appear in the report.
- **Stop:** stop once β2 is done and show the results to the user. Whether it goes into the package is for the user to decide.

## β0 test bed

**β grid:** 5e-5 × m, m ∈ {1/4, 1/2, 1, 2, 4, 8, 16, 32, 64}, 9 values in total, adjacent levels a factor of 2 apart.

| Set | Truth | Content |
|---|---|---|
| e1 | ×1 | passage culture, bottleneck 200 / 2000 / 20000, one time point every 6 or 60 days, 200 or 1000 cells, 4 or 8 time points (24 grid cells × 12 landscapes) |
| e2 | ×1 | clone replacement |
| e3 | ×1 | neutral culture (all karyotypes have the same true fitness) |
| e4 | odd landscapes ×4, even ×16 | the design of e1 |
| C2 grid | ×3, ×10, ×30 | the 6 settings of C2 round 1: β_sim ∈ {1.5e-4, 5e-4, 1.5e-3} × bottleneck {300, 3000}, 48 datasets each. One extra fit each at β_sim, for reference only |

- **Fitting:**
  - e sets: the first and last time points, as in `dev/R/eset.R`.
  - C2: as in `CE/c2_fit.R`.
- **Scoring:** the truth is each set's net growth rate at the true missegregation rate; scoring as in the original evaluation.
- **Real data:**
  - Salehi 277 trajectories (first and last passage, the setting of the paper benchmark); for every β, results are written in alfak format (minobs 5, 10, 20);
  - NUMBAT 43 inputs with the diploid added back, using the central division rate 0.141.

## β1 candidate estimates

Every fit records three things:
- **Marginal likelihood:** `diagnostics$objective`.
- **Data term:** the Poisson negative log-likelihood of the two time points on the model-domain states at the posterior mode, i.e. −Σ log Pois(y | expected).
- **Held-out prediction:**
  - the non-major states in the model domain are split at random into 5 folds with the dataset's seed;
  - for each fold, its count weights are set to 0 and the posterior mode is found with `sparse_joint_rebuild` (hyperparameters fixed at that fit);
  - take the expected counts of this fold's states at the two time points and compute the Poisson negative log-likelihood of their actual counts;
  - sum over the 5 folds. Majors are always kept.

| Estimate | Which β is taken |
|---|---|
| E1 | smallest marginal likelihood |
| E2 | smallest data term |
| E3 | smallest negative log-likelihood of the held-out prediction |
| E4 (Salehi only) | highest per-cell log score of the later-passage prediction (see below) |

**Two versions:**
- **Single:** each dataset chooses its own β.
- **Pooled:** the datasets of one group add up their criteria and choose one β together. How each set is grouped:

  | Set | Pooled group |
  |---|---|
  | e sets | same set, same design cell, same truth |
  | C2 | same setting |
  | Salehi | same cell line |
  | NUMBAT | same cancer type |

**E4:**
- The alfak results of each β are used to predict the later passages with `validation/sparse_joint/salehi_horizon/forecast_score.R` unchanged. The missegregation rate of the prediction is the same β.
- The per-cell log score on the common support, `ll_fc_s_e2`, is used. The common support is the union of the landscapes of all β for the same trajectory, computed on later passages ≤ 60 days.
- E4 has only the pooled version, and it is leave-one-out: each trajectory's β is chosen by the other trajectories of the same cell line, and this trajectory's own later passages do not take part.

## β2 verdict

**Split:**
- **Development:** landscapes L01–L06 of the e sets and C2, the Salehi development group (SA609, SA039).
- **Verdict:** L07–L12 of the e sets and C2, the Salehi confirmation group (SA906, SA532, SA535, SA1035).
- GBM and the other NUMBAT inputs are only reported.

**Procedure:**
- Each estimate × version is first checked against all criteria on the development set; only those that meet them enter the verdict set, where they are checked once more.
- Only meeting them both times counts as a pass.
- E4 has no simulations; only the Salehi criteria apply.

**Criteria (all must be met):**
- **N null hypothesis** (e1–e3, truth ×1; RMSE ratio = RMSE at the chosen β ÷ RMSE at the true β):
  - e1: take the dataset median per grid cell; ≥ 18 of the 24 grid cells are ≤ 1.03 at all of d0, d1, d2, and for each distance the median of the grid-cell medians is ≤ 1.02;
  - e2, e3: the median at each distance ≤ 1.03.
  - The development and verdict sets each have only 6 landscapes per grid cell; computed as above.
- **A1 chooses accurately (e4, C2):**
  - single version: ≥ 70% of choices within a factor of 2 of the truth;
  - pooled version: ≥ 80% of groups within a factor of 2 of the truth.
  - "Within a factor of 2" means |log₂(choice / truth)| ≤ 1. The truth for C2 is the actual β_sim (about ×3, ×10, ×30).
- **A2 more accurate (e4, C2):** for each truth, the dataset median of the RMSE ratio (chosen β ÷ fixed 5e-5) at d0, d1, d2 is < 1 when the truth is ≥ 4 times, and ≤ 1.02 for C2's ×3.
- **S1 Salehi later passages:** on later passages ≤ 60 days, the paired difference of `ll_fc_s_e2` between the chosen β and 5e-5 has mean ≥ 0.
- **S2 Salehi parallel lineages (single version only):** among trajectory pairs with the same cell line, the same treatment state and no overlapping input passages (as in B), in ≥ 70% of pairs the chosen β differ by at most a factor of 2.

**Report only:**
- the β chosen by each estimate on each Salehi cell line and each NUMBAT cancer type;
- the RMSE of C2 at the true β_sim;
- the curve of each criterion against β.

## Stopping rules

- a build or script fails for an unclear reason;
- > 5% of the tasks of an array fail, and it is neither a resource problem nor a model-domain error of the package itself;
- a situation the rules do not cover.

After all verdicts are made, write `report.md` and stop.

## Corrections
- **2026-11-20 01:34 EST (clarification, before any fit):** E4 and S1 use the minobs 10 alfak results, the same as W3 of D3. The prediction settings nb, n0 etc. are the same as in B.
