# D3: WGD (phases C–E)

Criteria: D3 in the "Phases C–E" section of `../PLAN.md` (fixed 2026-11-20 00:35 EST). W1–W3 are defined as in `dev/phase3/r10_spec.md`, and both values of ω (1e-4, 1e-3) must satisfy them. Scoring script `d3_eval.R`, results `dev/results/realdata/CE/d3/evaluation.rds`. The tables contain only summary numbers.

## Conclusion: not passed

At ω = 1e-3, W2 is not satisfied (the RMSE ratio at distance 0 is 1.031 > 1.01). By the rule both values of ω must satisfy the criteria, so D3 does not pass regardless of the W1 result.

| Criterion | ω = 1e-4 | ω = 1e-3 | Verdict |
|---|---|---|---|
| W1 C2 scenario W, scored karyotypes with mean copy number ≥ 3, median RMSE ratio ≤ 0.85 | pending | pending | C2's W scenario not yet accepted, not evaluated |
| W2 development GRF without WGD (96 in the dev split), median d0/d1/d2 RMSE ratio ≤ 1.01 | 1.003 / 1.001 / 1.002 | **1.031** / 1.014 / 1.018 | satisfied at ω = 1e-4, not satisfied at ω = 1e-3 |
| W2 C2 scenarios H, M ≤ 1.01 | pending | pending | C2 not yet accepted, not evaluated |
| W3 eligible SA906a ends trajectories, paired mean per-cell `ll_fc_e2` ≥ base | +0.173 (sd 0.40, 55% of rows better) | +0.246 (sd 0.51, 59% of rows better) | satisfied at both ω |

## Method

- **Patch** `dev/patches/D3_wgd.patch`: `R10_wgd.patch` moved onto cd342db with unchanged content (`fit_sparse_joint(wgd = ω)`: the doubled karyotypes of ploidy-2 majors enter the model domain, the single-daughter probability on the edge x → 2x is ω/2, the missegregation rate of ploidy-4 states is 2.5β, and doubled karyotypes take the parent's trend row, shell 1). cd342db changed the code near the return value of `fit_sparse_joint()`, so the one place that adds `wgd = settings$wgd` to the `trend` list had to be put back by hand; the other 8 places applied as is.
- **Build** HPC `karst_dev/lib/karst_d3` (job 22139515, passed, version 1.1.2).
- **Bit-for-bit identical to base without `wgd`**: on 6 datasets (4 development GRF, 2 RMF), D3 without `wgd` and the cd342db base (`ce_base`) have exactly the same near and far node tables and objective values (job 22139549).
- **W2**: `dev/phase2/far_fit.R`, methods `d3_w4` (ω = 1e-4, job 22139551) and `d3_w3` (ω = 1e-3, job 22139553), paired one by one with `ce_base` (22138400) on the same 96 datasets; all 96 fit successfully.
- **W3**:
  - Trajectories: the 12 SA906a ends trajectories in which cells with mean copy number ≥ 3 make up ≥ 20% of the last passage and which have descendants (the same list as for R10), minobs 10.
  - base: B's cd342db default fits and B's later-passage scores (same parameters: β 5e-5, `min_cn = 0`, nboot 500, n0 1e5, nb 1e7, correct_efflux).
  - D3: the same parameters plus `wgd` (`d3_salehi.R`, job 22139813), alfak format written out as in the paper benchmark, scored unchanged with `forecast_score.R` (`d3_horizon.R`, job 22140455).
  - Paired by (trajectory, descendant passage), 22 rows in total. The median landscape coverage is the same for both (0.906).

## Observations

- At ω = 1e-4, the WGD term does almost no harm on data without WGD (RMSE ratio ≤ 1.003); at ω = 1e-3, distance 0 gets about 3% worse. This is consistent with the reference results of phase 3 (1.003 / 1.000 / 1.000 and 1.032 / 1.010 / 1.011).
- On SA906a, the per-cell log score at both ω is higher than base (+0.17, +0.25), the same as the phase 3 reference results on b319cfd. But the standard deviations are large (0.40, 0.51), only 55–59% of rows are better, and the 12 trajectories share many overlapping passages, so there is little confidence in this improvement.

## Deviations and incidents

- **Scoring for W3**: `r10_spec.md` says "by the method of `forecast_score.R`". Phase 3's R10 used self-written scoring code and fits with nboot 50 and `max_cn` ≥ 10; this time, following PLAN, scoring uses `forecast_score.R` unchanged, the fit parameters are the same as B's base, and base is taken directly from B's results.
- **Task table written wrongly once**: when the W3 task table was first generated, zsh did not word-split the unquoted variable, so the method column became "d3_w4 1e-4" and ω was empty. All of the first fits (22139555) therefore errored; the horizon job that depended on them (22139558) had already started before it was cancelled, could not find fits according to the other table, and wrote 24 "no fit" outputs. All these outputs were deleted, and the table was generated with R and rerun: fits 22139813, horizon 22140455 (an intermediate run, 22139815, skipped everything because the old outputs were still there, and has also been discarded).
- **W1 and C2's W2 pending**: C2 is handled by another branch; at the time of writing this report (2026-11-20 00:53 EST) C2 was still at the tuning-grid stage and `results/realdata/CE/c2/regimes/` did not exist yet. These two items do not affect the conclusion (W2 is already not satisfied at ω = 1e-3). Once C2 is accepted, they can be evaluated afterwards with `c2_fit.R` (methods `d3_w4`/`d3_w3`, library `karst_d3`, `extra` set to `list(wgd = ω)`), for reference only.
