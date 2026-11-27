# Overcoming KARST's limitations: overall autorun report (phases 0–3)

Autorun approved by the user on 2026-11-19.

- **Criteria:** written in `PLAN.md` and in its Clarifications 2–5, all fixed before scoring and never changed.
- **Rules:** no commits during the run; the package code and defaults are not changed; all model changes are patches in `patches/`, built into a separate library; HPC uses xxlarge, with no task limit.
- **Reports:** details of each phase in `phase0/report.md` … `phase3/report.md`; progress and incidents in `README.md`.

## Conclusion

**No candidate passed; no change to the package is recommended.**

| Phase | Candidate | Content | Verdict |
|---|---|---|---|
| 1 | R1 | drift-aware second-time-point likelihood (NB1, F computed from the passage records) | not passed (all five) |
| 1 | R2(i) | put the hyperparameter uncertainty back into the sd | not kept (coverage almost unchanged) |
| 1 | R2(iii) | mixture over missegregation rates ×1/×4/×16 | not kept (width 1.9 times the current one, limit 1.25) |
| 1 | R4 | symmetric karyotype-miscall observation kernel | not passed (all three) |
| 2 | R3 | aggregated count likelihood for the second shell | not passed; out of memory on large datasets even with 256 GB |
| 2 | R9 | R3 + mean-reverting decay | not passed; helps only RMF, GRF λ 0.4 gets worse |
| 3 | R5 | random rare-karyotype cloud | not passed (RMSE passed, false positives and power not passed) |
| 3 | R6 | β prior | not passed |
| 3 | R6 | chromosome dosage prior | not passed |
| 3 | R7 | history fixed to the recorded population age | not passed |
| 3 | R8(a) | interval calibration by per-dataset parametric bootstrap | not passed (under-coverage at far distances on GRF) |
| 3 | R8(b) | cross-family conformal | not passed (scale difference between NK and GRF/RMF too large) |
| 3 | R10 | WGD; piecewise fitness | not passed: the WGD part is not evaluable because of a problem in the simulation specification, the piecewise fitness part did not pass |

## What was learned

1. **The interval shortfall comes from model misspecification, not from the approximate computation.** Four lines of evidence:
   - D5: Laplace's sd is only about 15% smaller than that of NUTS.
   - R8(b): to reach 95% coverage, the model sd has to be multiplied by 6–15, not by 1.96.
   - R8(a): a bootstrap on the fitted model inherits the same misspecification, and coverage at far distances is still insufficient.
   - R2: putting the hyperparameter uncertainty back has almost no effect.
2. **Calibration constants do not transfer between landscape families.** (On 2026-11-20 the user pointed out that karyotype fitness landscapes are not like RMF; RMF is not used after that. The verdicts below that involve RMF are all unchanged when only GRF is considered.) D1, R8(b) and R9 all show that constants or mean-reversion strengths that suit GRF do not suit NK or RMF.
3. **The upward bias of β comes from the prior and the Laplace determinant, not from the data or drift.** D3 shows that the data term itself favours the true β; none of R1, R5, R6 reduced the upward bias.
4. **Drift holds only for the major karyotypes.**
   - The variance of the major karyotypes agrees with the drift prediction; the variance of rare karyotypes is only 0.47 times the prediction, because missegregation inflow keeps replenishing them (D4).
   - Using the same F for all karyotypes (R1) weakens the constraint that rare karyotypes place on fitness.
   - R5 corrects this, but the model sd becomes smaller with it, and false positives increase.
5. **Far-distance ranking is far from the upper bound, and none of the changes improved it.** KARST reaches only 0.27–0.36 times the upper bound at λ 0.4 and 0.8 (D2). R3 improves ranking slightly only at λ 1.2 (distance 3 Spearman 0.28 → 0.32, pooled +0.001, interval includes 0); R9's decay suits only RMF.
6. **Choosing history by the marginal likelihood is more robust than fixing it to the recorded population age (R7).**
7. **Problems with the specifications themselves (recorded as they are; the rules were not changed):**
   - R7's ABM thresholds are too high: base itself reaches only 0.54 / 0.44, while the thresholds are 0.75 / 0.66.
   - R10's WGD simulation specification had the direction of the adjustment reversed: the gain should have been reduced, not increased, so the WGD part is not evaluable. The results, for reference only, are positive: SA906a per-cell log score +0.17 to +0.25.

## Relation to the real-data plan

The autorun was evaluated entirely on simulations, and neither the candidates nor the simulations were aimed at the problems of real data; R7 and R10 are examples. The subsequent real-data plan (`realdata/`, confirmed by the user 2026-11-19 14:25) instead diagnoses on real data first. The results of B (`realdata/B/report.md`) agree with the points above:
- The interval problem on real data cuts both ways: the reported sd is too wide and uninformative, while the model sd is too narrow.
- The symmetric miscall kernel (R4) does not match the one-directional missed calls in the real NUMBAT data.
- WGD matters only on a few trajectories of SA906 and SA609.

## Incidents during the run (details in `README.md`)

- **Wrong times:** the times stated in README and PLAN before 08:24 were originally estimates and have been corrected according to the session log; the ordering of the criteria is not affected.
- **Compilation cache conflict:** several tasks shared the Rcpp compilation cache, causing some fits of R8(a), base and R9 to fail; they were rerun after switching to one cache per process.
- **Memory for R3 and R9:** out-of-memory tasks were rerun step by step at 8 GB → 32 GB → 256 GB; on 2026-11-20 the user decided not to rerun them any further.
- **The add-hours monitor stopped first:** the monitor job for the R3 and R9 reruns was placed in xxlarge and was cancelled after 12 hours. Once this was noticed, hours were added by hand and the monitor was reattached in the small QOS; no task was lost to a timeout.
- **A problem in the package itself (unrelated to the candidates):** 1.1.0 reports `karst_inner_failure` on 2–6% of the datasets in the e sets (phase 1 report). The later 1.1.1 (50415cf) searches for the posterior mode once more before raising the error; 1.1.2, used by the real-data plan, did not produce this error in the B fits.
