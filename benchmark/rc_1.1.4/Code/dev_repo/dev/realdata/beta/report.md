# Determining β on real data: β0–β2 results

**Basis:**
- The plan and criteria are in `PLAN.md`, fixed 2026-11-20 01:34 EST and not changed before scoring.
- The code is in `dev/realdata/beta/`; the package is cd342db, unchanged.
- The numbers come from `beta_eval.R` (verdict) and `beta_summary.R` (report-only tables); results are in `results/realdata/beta/evaluation.rds`, `summary.rds`.

**Notation:**
- β is written as a multiple of the default 5e-5, ×1 = 5e-5.
- RMSE ratio = RMSE under the selected β ÷ RMSE under the reference β; the three numbers are distances d0/d1/d2 in order.

## Conclusion: E1–E3 failed; E4 passed by the criteria set in advance, but it is not an estimate usable in the package

> **2026-11-20 follow-up:** following the user's "attempt 2", E4 was decomposed (`decomp/report.md`). E4's gain comes entirely from the missegregation rate used for prediction (changing only that gives +0.45 in the confirmation group); changing only the β used for fitting gives no gain (confirmation group +0.03, interval contains 0), and the high-β fit itself is worse. So the high β chosen by E4 is not evidence that β is higher.

**Verdicts:**
- **E1–E3 failed.** The 6 versions (E1 marginal likelihood, E2 data term, E3 held-out prediction, each in a single and a pooled version) all fail N, A1, A2 on the development set, and did not enter the verdict set.
- **E4 passed.** E4 chooses β by the prediction of later passages of the other trajectories of the same cell line (leave-one-out). By the plan, it looks only at S1:
  - in the development group (SA609, SA039), the per-cell log score is +0.18 higher than at 5e-5;
  - in the confirmation group (SA906, SA532, SA535, SA1035) it is +0.44 higher.
  - Both times the criterion is met.

**Can β be determined from the counts of one dataset:** no.
- In the counts of two time points the information on β is weak. Between ×1/4 and ×4, the criteria differ by only one or two log-likelihood units, and the choice is mainly determined by noise.
- Moreover, the cost of a wrong choice is asymmetric: the loss from β too high is much larger than from β too low. On e1 with truth ×1:
  - fitting with ×4, the RMSE is 3.1–3.6 times that with the true value;
  - fitting with ×1/4, it is only 1.2–1.5 times.
- So choosing β from the counts makes accuracy clearly worse near the truth. On e1, after choosing β the RMSE is 1.35–2.15 times that with the true value fixed (median over grid cells).

**What E4 shows and what it does not show:**
- **It shows:** on Salehi SA532, SA535, SA609, SA1035, fitting and predicting with a β much higher than 5e-5 (×8–×64) makes the prediction of later passages clearly better. SA906 is best near ×1; SA039 has only 2 trajectories, and the differences between grid levels are small.
- **It cannot show "β itself":**
  - the missegregation rate used for prediction also takes the same β, and the higher β, the more dispersed the predicted distribution. This favours the log score and is unrelated to whether the fitness estimates improve.
  - At minobs 10, the landscapes for all β contain the same set of karyotypes; they differ only in fitness values and in the predicted missegregation rate. The two cannot be separated in this design.
- **The choice hits the upper end of the grid:** 3 of the 4 cell lines in the confirmation group chose ×64; the real optimum may lie outside the grid.
- **E4 is not a method usable on a single dataset:** it needs later passages of the same cell line, so it cannot be put into `fit_sparse_joint` as an automatic choice.

**Recommendations:**
- Keep β fixed at 5e-5 in the package, with no automatic choice. This agrees with the conclusion of the beta auto study of the alfak2 period.
- The E4 result suggests that for some cell lines 5e-5 may be too low. This points in the same direction as two observations of the C–E phase: first, C2 refitted at β_sim 1.5e-4 and 5e-4 had clearly lower errors (at 1.5e-3 it got worse instead, see Section 4); second, the excess in O comes from rare karyotypes that appear.
- To confirm this, "the β used for fitting" and "the missegregation rate used for prediction" need to be tested separately (see Section 8).
- Whether it goes into the package and whether to follow up is for you to decide.

## 1. Criteria on the development set

**Criteria (details in `PLAN.md`):**
- **N:** ≥ 18 of the 24 grid cells of e1 are ≤ 1.03 at all of d0–d2, and the median of the grid-cell medians ≤ 1.02; the medians of e2, e3 ≤ 1.03.
- **A1:** ≥ 70% of choices for single versions, and ≥ 80% of groups for pooled versions, fall within a factor of 2 of the truth.
- **A2:** for each true value, the median RMSE ratio against fixed 5e-5: < 1 when the truth is ≥ ×4, ≤ 1.02 at C2 ×3.
- **S1:** the per-cell log score on Salehi later passages is not worse than at 5e-5.
- **S2:** in ≥ 70% of trajectory pairs of Salehi parallel lineages, the chosen β differ by no more than a factor of 2.

| Estimate | N: e1 grid cells passed (/24) | N: e1 grid-cell median d0/d1/d2 | N: e2 / e3 median (d1) | A1 | A2 true values not met | S1 | S2 |
|---|---|---|---|---|---|---|---|
| E1 single | 0 | 1.90 / 1.83 / 1.56 | 1.35 / 1.00 | 0.46 | e4 ×4 (1.10/1.07/1.00) | +0.14 yes | 0.37 |
| E1 pooled | 7 | 1.35 / 1.44 / 1.44 | 1.25 / 0.63 | 0.44 | e4 ×4 (1.08/1.03/1.00) | +0.27 yes | — |
| E2 single | 0 | 1.48 / 1.65 / 1.64 | 1.42 / 0.63 | 0.18 | C2 ×3, e4 ×4, e4 ×16, C2 ×30 | −0.12 no | 0.57 |
| E2 pooled | 5 | 1.54 / 1.51 / 1.52 | 1.31 / 0.58 | 0.09 | C2 ×3, e4 ×4, e4 ×16 | −0.26 no | — |
| E3 single | 0 | 2.11 / 2.15 / 2.07 | 1.63 / 0.65 | 0.38 | e4 ×4, e4 ×16, C2 ×30 | +0.14 yes | 0.35 |
| E3 pooled | 4 | 1.49 / 1.48 / 1.55 | 1.26 / 0.63 | 0.33 | e4 ×4, e4 ×16 | −0.0003 no | — |
| E4 (Salehi only, pooled, leave-one-out) | — | — | — | — | — | development +0.18 yes; confirmation +0.44 yes | — |

**How to read this table:**
- **N fails for all, and by a wide margin.** With truth ×1, the chosen β makes the median RMSE of e1 1.35–2.15 times larger.
- **The e3 median < 1 is an artefact.** e3 is a neutral culture in which all karyotypes have the same true fitness; the lower β, the closer the prediction is to a constant, and the RMSE is naturally smaller (see "the best β" in Section 3). This does not mean the choice is right.
- **S1 (computed only on the development group SA609, SA039):** both versions of E1 and single E3 are not worse than 5e-5, but they already fail N, A1, A2.
- **A1 is far below 70%/80% for all.**
- **The A2 situation:**
  - both versions of E1 meet A2 on the three true values of C2 and on e4 ×16, missing only e4 ×4;
  - pooled E3 meets it on C2, but not on e4 ×4, ×16;
  - on C2 ×10, the RMSE ratios of these three versions are 0.55–0.86, better than fixed 5e-5;
  - but on e4 ×4, no version is better than fixed 5e-5.

## 2. What each estimate chose

Single versions, all landscapes (development and verdict together, only to see the distribution). Rows are set and true value, columns are the chosen multiple, numbers are counts of datasets.

**E1 marginal likelihood:** the choice rises with the true value, but is widely spread.

| Set (truth) | ×1/4 | ×1/2 | ×1 | ×2 | ×4 | ×8 | ×16 | ×32 | ×64 | n |
|---|---|---|---|---|---|---|---|---|---|---|
| e1 (×1) | 48 | 51 | 42 | 67 | 50 | 18 | 7 | 5 | 0 | 288 |
| e2 (×1) | 41 | 38 | 22 | 27 | 19 | 8 | 0 | 0 | 0 | 155 |
| e3 (×1) | 102 | 37 | 36 | 34 | 47 | 23 | 9 | 0 | 0 | 288 |
| e4 (×4) | 16 | 15 | 20 | 25 | 27 | 22 | 11 | 7 | 0 | 143 |
| e4 (×16) | 16 | 6 | 9 | 13 | 22 | 15 | 26 | 10 | 6 | 123 |
| C2 (×3) | 3 | 8 | 23 | 30 | 24 | 7 | 1 | 0 | 0 | 96 |
| C2 (×10) | 2 | 1 | 7 | 14 | 26 | 25 | 10 | 8 | 3 | 96 |
| C2 (×30) | 5 | 0 | 6 | 7 | 7 | 17 | 17 | 10 | 9 | 78 |

**E2 data term:** whatever the true value, it leans to the lowest β. On e1 138/288 and on C2 ×30 37/78 chose ×1/4. The reason is that the data term at the posterior mode gets smaller as β decreases: when β is low, the count of each non-major karyotype is explained mainly by its own fitness, and the fit is tighter.

**E3 held-out prediction:** similar to E1, it also rises with the true value, but is more spread than E1. On C2 ×30, every level from ×1/4 to ×64 has choices.

**How flat each criterion is over the grid:** the table below is the dataset median of "criterion value − its minimum over the grid", in units of negative log-likelihood.

| Criterion | Set | ×1/4 | ×1/2 | ×1 | ×2 | ×4 | ×8 | ×16 |
|---|---|---|---|---|---|---|---|---|
| E1 | e1 (×1) | 3.2 | 2.1 | 2.4 | 2.7 | 4.1 | 10.6 | 27.1 |
| E1 | C2 (×10) | 8.3 | 7.0 | 5.5 | 2.5 | 1.4 | 2.0 | 5.7 |
| E3 | e1 (×1) | 2.3 | 1.9 | 1.9 | 2.5 | 3.3 | 5.5 | 16.2 |
| E2 | e1 (×1) | 0.2 | 3.2 | 4.9 | 7.0 | 10.6 | 16.2 | 25.9 |

- On e1, across the four levels ×1/4 to ×2, the medians of E1 differ by less than 1.2 units, and those of E3 by less than 0.7 units.
- Only when far off (above ×8 on e1, ×1/4 to ×1 on C2 ×10) does the criterion get clearly worse.
- That is, the data can only rule out very wrong values of β, and cannot pin down a value within a factor of 2.

## 3. How much a wrong β costs: accuracy curves

The table below is the dataset median of RMSE(fit β) ÷ RMSE(true β), distance d1. The true values of C2 are not on the grid; the nearest level (×3 → ×4, ×10 → ×8, ×30 → ×32) is used as reference.

| Set (truth) | ×1/4 | ×1/2 | ×1 | ×2 | ×4 | ×8 | ×16 | ×32 | ×64 |
|---|---|---|---|---|---|---|---|---|---|
| e1 (×1) | 1.46 | 1.18 | 1.00 | 1.68 | 3.60 | 6.30 | 9.26 | 14.3 | 14.9 |
| e2 (×1) | 1.23 | 1.11 | 1.00 | 1.08 | 1.66 | 3.26 | 4.96 | 5.87 | 5.63 |
| e3 (×1) | 0.42 | 0.61 | 1.00 | 1.71 | 2.94 | 5.13 | 8.59 | 11.9 | 12.1 |
| e4 (×4) | 1.63 | 1.59 | 1.35 | 1.10 | 1.00 | 1.72 | 3.24 | 4.34 | 4.08 |
| e4 (×16) | 1.16 | 1.14 | 0.99 | 0.96 | 1.04 | 0.98 | 1.00 | 1.14 | 1.13 |
| C2 (×3) | 2.00 | 1.83 | 1.43 | 1.08 | 1.00 | 2.12 | 4.85 | 8.44 | 11.2 |
| C2 (×10) | 1.97 | 1.92 | 1.93 | 1.74 | 1.30 | 1.00 | 1.69 | 3.42 | 5.19 |
| C2 (×30) | 0.91 | 0.68 | 0.58 | 0.65 | 0.82 | 0.68 | 0.76 | 1.00 | 0.98 |

- **The cost is asymmetric:** in each row, the side above the truth rises much faster than the side below. On e1, ×4 is 3.6 times, ×1/4 only 1.46 times.
- **The curves of e4 ×16 and C2 ×30 are very flat, and even best at ×1–×2.** When the true missegregation rate is very high, the β used for fitting has little effect on accuracy, and the true value is not the best.
- **"The best β" is not always the true value.** Take, for each dataset, the level with the lowest d1 RMSE:
  - on e1 only 149/288 (52%) fall at ×1;
  - on e3, 201/288 fall at ×1/4;
  - on C2 ×30, every level from ×1/4 to ×64 occurs.
  - Even if the true value could be chosen perfectly, accuracy would not be the best on part of the data.

## 4. C2: fitting with the true β_sim (report only)

RMSE ratio = fit with β_sim ÷ fit with ×1, dataset median, d0/d1/d2. "Grid best" is the level with the lowest RMSE on the grid for each dataset, i.e. the upper bound of picking after the fact.

| Setting | n | With β_sim | Grid best |
|---|---|---|---|
| β_sim 1.5e-4, bottleneck 300 | 48 | 0.85 / 0.76 / 0.60 | 0.82 / 0.70 / 0.64 |
| β_sim 1.5e-4, bottleneck 3000 | 48 | 0.73 / 0.64 / 0.63 | 0.69 / 0.57 / 0.52 |
| β_sim 5e-4, bottleneck 300 | 48 | 0.64 / 0.67 / 0.56 | 0.62 / 0.48 / 0.46 |
| β_sim 5e-4, bottleneck 3000 | 48 | 0.65 / 0.44 / 0.44 | 0.54 / 0.39 / 0.36 |
| β_sim 1.5e-3, bottleneck 300 | 43 | 0.94 / 1.37 / 1.13 | 0.61 / 0.34 / 0.31 |
| β_sim 1.5e-3, bottleneck 3000 | 35 | 1.35 / 2.08 / 2.09 | 0.98 / 0.35 / 0.34 |

- At β_sim 1.5e-4 and 5e-4, fitting with the true β is 15–56% better than with 5e-5, consistent with the results of the C–E phase.
- At β_sim 1.5e-3 (×30), using the true β is worse instead; at bottleneck 3000 the RMSE at d1, d2 is about 2 times. The grid does have a much better level, but it differs between datasets, and no criterion can pick it out (Section 2).
- **This means that when the missegregation rate is very high, the β in the model and the missegregation rate in the simulator are no longer the same quantity.** One possibility: the model domain includes only the majors and their one-step neighbours, whereas when the missegregation rate is very high, each division often changes several chromosomes at once. This has not been verified separately.

## 5. β chosen on real data (report only)

Since all the estimates above are unreliable on simulations, the choices below cannot be taken as estimates of β; they only show what would be chosen on real data.

**Salehi (single versions: median of the distribution; pooled: one value per cell line):**

| Cell line | n | E1 single median | E1 pooled | E2 single median | E2 pooled | E3 single median | E3 pooled | E4 pooled |
|---|---|---|---|---|---|---|---|---|
| SA039 | 6 | ×1 | ×1 | ×1/4 | ×1/4 | ×0.6 | ×1/4 | ×2, ×8 (2 trajectories) |
| SA609 | 111 | ×16 | ×8 | ×1/4 | ×1/4 | ×8 | ×1 | ×8 (47), ×64 (2) |
| SA1035 | 20 | ×2 | ×4 | ×1 | ×1 | ×2 | ×2 | ×64 |
| SA532 | 28 | ×6 | ×8 | ×1/4 | ×1/4 | ×3 | ×4 | ×64 |
| SA535 | 39 | ×8 | ×8 | ×1/2 | ×1 | ×4 | ×1 | ×64 |
| SA906 | 49 | ×1/2 | ×1/2 | ×1/4 | ×1/4 | ×1/2 | ×1 | ×1/2 (24), ×1 (12) |

- E4 is leave-one-out: the choice for each trajectory is determined by the other trajectories of the same cell line; the numbers in brackets are counts of trajectories. E4 includes the 136 trajectories that were fitted successfully and have later passages ≤ 60 days.
- **Pooled E1 and E4 agree in direction at the cell-line level:** both put SA906 lowest (×1/2), and SA609, SA532, SA535 high (×8 or above). But E1 is unreliable on simulations, so this is only an observation.
- The 24 datasets of SA609 are model-domain errors at every β and are not in the table (as in B).
- Within the same cell line, the choices of single versions spread from ×1/4 to ×64 (e.g. SA609 E1), and the agreement rate of parallel lineages is only 0.35–0.57 (S2).

**NUMBAT (diploid added back, division rate 0.141/day; pooled versions):**

| Cancer type | n | E1 | E2 | E3 |
|---|---|---|---|---|
| GBM | 20 | ×1/4 | ×1/4 | ×1 |
| HNSC | 4 | ×1 | ×1/4 | ×4 |
| Melanoma | 4 | ×1/4 | ×1/4 | ×1/4 |
| Ovarian | 13 | ×1 | ×1/2 | ×8 |

In addition, 1 NUMBAT input is a model-domain error at every β.

## 6. E4: Salehi later passages (S1)

**Method (as in `PLAN.md`):**
- For the alfak results of each β (minobs 10), predict later passages with `forecast_score.R` as is. The missegregation rate for prediction takes the same β.
- Scoring uses the per-cell log score on the common support, `ll_fc_s_e2`, looking only at later passages ≤ 60 days.
- There are 136 trajectories and 301 later passages in total. A further 18 SA609 trajectories (dataset names starting with SA000 or SA609) are model-domain errors at every β and were not fitted; and 1 SA039 trajectory (SA039U_X30_l_2) has all its later passages after 60 days.

**For each cell line and each β, how much higher than ×1:** the table below is the difference of `ll_fc_s_e2` from ×1 on later passages, mean first, median in brackets, per cell.

| Cell line | Trajectories / later passages | ×1/4 | ×1/2 | ×2 | ×4 | ×8 | ×16 | ×32 | ×64 |
|---|---|---|---|---|---|---|---|---|---|
| SA039 | 2 / 2 | −0.02 | −0.03 | +0.02 | +0.03 | +0.02 | −0.01 | −0.09 | −0.23 |
| SA609 | 49 / 134 | −0.27 (−0.28) | −0.12 (−0.14) | +0.15 (+0.16) | +0.26 (+0.28) | +0.27 (+0.37) | +0.18 (+0.37) | +0.09 (+0.39) | +0.25 (+0.48) |
| SA1035 | 12 / 20 | −0.07 (−0.11) | +0.02 (−0.06) | +0.05 (+0.05) | −0.01 (+0.08) | −0.06 (+0.09) | −0.11 (+0.16) | −0.18 (+0.20) | +0.50 (+0.31) |
| SA532 | 21 / 46 | −0.21 (−0.20) | −0.23 (−0.14) | +0.19 (+0.13) | +0.31 (+0.19) | +0.31 (+0.22) | +0.44 (+0.27) | +0.64 (+0.37) | +0.82 (+0.55) |
| SA535 | 16 / 49 | −0.21 (−0.27) | −0.16 (−0.16) | +0.18 (+0.16) | +0.24 (+0.26) | +0.35 (+0.36) | +0.39 (+0.35) | +0.41 (+0.14) | +0.53 (+0.30) |
| SA906 | 36 / 50 | −0.02 (−0.02) | +0.00 (−0.00) | −0.03 (−0.00) | −0.12 (−0.04) | −0.14 (−0.08) | −0.15 (−0.05) | −0.07 (−0.01) | −0.16 (−0.06) |

- **SA609, SA532, SA535:** as β goes up from ×1/4, the prediction improves all the way, and is still rising above ×8. At ×8, 76–87% of later passages are better than ×1.
- **The non-monotone mean of SA1035 is caused by one trajectory.** SA1035T_X5_l_2 is worse than ×1 by 2.1 and 3.3 at ×16 and ×32, and better by 0.47 at ×64. The median rises steadily from ×1/2 to ×64.
- **SA906 is best near ×1, and the higher β, the worse.** SA039 has only 2 trajectories, and nothing can be seen.

**E4's choice and S1:**
- **Development group:** of the 49 SA609 trajectories, 47 chose ×8 and 2 chose ×64; the two SA039 trajectories chose ×2 and ×8. S1 = +0.18, met.
- **Confirmation group:** SA1035, SA532, SA535 all chose ×64; SA906 chose ×1/2 (24 trajectories) or ×1 (12 trajectories). S1 = +0.44, met.

**Why this cannot be read directly as "these cell lines have a higher β":**
- **The predicted missegregation rate changes together with β.** The higher β, the more dispersed the predicted distribution, and the less likely a very low log score when new karyotypes appear. This effect is unrelated to how well the landscape is estimated.
- **The fitness values also change with β.** In this design the two change at the same time and cannot be separated.
- **The grid is not wide enough.** 3 cell lines of the confirmation group chose the upper end ×64; the optimum may be higher.

## 7. Problems during the run and deviations from the plan

- **The first prediction of later passages was all "no fit".**
  - Cause: `benchmark/paper/R/fit.R` also defines `beta_label()`, producing labels like x1_4 … x64. `beta_salehi.R` loads it after `R/beta.R`, which overrode my function. So the directory names of the alfak-format results are x labels, while the prediction task table used b labels.
  - The fit itself is not affected: `multiple` in the records is correct.
  - Fix: the task table was changed to x labels (the predicted missegregation rate unchanged), the old output moved to `horizon_nofit_blabels/`, and resubmitted (22152606).
- **Errors in how `beta_eval.R` reads data, corrected before looking at results:**
  - e1–e4 reuse the same dataset names; merging by name would have mixed the four sets, so identification was changed to kind/set/id;
  - the output for later passages is `.Rds`, with the contents in `$rows`, and the labels are x labels.
  - None of these touch the criteria themselves.
- **Failed fits (allowed by the stopping rules):**
  - 414 are model-domain errors, with the same dataset failing at every β: e4 22 datasets, C2 β_sim 1.5e-3 18, Salehi SA609 24, NUMBAT 1.
  - 46 are "no candidate history" of e2 at ×32, ×64, involving 37 datasets, 2.8% of the simulation tasks, below 5%.
  - The evaluation uses only datasets fitted successfully at all 9 β (1267/1344), so e2 lacks these 37 (155/192 left). Even with them added back, N already passes 0/24 grid cells on e1, and the verdict would not change.
- **NODE_FAIL:** the simulation arrays had several NODE_FAIL, which Slurm requeued automatically and which then completed; nothing was submitted twice.
- **File names:** the β labels in the file names of the simulation results contain spaces (e.g. `b 1.rds`), and the Salehi results have x labels. The evaluation always reads by `multiple` in the records, so it is not affected.

## 8. What could be done next (for you to decide)

1. **Keep β fixed at 5e-5 in the package (recommended).** The package is not changed. None of E1–E3 can determine β from a single dataset.
2. **(Done, see `decomp/report.md`: verdict "mainly a more dispersed prediction"; the fitted β itself contributes nothing.) Decompose the effect of E4 (if you want to know whether β is really higher for some Salehi cell lines):** do a 2 × 2 on Salehi with the existing fit and prediction scripts, which only needs the prediction recomputed:
   - fit β ∈ {×1, ×chosen value};
   - prediction missegregation rate ∈ {×1, ×chosen value}.
   - If changing only the prediction missegregation rate gives most of the gain, it is an effect of prediction dispersion, not an estimate of β.
   - At the same time, extend the grid upward to ×128, ×256.
3. **Use an externally measured missegregation rate.** If a cell line has a direct measurement such as live-cell imaging, it can be passed in as the β for those data (`fit_sparse_joint(beta = ...)` already supports this). But Section 4 shows that when the missegregation rate is very high, using the true value directly is not necessarily better.
4. **Multi-time-point data.** This study uses only the first and last time points. Whether more time points can determine β has not been tested.
