# Why the forecast is too concentrated: results

**Basis:**
- The plan and the decision rules set in advance are in `PLAN.md`, fixed 2026-11-20 08:52 EST, with no corrections.
- Numbers come from `disp_eval.R`; results are in `results/realdata/beta/disp/evaluation.rds`. There are also three report-only tables (section 4), computed from the same batch of outputs.
- **Subjects:** 137 trajectories, 301 later passages ≤ 60 days, fit fixed at 5e-5.
- **Check:** base at F = ×1, ×16 is identical to the scores of the split study (maximum difference 0).
- **What "landscape" means:** the alfak-format karyotype list into which the KARST fit is converted by `as_alfak()` (minobs 10): the frequent karyotypes and their neighbours within two steps. The forecast uses alfakR's ABM, which only knows this list. KARST's `predict()` and alfak's kriging can both give fitness for karyotypes outside the list; the limit here comes from the forecast step (`predict_evo()` only accepts a fixed list), consistent with what the paper benchmark does.

## Conclusion

**Neither intervention passed the decision rules set in advance.**

| Intervention | Confirmation group gain at F = ×1 [95% interval] | Development group | Verdict |
|---|---|---|---|
| I1 parent frequencies plus the Good–Turing unseen part | −0.02 [−0.16, 0.09] | +0.17 [−0.07, 0.43] | Fail (threshold +0.23) |
| I2 ABM bottleneck 1e3 | −0.48 [−0.60, −0.36] | −0.53 | Fail, worse |
| I2 ABM bottleneck 1e4 | −0.10 [−0.14, −0.07] | −0.21 | Fail, worse |

For reference, only raising the forecast missegregation rate to ×16 (base F = ×16): confirmation group +0.43 [0.36, 0.51], development group +0.78.

**Causes given by the diagnostics (in order of importance):**

1. **Later passages contain many karyotypes far from the parent that fall outside the landscape (H3, H4).** In the confirmation group 19% and in the development group 41% of cells are ≥ 2 steps from all parent sample karyotypes. Almost all of these cells are outside the landscape; the forecast cannot give them any probability, and scoring simply skips them.
   - SA535 and SA609 are the worst: 51–100% of the later-passage cells are outside the landscape.
   - That is, the log score only evaluates the part of the cells that can be evaluated.
2. **Inside the landscape, the gain comes mainly from new karyotypes one step from the parent (class C).**
   - Of the confirmation group gain from F ×1 → ×16, class C contributes +0.28, class A (already in the parent) −0.10.
   - The 5e-5 forecast puts only 3.5–4.3% of the probability at distance one, while 16–34% of the cells are actually there. Converted per cell line, the generation rate of one-step new karyotypes is 2.4–11 times what 5e-5 gives.
   - At distance two the gap is larger, converting to 4–17 times. This indicates that single changes spanning several steps are more common than a single missegregation rate gives.
3. **New karyotypes increase over time; they are not only ones the parent sample missed (H1 explains only part).**
   - Within the same trajectory, later passages have more new karyotypes than earlier ones: this holds for 66% of trajectories, and for the part ≥ 2 steps from the parent, for 72% of trajectories.
   - The Good–Turing unseen part of the parent sample is large (median 0.37, 34% of trajectories at the cap 0.5). So rare karyotypes missed in the parent do exist, but they are not enough to explain everything.
   - Why I1 failed: it puts the unseen part on karyotypes seen in earlier passages (120/136 trajectories), and this class of karyotypes (class B) makes up only 3–4% of the later-passage cells. It was put in the wrong place; this does not show that H1 is absent.
4. **The frequencies of karyotypes already in the parent change far more than sampling noise (a sign of H2).**
   - The Pearson dispersion φ has a median of 5.6 (15 days) to 107 (60 days), giving a rough effective bottleneck per passage of only 25–130 cells.
   - But error in the selection estimate also raises φ, so this is only a lower bound.
   - Adding a small bottleneck in the ABM (I2) made it worse. Possible reason (not verified separately): there are only 5 replicates, the bottleneck makes each replicate lose rare karyotypes, and after averaging the tail is thinner. So I2 only shows that adding drift this way cannot correct the tail; it does not rule out drift affecting the frequencies of existing karyotypes.

**In one sentence:** the forecast at 5e-5 is too concentrated mainly because in the real later passages, new karyotypes arise several times faster than 5e-5, and often span two or more steps at once. Most of those that span far fall outside the landscape; the forecast cannot cover them at all, and scoring cannot see them. Raising the forecast missegregation rate only fills in the one-step and two-step new karyotypes inside the landscape.

## 1. D1: which class of cells the gain from F ×1 → ×16 comes from

Only cells inside the landscape. Contribution = fraction of cells in the class × change in log score per cell, pooled over cells. The mean over later passages (+0.43) and the pooled value here (+0.28) use different weights.

| Group | Class | Fraction of cells | Log score per cell, F ×1 | F ×16 | Contribution to gain |
|---|---|---|---|---|---|
| Confirmation | A already in parent | 0.82 | −2.86 | −2.98 | −0.10 |
| Confirmation | B seen earlier | 0.04 | −10.99 | −8.76 | +0.08 |
| Confirmation | C new, one step | 0.13 | −10.78 | −8.59 | +0.28 |
| Confirmation | D new, ≥ 2 steps | 0.02 | −13.16 | −11.97 | +0.02 |
| Development | A | 0.63 | −3.63 | −3.86 | −0.15 |
| Development | B | 0.03 | −10.78 | −8.54 | +0.08 |
| Development | C | 0.30 | −10.93 | −8.81 | +0.64 |
| Development | D | 0.03 | −12.90 | −11.94 | +0.03 |

- Class B, C and D cells have only about −11 to −13 per cell at F ×1, close to the uniform lower bound (about −13), i.e. the forecast gives them almost no probability.
- F ×16 raises them by about 2 units, at the cost of a small drop for class A.

## 2. D2: new karyotypes and parent sampling

"Fraction of new karyotypes" is the fraction of cells in the later passage with karyotypes absent from the parent sample, including cells outside the landscape.

| Group | Days | Number of later passages | Fraction of new karyotypes | Parent Good–Turing unseen part | Fraction outside landscape |
|---|---|---|---|---|---|
| Confirmation | 15 | 55 | 0.46 | 0.36 | 0.36 |
| Confirmation | 30 | 55 | 0.39 | 0.27 | 0.28 |
| Confirmation | 45 | 22 | 0.57 | 0.31 | 0.40 |
| Confirmation | 60 | 33 | 0.29 | 0.18 | 0.15 |
| Development | 15 | 65 | 0.64 | 0.48 | 0.51 |
| Development | 30 | 38 | 0.74 | 0.49 | 0.66 |
| Development | 45 | 21 | 0.63 | 0.43 | 0.51 |
| Development | 60 | 12 | 0.84 | 0.50 | 0.69 |

- In every cell of the table, the fraction of new karyotypes is higher than the Good–Turing unseen part.
- The groups at each day are made up of different cell lines, so comparing time trends across rows is not reliable. See section 4 for paired comparisons within the same trajectory.

**All cells (including outside the landscape) by class:**
- Confirmation group: A 0.64, B 0.03, C 0.14, D 0.19;
- Development group: A 0.33, B 0.03, C 0.23, D 0.41.

## 3. D3: distance from the parent

The distance is the copy-number L1 distance to the nearest parent sample karyotype.

| Group | Distance | Observed: all cells | Observed: inside landscape | Forecast F ×1 | Forecast F ×16 |
|---|---|---|---|---|---|
| Confirmation | 0 | 0.64 | 0.82 | 0.957 | 0.734 |
| Confirmation | 1 | 0.17 | 0.16 | 0.043 | 0.241 |
| Confirmation | 2 | 0.08 | 0.021 | 0.0006 | 0.025 |
| Confirmation | ≥ 3 | 0.11 | 0.0008 | ≈ 0 | ≈ 0 |
| Development | 0 | 0.33 | 0.63 | 0.965 | 0.765 |
| Development | 1 | 0.26 | 0.34 | 0.035 | 0.223 |
| Development | 2 | 0.17 | 0.035 | 0.0001 | 0.011 |
| Development | ≥ 3 | 0.23 | 0.0001 | 0 | ≈ 0 |

- **Inside the landscape,** the one-step and two-step probabilities of F ×16 are close to the observed ones (confirmation group 0.24 vs 0.16, 0.025 vs 0.021); this is why it scores higher.
- **The part outside the landscape:** cells at distance ≥ 2 are 19% (confirmation) and 41% (development) of all cells, but there are almost none inside the landscape. The landscape (`as_alfak`, minobs 10) contains only the frequent karyotypes (≥ 10 cells summed over the two input time points) and their neighbours within two steps. Most karyotypes in the parent sample have only one cell and are not frequent karyotypes, so karyotypes two steps from them are often further from the frequent karyotypes and are not in the list.

## 4. Report-only supplementary tables

**Fraction of new karyotypes by cell line and day** (new = absent from the parent sample; far = new and ≥ 2 steps; outside = outside the landscape):

| Cell line | Days | New | Far | Outside |
|---|---|---|---|---|
| SA1035 | 15 / 30 / 45 | 0.20 / 0.19 / 0.20 | 0.06 / 0.04 / 0.04 | 0.03 / 0.03 / 0.04 |
| SA532 | 15 / 30 / 45 / 60 | 0.34 / 0.41 / 0.52 / 0.51 | 0.13 / 0.20 / 0.21 / 0.21 | 0.16 / 0.18 / 0.22 / 0.25 |
| SA535 | 15 / 30 / 45 / 60 | 0.65 / 0.80 / 0.95 / 1.00 | 0.39 / 0.58 / 0.85 / 0.99 | 0.65 / 0.76 / 0.93 / 1.00 |
| SA609 | 15 / 30 / 45 / 60 | 0.64 / 0.74 / 0.70 / 0.83 | 0.37 / 0.48 / 0.44 / 0.54 | 0.51 / 0.66 / 0.57 / 0.69 |
| SA906 | 25 / 35 / 50 | 0.15 / 0.25 / 0.17 | 0.05 / 0.07 / 0.06 | 0.07 / 0.07 / 0.07 |
| SA039 | 45 | 0.02 | 0.00 | 0.01 |

**Within the same trajectory, the earliest and latest later passages** (76 trajectories with more than two later passages):

| Cell line | Trajectories | New: early → late | Far: early → late |
|---|---|---|---|
| SA535 | 9 | 0.61 → 0.85 | 0.34 → 0.66 |
| SA609 | 32 | 0.65 → 0.73 | 0.38 → 0.48 |
| SA532 | 15 | 0.40 → 0.45 | 0.16 → 0.21 |
| SA1035 | 6 | 0.16 → 0.22 | 0.04 → 0.05 |
| SA906 | 14 | 0.14 → 0.13 | 0.04 → 0.05 |

- In 66% of trajectories the fraction of new karyotypes rises, and in 72% of trajectories the "far" part rises.
- At 60 days almost all cells (0.99) of SA535 are new karyotypes ≥ 2 steps from the parent. It looks like a clonal replacement, and the new clone is not in the landscape.

**Rate multiples converted from the one-step and two-step fractions inside the landscape:**
- One step: observed fraction ÷ forecast probability at F ×1.
- Two steps: the square root of the same ratio, because the two-step probability is roughly proportional to the square of the rate.

| Cell line | One step | Two steps |
|---|---|---|
| SA039 | 2.4 | — |
| SA532 | 4.2 | 4.2 |
| SA535 | 4.3 | 12.6 |
| SA1035 | 7.5 | 16.5 |
| SA609 | 9.7 | 17.0 |
| SA906 | 10.8 | 15.2 |

- Of the 5 cell lines with two-step data, all except SA532 give a clearly higher multiple from two steps than from one step.
- If the missegregation rate were simply higher, the two should be similar. Two steps being higher indicates that single changes spanning several steps (or selection for them) are more frequent than a single-rate model gives.

**Where I1 puts the unseen part:** 120 trajectories put it on karyotypes seen in earlier passages, 16 share it equally among one-step neighbours.

## 5. Assessment of the four candidate causes

| Cause | Evidence | Assessment |
|---|---|---|
| H1 incomplete sampling of the parent | The unseen part is large (median 0.37), but the fraction of new karyotypes is larger and increases over time; I1 fails, but it put the mass in the wrong place | Present, but not the main cause; the intervention in this study cannot test it |
| H2 passage bottleneck | φ of existing karyotypes reaches 5–107, rough N_e 25–130; but adding a bottleneck to the ABM is worse | The frequencies of existing karyotypes may be affected by drift; it is not the cause of the too-thin tail |
| H3 multi-chromosome changes | The rate converted from two steps is higher than from one step; 19–41% of later-passage cells are ≥ 2 steps | Supported |
| H4 outside the landscape | By group and day, 15–69% of later-passage cells are outside the landscape (SA535 up to 100%), skipped by scoring | Supported, and it is the largest coverage gap |

## 6. What could be done next (for you to decide)

1. **The fit β stays at 5e-5.** This study does not change that.
2. **When reporting later-passage forecasts, also report the fraction of cells outside the landscape.** The current log score only evaluates cells inside the landscape, and on SA535 and SA609 misses more than half.
3. **Widen the karyotype range in the forecast.**
   - KARST itself is not limited by this list: `predict()` gives fitness for any karyotype. Inside the model domain it uses the posterior; observed karyotypes outside the domain are placed at the level of the non-major observed karyotypes inside the domain; the remaining karyotypes are given by the trend, with extra uncertainty added.
   - The karyotype list passed to the ABM could be extended outwards, with the fitness of the extended part from KARST's `predict()`, and then see whether cells outside the list can get a probability.
   - The extension can only be set by distance around the parent, not by the karyotypes that actually appear in the later passages; otherwise future data would be used.
   - The scoring support also has to change accordingly, to be comparable with the current results.
   - alfak itself can also extrapolate: it fits fitness with kriging (`fields::Krig`), and the landscape is just this model evaluated on a set of karyotypes. With `landscape_data_output = TRUE`, alfak saves the kriging model (`landscape_data.Rds`), which can predict any karyotype. `as_alfak()` has an option of the same name that does the same with KARST's `predict()`.
   - So the limit is not in the inference of alfak or KARST, but in the forecast step: the ABM of `alfakR::predict_evo()` only accepts a fixed karyotype list and does not extrapolate on the fly.
   - A fair extension should treat both methods the same way: alfak with kriging, KARST with `predict()`.
   - Correction (11-20 10:26): the GRF branch of the ABM in C++ (`grf_centroids`, `grf_lambda`) is the true GRF landscape used for simulation (a sum of sinusoids with several centres as wave sources), not kriging extrapolation, and is unrelated to this.
   - This needs a new plan.
4. **Set the forecast missegregation rate separately.** The split study already showed that around ×16 is better inside the landscape; this study shows that what it fills in is the one-step and two-step new karyotypes. To become a recommendation, it still needs to be tested in advance on the ABM or other data.
5. **Further interventions in the drift direction are not recommended.** H2 can only be tested fairly if the forecast is changed to a mixture distribution over many replicates.
