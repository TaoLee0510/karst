# Pre-fit preprocessing v7: redo the v6 tests on new, aligned simulations

**Fixed:** 2026-11-20, before any v7 data were generated.

**User confirmation:** 2026-11-20 about 21:45, "Confirmed, go ahead". The four items confirmed:
1. Cell numbers use the 5 combinations among the 8 core cases of the own benchmark;
2. The three designs, selection, clonal turnover and neutral, are each run at ×1, ×4 and ×16, every landscape is run at all three missegregation rates, and likewise for C2;
3. C2 uses `dev/abm/abm_c2.cpp`;
4. The verdict numbers carry over from v6, and step 3 is changed so that the confirmation halves of set e and of C2 must each pass; the paper simulations are report only; steps 0–3 are automatic; if they pass, continue with 4 and 5; if not, stop and report.

**The user's requirements:**
- Entirely new data, generated with the current 1.1.3 code.
- Set e keeps its original structure, with wavelength 0.2 added.
- The cell number per time point matches the own benchmark.
- C2 keeps its Salehi-like structure, with wavelength, cell number and missegregation rate aligned with set e.
- The new data are strictly isolated from earlier simulations.

**Method:** exactly the same as v6 (`PLAN_v6.md`): statistics A and B, calibration, the 25% prediction quantile, upward adjustment only, capped at 10 times, and each verdict number. Only the differences are written below.

**Not done:** the KARST package and the repository code are not changed; no commit.

## Data

### Landscapes (shared by set e and C2)

- 20: wavelengths 0.2, 0.4, 0.8, 1.2, 1.6, 4 each.
- Numbered `v7_L01`–`v7_L20`, ordered by wavelength: 0.2 is L01–L04, 0.4 is L05–L08, and so on.
- Seeds 131001–131020 (base 131000).
- Generation: `set.seed(landscape_seed)`, then `stt_generate_centroids(rep(2, 22), 10, wavelength)`, the same as the original practice for set e and C2. With the same seed and wavelength, the landscapes of the two are identical.
- **Split into halves:** the 1st and 2nd of each wavelength go into the development half, the 3rd and 4th into the confirmation half.

### Missegregation rates

×1, ×4, ×16, i.e. 5e-5, 2e-4, 8e-4. Every landscape is run at all three.

### Cell numbers (first / second time point)

Taken from the 8 core cases of `benchmark/sparse_two_timepoint/config.R`, 5 combinations in total: 50 / 50, 200 / 200, 500 / 700, 100 / 3000, 3000 / 100.

### Set e

**truth:**
- One per landscape × missegregation rate.
- `mt_config(…, "series")`, with the missegregation rate set to the corresponding multiple. `stt_generate_truth()` grows from day 0 at that missegregation rate, the same as the original e4.

**Three designs (window functions unchanged):**

| design | window function | settings | windows |
|---|---|---|---|
| e1 selection | `mt_make_window()` | from day 1440, bottleneck 200 / 2000 / 20000, spacing 6 or 60 days | 6 |
| e2 clonal turnover | `sm_make_window()`, set e2 | from day 240 or 480, spacing 120 days, bottleneck 2000 / 20000 | 4 |
| e3 neutral | `sm_make_window()`, set e3 | as e1; from the first time point on, all karyotypes divide at the same rate | 6 |

- Inside the windows, the simulation continues at the truth's own missegregation rate.
- The window seeds keep each design's original formula, plus missegregation-rate index × 10⁶, so that the random numbers of the three missegregation rates differ.

**Sampling:**
- Each window takes 4 or 8 time points, of which only the first and the last are used. Intervals: e1 / e3 18, 42, 180, 420 days; e2 360, 840 days.
- The two time points are sampled with `stt_sample_counts()` separately for each of the 5 cell-number combinations; the seed is determined by the dataset name.
- Number of datasets: (6 + 4 + 6) windows × 2 × 5 × 3 missegregation rates × 20 landscapes = **9600**.

### C2

- truth: one per landscape × missegregation rate × bottleneck (300, 3000), `c2_simulate()` called as is.
  - Starts from a bottleneck-sized number of diploid cells, passaged every 15 days, cap 2e5 between two passages, first evolves for 720 days, copy number 1–8.
  - The ABM is `dev/abm/abm_c2.cpp`.
- Sampling: t₀ and t₀ + 15, 30, 45, 60 days, 5 cell-number combinations, `stt_sample_counts()`, seed determined by the dataset name.
- Number of datasets: 2 × 4 × 5 × 3 × 20 = **2400**.

### Data for reporting only

Not part of the verdict; recomputed within v7, without reading v6 results:
- the paper's 125 missegregation-rate simulations (750 datasets, taken as in v6);
- the 155 Salehi trajectories;
- the 42 NUMBAT datasets.

## Isolation from earlier data

1. **Seeds:**
   - Base 131000, not used anywhere in the repository or dev code.
   - Already used: 63000, 64000, 71000, 82000 (including 82016–82500), 93000–99000, 101000–103000, 111000, 112000.
   - Sampling seeds are derived from the dataset names.
2. **Checks of step 0:**
   - Read the earlier truth files (read only) and confirm that the seeds and centroid matrices of the new landscapes all differ from them.
   - Scope: this session's results directory, the results directory of the other session's karst clone, and the `true_landscape` of the paper simulations.
   - Stop on any overlap.
3. **Directories:**
   - Everything goes into `results/realdata/noise/v7/`.
   - v7 scripts read and write only this directory, plus the raw inputs of the paper, Salehi and NUMBAT; a path that lands anywhere else raises an error.
4. **No reuse:** earlier truths, datasets, scoring nodes, screening records, statistics and fits are never used. The T0 fits are also rerun.
5. **Version record:** each truth file records the karst version (library karst_6ae2a24, 1.1.3), the local repository commit (942182c; since 1.1.3 the simulation code has changed only in comments and in the held-out landscape list), and the md5 of the simulation source files used.
6. **Order:** calibration uses only the development half; the confirmation half is evaluated only once.

## Steps

**Step 0:** generate the truths and datasets, and run the isolation checks.

**Step 1:** compute S_A and S_B with `R/cloud_level.R`, ρ = 1 (Salehi 1, NUMBAT 0.141).

**Step 2:** development half (6000 datasets, both kinds together).
- Of A and B, choose the one with AUC ≥ 0.80 and the higher AUC; on ties, A.
- Calibrate log₂ S = c + d · log₂(multiple) + ε, requiring the 95% interval of d > 0.

**Step 3:** the set e confirmation half (4800) and the C2 confirmation half (1200) **must each pass**:

| criterion | requirement |
|---|---|
| false raise | among true ×1, β_used ≥ 2 times in ≤ 10% |
| useful | among true ≥ ×4, β_used ≥ 2 times in ≥ 50% |
| no overshoot | β_used > 2 × true value in ≤ 5% |

**Report only:**
- intercept and slope calibrated separately for set e and C2, and their difference (bootstrap interval, clustered by landscape);
- cross transfer: calibrate on the set e development half and test on the C2 confirmation half, and the reverse;
- fraction raised by design, cell number and wavelength;
- β_used of the paper simulations, Salehi and NUMBAT.

**Step 4:** fits (HPC), datasets in the confirmation half with β_used > 5e-5, T0 vs T6.
- Library karst_6ae2a24.
- Fit settings for set e as in v4: two time points, birth rate 1, copy number 1–8. C2 as in the C2 settings of v4.
- Scoring nodes `stt_scoring_nodes()`. Truth: set e uses `sm_true_net_rate()` (constant fitness for e3), C2 uses `c2_true_rate()`.
- Criteria as in v6 (0.95 and cluster upper limit < 1; raised ×1 ≤ 1.02, report only if fewer than 10), computed separately for set e and C2.

**Step 5:** Salehi, as in v6: parallel lineages ≥ 60%; prediction (i) gain ≥ 0, lower limit > −0.1.

## Outputs

- Code: `dev/realdata/noise/v7/`.
- Results: HPC `results/realdata/noise/v7/`, local copy in `dev/results/realdata/noise/v7/`.
- Report: `dev/realdata/noise/report_v7.md`.
