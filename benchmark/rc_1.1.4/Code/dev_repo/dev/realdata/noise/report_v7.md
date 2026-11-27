# Pre-fit preprocessing v7: redoing the v6 test on new, aligned simulations

**Basis:**
- Plan `PLAN_v7.md`, fixed 2026-11-20 before any v7 data were generated. The user confirmed at about 21:45: "Confirmed, go ahead".
- The verdict numbers come from `v7/calib_eval.R`, with output stored in `dev/results/realdata/noise/v7/calib_eval.log`.
- Section 3 is an after-the-fact description that **uses only the development half**, from `v7_posthoc.R` → `posthoc.log`. Not a single number of the confirmation half was looked at; it is left for a possible new plan later.

## Conclusion: step 2 did not pass; stopped as the plan says

KARST keeps the default missegregation rate 5e-5. No calibration, no steps 3–5, and no fits of any kind were done.

| Statistic | AUC (development half, ×1 vs ≥ ×4; required ≥ 0.80) |
|---|---|
| A (neighbour cloud of high-frequency karyotypes) | 0.551 |
| B (adjacent cell pairs ÷ identical cell pairs) | 0.667 |

- The development half has 6000 datasets in total (e sets 4800, C2 1200): ×1 2000, ≥ ×4 4000.
- A cannot be computed for 1079 (no high-frequency karyotype), of which 1013 are ≥ ×4.

## 1. Data (step 0)

- **Generation:** jobs 22236745 (e sets, 60 tasks) and 22236747 (C2, 120), all COMPLETED, the longest about 10 minutes.
  - 9600 e-set and 2400 C2 datasets, as in the plan.
  - The results of the trial runs 22236731 and 22236733 are reused in the main jobs (existing files are skipped).
- **Isolation from earlier data (job 22236742): passed.**
  - Seeds: of the 2195 landscape seeds registered in the code, none equals any of 131001–131020.
  - centroid matrices: 2025 earlier files were read and 2022 of their matrices compared; none is the same as a new landscape, and no file was unreadable. Sources of the files:
    - 390 from this session's results directories (excluding v7);
    - 491 from the other session's karst clone, read-only;
    - 618 from the alfak2 backup;
    - 526 paper simulations.
  - **Correction:** in the first check (job 22236735), `find` did not run because the brackets were not quoted, and only the paper simulations were read, so the "isolated" result of that time is void. After the fix, a guard was added: if fewer than 100 earlier truth files are found, the check is judged failed.
- **Within this round (job 22236751):** all 180 truth files used the planned centroids, and the e sets and C2 share the same set of landscapes.
- **Statistics (job 22236749):** 187 tasks, all 12947 records successful: v7 simulations 12000, paper simulations 750, Salehi 155, NUMBAT 42.

## 2. Comparison with v6

| | v6 (old data) | v7 (new data, development half) |
|---|---|---|
| AUC of B | 0.838 | 0.667 |
| AUC of A | 0.525 | 0.551 |

Part of v6's 0.838 came from "missegregation rate tied to the source of the simulation": high multiples were mostly C2, and S_B of C2 is high to begin with. With v7 fully crossing missegregation rate with design and source, this part is gone.

## 3. After-the-fact description (development half, not a verdict)

### 3.1 AUC of B broken down by factor

| Factor | AUC |
|---|---|
| Design: C2 / e1 / e2 / e3 | 0.678 / 0.709 / 0.669 / 0.642 |
| Number of cells: 50/50, 200/200, 500/700, 100/3000, 3000/100 | 0.654, 0.670, 0.666, 0.640, 0.704 |
| Bottleneck: 200 / 300 / 2000 / 3000 / 20000 | 0.68 / 0.68 / 0.66 / 0.69 / 0.68 |
| **Wavelength: 0.2 / 0.4 / 0.8 / 1.2 / 1.6** | **0.783 / 0.949 / 0.953 / 0.714 / 0.339** |

Design, number of cells, bottleneck and interval have little effect; **wavelength has the largest effect**.

### 3.2 Median of log₂ S_B, by wavelength and missegregation rate

| Wavelength | ×1 | ×4 | ×16 |
|---|---|---|---|
| 0.2 | −11.32 | −9.84 | −7.86 |
| 0.4 | −9.74 | −7.78 | −5.91 |
| 0.8 | −6.66 | −5.32 | −4.37 |
| 1.2 | −6.13 | −5.13 | **−5.98** |
| 1.6 | −5.50 | −5.36 | **−6.39** |

- **At ×1, the gap between wavelengths is about 6 log₂ units, while within one wavelength going from ×1 to ×16 moves only 1–3.5 units.**
- At wavelengths 0.4 and 0.8 it is monotone and separates well (AUC 0.95).
- At wavelengths 1.2 and 1.6, ×16 is lower than ×4 instead; at 1.6 it is even lower than ×1, so the AUC is 0.34. I did not look into the reason.
- The paper simulations have only wavelength 0.8; their AUC of 0.968 in v6 agrees with the 0.953 at wavelength 0.8 here. The paper simulations separated well in v6 because there was only one wavelength.

### 3.3 Paired comparison: the effect of the missegregation rate with all other conditions equal

Same landscape, design, window and number of cells, only the missegregation rate changed:

| Comparison | Pairs | Fraction with higher S_B | Median log₂ difference |
|---|---|---|---|
| ×4 vs ×1 | 2000 | 80.2% | +1.17 |
| ×16 vs ×1 | 2000 | 73.5% | +1.66 |

- By design (×16 vs ×1): C2 0.79, e1 0.78, e2 0.74, e3 0.66.
- So the signal of the missegregation rate exists: under equal conditions, a higher missegregation rate mostly gives a larger neighbour cloud. But it is much smaller than the differences between landscapes.
- Interquartile range grouped by design and missegregation rate: about 2.5–4 log₂ units at ×1, 1.2–1.6 at ×16.

### 3.4 After aligning the e sets and C2

Medians at ×1: C2 −6.38, e1 −7.55, e2 −6.34, e3 −6.81. The difference between the two kinds of simulation (about 1 unit) is much smaller than the differences between wavelengths.

**Conjecture from this (not verified):**
- The level of the neighbour cloud is determined mainly by the landscape: roughness determines the fitness disadvantage of neighbours. The missegregation rate comes second.
- This is the "missegregation rate ÷ neighbour disadvantage" confounding mentioned in the plan.
- The landscape roughness of real data is unknown, so whatever kind of simulation is used, S_B cannot be converted into a missegregation rate. On a smooth landscape, even the direction "above or below the default" may be called the wrong way round.

## Next steps (waiting for your decision; I have not done them)

1. **Recommendation: end the line of "estimating the missegregation rate from the data themselves".**
   - v5 (growth of the neighbour cloud), v6, v7 (level of the neighbour cloud) and the earlier β studies (marginal likelihood, held-out prediction, rate of appearance of new karyotypes) all failed.
   - The missegregation rate stays fixed at 5e-5.
   - If an externally measured missegregation rate is available, input it with `beta` or `fit_sparse_joint_rates(beta_range = …)`.
2. **If it is to be tried again, there must first be independent information on landscape roughness.**
   - For example, a control of the same cell line at a known missegregation rate, or an experimentally measured fitness disadvantage of neighbours.
   - Estimating roughness from the fitness differences fitted by KARST itself is circular: the fit itself depends on the missegregation rate. This is why choosing the missegregation rate by marginal likelihood failed before.
3. The confirmation half of v7 (6000 datasets) has not been used. If there is a new plan fixed in advance later, it can be tested on it directly.

## Files

- Code:
  - `dev/realdata/noise/v7/`: `common_v7.R`, `sim_eset.R`, `sim_c2.R`, `isolation_check.R`, `stats_task.R`, `calib_eval.R`;
  - `fit_task.R` was prepared in advance and not used;
  - task tables `v7/tasks/`.
- After-the-fact description: `dev/results/realdata/noise/v7/v7_posthoc.R` → `posthoc.log`, in the results directory, not in the repository.
- Results:
  - HPC `results/realdata/noise/v7/`: `sims/` (truth, datasets), `stats/`, `isolation_*.rds`;
  - local `dev/results/realdata/noise/v7/`: copies of `stats/` and `isolation_*.rds`, plus the locally computed `calibration.rds`, `calib_eval.log`, `posthoc.log`.
- Jobs:

  | Content | Job ID |
  |---|---|
  | Trial runs | 22236731, 22236733 |
  | Isolation checks | 22236735 (void), 22236742, 22236751 |
  | Generation | 22236745, 22236747 |
  | Statistics | 22236749 |

  Library karst_6ae2a24 (1.1.3).

**User decision (2026-11-20 22:46):** confirmed ending this line. The missegregation rate is fixed at 5e-5 and is no longer estimated from the data themselves.
