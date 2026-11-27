# Per-chromosome marginal pseudo-likelihood: results (stopped at step 0, as planned)

2026-11-26. Plan and pre-registered criteria: `PLAN.md`. KARST was not
modified; nothing was updated in any fit.

## Verdict

**Step 0: S0a passed, S0b failed, so the study stopped.** The marginal
distributions do carry information about selection, but KARST's fit already
explains the shift between its two time points as well as the true landscape
does. No step-1 update was judged.

## Step 0 (development half; `results/step01/`)

PL per cell at the last time point, medians over datasets:

| Family | Datasets | Truth − neutral (S0a) | Truth − KARST (S0b) | KARST − neutral | Truth above KARST |
|---|---|---|---|---|---|
| Own (dev) | 168 | +0.102 | **−0.017** | +0.250 | 32% |
| Paper (odd replicates) | 925 | +0.144 | **−0.009** | +0.154 | 36% |

- **S0a:** the true landscape explains the change of the chromosome
  distributions better than neutral, so the marginals carry selection
  information.
- **S0b:** KARST's rates explain the same change at least as well as the
  truth's, in about two thirds of the datasets.
  - KARST was fitted to these two samples, so its PL is in-sample.
  - The information the marginals hold about this interval is already
    absorbed by the fit.
  - Pushing the PL further would fit the same two samples harder; it would
    not add information.

**The exception is the dispersed simulations with β too low** (paper,
medians):

| Group | Truth − neutral | Truth − KARST | KARST − neutral |
|---|---|---|---|
| vM 2e-4, β fixed | 0.178 | +0.004 | 0.143 |
| vM 4e-4, β fixed | 0.089 | +0.009 | 0.074 |
| vM 8e-4, β fixed | 0.022 | **+0.028** | **−0.002** |
| vM 8e-4, β true | 0.090 | +0.016 | 0.018 |
| vM 4e-4, β true | 0.099 | −0.014 | 0.119 |

- In the SA609-like regime (8e-4 with β fixed), KARST explains the marginal
  shift worse than neutral. There, the unused marginal information is real.
- Most of the gap closes when the fit uses the true β (4e-4: −0.014). It
  remains at 8e-4 (+0.016).
- So what is left in the marginals sits in the same place as before: β is
  too low.

**Salehi (descriptive; in-sample, no truth):**
- KARST's PL is above neutral in 27 of the 34 windows. Medians: SA609
  +0.38, other +0.55, SA906 +0.09.
- In SA609, only 7.7% of the first time point's cells are in the domain
  (other 82%, SA906 96%), so the fit's own likelihood sees few of them.
- S0b cannot be computed on Salehi.

## Notes

- 7 of the 1,875 paper inputs failed in the analysis (main λ 0.2, almost
  monomorphic inputs). They are left out.
- **Cost:** 7.0 CPU-hours (own development refits 4.9, paper analysis 2.1).
  The held-out half (3,500 refits) was not run.
- **Files:**
  - `R/mpl.R` (with a gradient check in the session);
  - `scripts/own_task.R`, `scripts/paper_task.R`, `scripts/step01_analyse.R`,
    `scripts/salehi_pl.R`;
  - `results/step01/`: dev_rows, step0, verdict, salehi_pl;
  - `results/step01_log.txt`;
  - on the HPC: `/share/lab_crd/taoli/Project/karst_patch/marginal_pl/`.
