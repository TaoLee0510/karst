# Dispersion gate for the missegregation rate: plan

Approved by the user on 2026-11-23 at about 18:44 EST ("Proceed as recommended, confirmed"):
D1a, D2a, D3a, D4a, D5a. Steps 0–2 run on their own and stop at the first
failed criterion; the run stops before step 3 (real data) in any case.

## Idea

KARST fits with a fixed missegregation rate, β = 5e-5. On ABM_varyMisseg
(`benchmark/paper/docs/fig2d_misseg_6ae2a24.md`):
- a β 8–16 times below the true rate makes KARST's landscape much worse;
- a β 2–4 times too low costs little;
- an earlier check (`../beta_choice/README.md`, step 2) found that a β too
  high costs more than one too low.

The gate does not estimate β. It measures how spread out the two input time
points are, compares that with a band calibrated on simulations with known
rates, and raises β by one conservative step only when the data are clearly
more dispersed than the default rate produces. Otherwise β stays 5e-5.

This reopens the "β from the data" line, which was closed on 2026-11-20, by
the user's decision. Earlier attempts v6/v7 (`dev/realdata/noise`) failed
because their statistic depended on the landscape (wavelength) and on the
simulation source more than on β. Step 0 tests that first, without fits.

## Statistics

Both use only the two input time points, t0 and t1, of a dataset. Each time
point is subsampled without replacement to n cells, for n = 100 and n = 300
(10 draws, seeds 1–10, mean over draws). A dataset with fewer than n cells
at either time point has no value at that n.

- **once**: the share of the 2n subsampled cells whose karyotype occurs
  exactly once in its own subsample (t0 and t1 pooled).
- **new**: the share of the n t1 cells whose karyotype does not occur in the
  t0 subsample.

Simpson diversity and the share of the most frequent karyotype are reported,
not used: on the paper's ABM at β = 5e-5 they vary with the wavelength as
much as with the rate.

## Step 0: does dispersion line up across sources? (no fits)

**Data** (all existing, read only):

| source | datasets | rates (multiples of 5e-5) | wavelengths |
| --- | --- | --- | --- |
| paper | main set 400 sims + ABM_varyMisseg 125 sims; inputs w2, w4, w8 ends (1000 cells stored per time point) | main ×1; varyMisseg ×1, ×2, ×4, ×8, ×16 | main 0.2–1.6; varyMisseg 0.8 |
| own | `sparse_two_timepoint_6ae2a24` round, 7000 datasets less the 50/50 scenario | ×1 | 0.2–1.6 |
| v7_e1, v7_e2, v7_e3, v7_c2 | `dev/realdata/noise/v7` simulations (both halves), less the 50/50 depth pair | ×1, ×4, ×16 | 0.2–1.6 |

A **cell** is source × wavelength (the paper's varyMisseg sims are their own
source, `paper_vm`, with one wavelength).

**Criterion S0** (per statistic and n): let τ be the smallest value such
that, in every cell, at most 5% of the ×1 datasets have a statistic > τ.
S0 passes for that statistic and n if, pooled over all sources, at least 50%
of the ×16 datasets exceed τ. Step 0 passes if S0 passes for at least one
statistic and n. Every cell's share above τ at ×1, ×4 and ×16 is reported.

If step 0 fails, the run stops with a report; nothing is fitted.

## Step 1: calibration (development data)

Fixed before step 1 is run, after step 0:
- Data: ABM_varyMisseg (all 125, inputs subsampled to the Salehi cell
  counts as in `fig2d_misseg`), and the v7 dev half at depth pairs 200/200
  and 500/700.
- Fits: KARST only (round library 1.1.3) at β = 5e-5 × {1, 2, 4, 8, 16}.
- Score: landscape accuracy against the truth (r over all karyotypes, and
  at distance 0, 1, 2), and RMSE at distance 1 against the net rate at the
  true rate (as `../beta_choice` after its 2026-11-18 scoring fix).
- Rule: statistic ≥ τ → β = ×k, else ×1. τ comes from step 0, k is chosen on
  this data among {2, 4, 8} (D3a: one step), the largest k whose median
  error at ×1 truth does not exceed the default's by more than 2%.
- Pick step (D2a): for datasets the gate raises, compare the default and
  the raised fit by KARST's cross-validation score, by the objective, and by
  a cell-half hold-out (fit on half of each time point's cells, score the
  other half). A criterion passes if, among datasets whose two fits differ
  in accuracy by more than 0.02 in r, it picks the more accurate fit in at
  least 70%, and its picks are on average not less accurate than the gate
  alone. If none passes, the gate alone decides.

### Step 1 rules, fixed on 2026-11-23 at about 19:10 EST, before any step-1 fit was read

Step 0 passed on `once` only (`new` fails: τ = 1). The gate is:
- n = 300 when both time points have at least 300 cells, else n = 100 when
  both have at least 100, else no gate (default β);
- raise when once ≥ τ_n, with τ_300 = 0.4015 and τ_100 = 0.6020 (step 0).

Errors:
- v7: RMSE at distance 1 against the net rate at the true rate
  (`stt_score_predictions`), per dataset; ratio_k = RMSE(×k) / RMSE(×1).
- ABM_varyMisseg: 1 − r over the karyotypes of KARST's landscape (as
  `fig2d_misseg`), per input and minobs; ratio_k = (1 − r_k) / (1 − r_1).

Choice of k ∈ {2, 4, 8}, on the datasets the gate raises:
- H(k) = median ratio_k at truth ×1 (false raises), v7 and varyMisseg pooled;
- G(k) = median ratio_k at truth ≥ ×4, v7 only (the package's metric);
  varyMisseg's G is reported and must be below 1 for the chosen k.
- k = the value with the lowest G(k) among those with H(k) ≤ 1.25.

Step 1 fails, and the run stops, if any of these holds:
- no k has H(k) ≤ 1.25 and G(k) ≤ 0.90;
- for the chosen k, varyMisseg's G is not below 1;
- the gate raises more than 10% of the ×1 datasets in any development
  cell (source × wavelength; varyMisseg at the Salehi cell counts included).

Pick step (D2a), on raised datasets, default fit against ×k fit:
- criteria: the objective, the posterior objective (both all data, lower
  wins), and the cell-half hold-out (`step1_pick_v7.R`: fit on half A of each
  time point's cells, multinomial NLL of half B over the fit's domain; only
  where both fits have the same domain);
- a criterion passes if, among raised v7 datasets whose two fits differ by
  more than 2% in RMSE at distance 1, it picks the more accurate fit in at
  least 70%, and the median RMSE ratio (pick / gate alone) is at most 1;
- if none passes, the gate alone decides (no pick step in step 2).

"Is closer to the truth better?" is answered descriptively: median ratio by
truth multiple and fitted multiple, v7 and varyMisseg.

## Step 2: confirmation (new simulations, read once)

- Simulator code copied from `dev/` into this folder (`dev/` untouched);
  seeds from 161001; v9 left unused.
- Own simulator and C2 passaging; wavelengths 0.2, 0.4, 0.8, 1.2, 1.6;
  rates ×1, ×4, ×16; 100, 300 and 1000 cells per time point.
- A report-only arm with 1% of cells given a random one-chromosome error.

| | criterion |
| --- | --- |
| C1 | ×1 raised in ≤ 5% overall and ≤ 10% in every cell |
| C2 | at ×1, median error ratio (gate / default) ≤ 1.02 |
| C3 | at ×16, raised in ≥ 50% and median error ratio ≤ 0.8 |
| C4 | β used > 2 × true rate in ≤ 5% |
| C5 | the pick step, if kept, meets its step-1 criterion |

## Step 3 (only after a stop and the user's approval)

Salehi: gate decisions per line and trajectory; agreement of parallel
lineages ≥ 70%; later passages as truth (log score with coverage) in a 2×2
of fit β × forecast rate. GBM / NUMBAT: decisions reported only.

## Rules

- Stop and report at the first failed criterion.
- No package change, `dev/` read only, no commit.
- HPC: one fit per task, qos xxlarge, the excluded nodes (with
  hpctpa3pc0003), an add-hours monitor per job.
- At most 3 heavy local processes.
- Results under `results/step<N>/` here (not tracked).
- Budget about 650 CPU-h; ask before going over 1,000.
