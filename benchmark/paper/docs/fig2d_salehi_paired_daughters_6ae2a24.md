# Forecasting daughter karyotypes in Salehi Fig. 2d

2026-11-24 · KARST benchmark `6ae2a24` · Default missegregation rate $\beta=5\times10^{-5}$

Across the 45 paired daughters from the original Fig. 2d, KARST predicted the observed karyotype distribution with a **smaller Wasserstein distance in 37** and **higher karyotype-frequency overlap in 35** than published ALFA-K. The difference was largest in the SA609 family. Of the 18 daughters in that family for which KARST beat ALFA-K on both metrics, **13 also beat the parental-frequency baseline on both metrics and five scored worse on both**.

## Cohort and scoring

The original Fig. 2d uses Spearman correlation to assess whether a model ranks daughters by the **magnitude of change from their parents**. It compares parent-to-forecast distance with parent-to-observation distance. Here we measure the **distance from forecast to observation** and the frequency mass assigned to the same karyotypes. Ranking changes across daughters and predicting each daughter’s composition can therefore favor different models. [Missegregation grid and Fig. 2d definition][grid]

**Pairing.** We retained each original training window, `minobs` threshold, and forecast day, and verified each parent–daughter pair by UID. Of the original 46 daughters, UID 75 in window `SA609UnBU_X4_l_3_d1_1_d2_2` lacked a paired KARST forecast. The remaining 45 unique daughters came from 34 parents and 34 windows. We scored 38 daughters on day 15 and seven on day 25. The original `lineage_type` labels assign 19 daughters to the SA609 family and 26 to other lineages. Each daughter counts once, including daughters that share a parent. [Original 46-point list][cohort]

**From inputs to scores.** Each row of the Salehi input matrix `x` is a complete 22-autosome copy-number vector. Each column contains cell counts at one passage. Dividing each daughter count by the total number of daughter cells gives its observed frequency distribution $Q$. For each model, we averaged the karyotype frequencies from five simulations at the selected forecast day with equal weight to obtain $P$. All observed daughter cells enter the evaluation, with zero forecast frequency for an unpredicted karyotype. We evaluated both the archived published ALFA-K forecasts and the dev2 version (`f93b2d2`) used in the missegregation grid. [Model versions][versions] · [Forecast extraction script][export]

| Metric | Calculation | Interpretation |
|---|---|---|
| Wasserstein $W_1$ | Minimum transport cost $\sum_{i,j}\pi_{ij}\lVert k_i-k_j\rVert_2$ | Transport preserves the forecast and observed masses. The cost uses Euclidean distance between 22-dimensional karyotypes. Lower is better. |
| Frequency overlap $O$ | $O(P,Q)=\sum_k\min(P_k,Q_k)$ | Only an exact match of all 22 entries counts. Values range from 0 to 1. Higher is better. |

## 1. Which model better predicts daughter composition?

For each daughter, we compared the KARST and ALFA-K forecasts against the same observed distribution $Q$. Each cell gives **KARST wins / ALFA-K wins**. The final column counts daughters won on both metrics. All comparisons used unrounded scores with a tie tolerance of $10^{-8}$. There were no ties.

| Lineages | N | Comparator | Lower $W_1$ | Higher overlap | Won both |
|---|---:|---|---:|---:|---:|
| SA609 family | 19 | Published ALFA-K | 18 / 1 | 18 / 1 | 18 / 1 |
| SA609 family | 19 | ALFA-K dev2 | 18 / 1 | 18 / 1 | 18 / 1 |
| Other lineages | 26 | Published ALFA-K | 19 / 7 | 17 / 9 | 15 / 5 |
| Other lineages | 26 | ALFA-K dev2 | 17 / 9 | 17 / 9 | 13 / 5 |
| All | 45 | Published ALFA-K | 37 / 8 | 35 / 10 | 33 / 6 |
| All | 45 | ALFA-K dev2 | 35 / 10 | 35 / 10 | 31 / 6 |

**The same 18 SA609 daughters favored KARST on both metrics.** The exception was SA609 X9 → X10, with $W_1$ of 0.759 for published ALFA-K and 0.869 for KARST, and overlap of 56.2% and 52.8%, respectively. Among the other 26 daughters, KARST won both metrics in 15, published ALFA-K won both in five, and the two metrics favored different methods in six. [Paired results for all 45 daughters][paired45]

**Why do the archived scores give different counts?** The original scoring filtered observed cells through each model’s own karyotype list before renormalization. In 19 of 45 daughters, the two lists retained different numbers of cells. Under archived scoring, KARST versus published ALFA-K gave 16 / 3 for both metrics in SA609, 15 / 11 for $W_1$ and 17 / 9 for overlap elsewhere, and 31 / 14 and 33 / 12 overall. The table above uses the full daughter distribution as a common target. [Win counts under each scoring convention][paircounts]

## 2. Do the forecasts beat carrying parent frequencies forward?

For the same 45 pairs, we read the count matrix for the two adjacent passages and verified that the parental counts matched the final training column karyotype by karyotype. The no-change forecast $B$ carries the parent’s full observed frequencies forward and is scored against the same $Q$. The table counts strictly lower $W_1$ or strictly higher overlap than baseline. Each model is compared independently with baseline, so multiple models can improve on the same daughter. Neither metric had ties.

| Lineages | N | Model | Lower $W_1$ | Higher overlap | Both better |
|---|---:|---|---:|---:|---:|
| SA609 family | 19 | KARST | 14 | 14 | 14 |
| SA609 family | 19 | Published ALFA-K | 7 | 4 | 4 |
| SA609 family | 19 | ALFA-K dev2 | 8 | 6 | 5 |
| Other lineages | 26 | KARST | 10 | 10 | 6 |
| Other lineages | 26 | Published ALFA-K | 6 | 8 | 1 |
| Other lineages | 26 | ALFA-K dev2 | 6 | 6 | 2 |
| All | 45 | KARST | 24 | 24 | 20 |
| All | 45 | Published ALFA-K | 13 | 12 | 5 |
| All | 45 | ALFA-K dev2 | 14 | 12 | 7 |

KARST beat the baseline on both metrics in 14 of all 19 SA609-family daughters. Thirteen came from the 18 daughters selected in Analysis 1. The fourteenth was SA609 X9 → X10, where published ALFA-K ranked first, KARST second, and no change third on both metrics (baseline $W_1=0.931$, overlap 48.2%). In the other lineages, KARST improved each metric in 10/26 daughters, with 6/26 improving on both. [Per-daughter baseline scores][base45] · [Grouped win counts][basecounts]

## 3. How large were the 13 gains and five losses in SA609?

We used the 18 daughters for which KARST beat published ALFA-K on both metrics in Analysis 1. For each daughter, $\Delta W_1$ is baseline distance minus KARST distance, and $\Delta O$ is KARST overlap minus baseline overlap multiplied by 100. Positive values indicate improvement on both scales. The former is an absolute distance reduction and the latter a percentage-point gain (pp). Arrows show baseline score → KARST score. Differences use unrounded values.

**Thirteen daughters improved on both metrics.** Median distance reduction was **0.293** (range 0.003–1.331), and median overlap gain was **3.27 pp** (0.95–12.66).

| Parent | Daughter | $W_1$ baseline → KARST | $\Delta W_1$ | Overlap (%) baseline → KARST | $\Delta O$ (pp) |
|---|---|---:|---:|---:|---:|
| SA609R2 X6 | SA609H2 X7 | 2.458 → 1.127 | +1.331 | 13.48 → 26.14 | +12.66 |
| SA609R2 X6 | SA609R2 X7 | 2.236 → 1.321 | +0.916 | 17.56 → 30.04 | +12.48 |
| SA609R2 X6 | SA609R2b X7 | 2.331 → 1.654 | +0.677 | 12.60 → 15.88 | +3.27 |
| SA609 X3 | SA609 X4 | 2.164 → 1.768 | +0.396 | 12.36 → 19.23 | +6.86 |
| SA609 X3 | SA609UnBU X4 | 2.491 → 2.151 | +0.340 | 5.66 → 7.09 | +1.44 |
| SA000 X5 | SA001 X6 | 1.652 → 1.335 | +0.317 | 16.41 → 25.73 | +9.33 |
| SA000 X5 | SA609H1 X6 | 1.672 → 1.379 | +0.293 | 17.39 → 29.91 | +12.51 |
| SA000 X5 | SA000 X6 | 1.851 → 1.622 | +0.229 | 10.88 → 12.23 | +1.35 |
| SA609UnBU X6 | SA609UnBU X7 | 2.256 → 2.109 | +0.146 | 7.27 → 8.23 | +0.95 |
| SA609 X6 | SA609 X7 | 1.538 → 1.425 | +0.112 | 23.64 → 25.81 | +2.16 |
| SA609 X7 | SA609 X8 | 1.383 → 1.280 | +0.104 | 28.97 → 33.94 | +4.97 |
| SA609 X4 | SA609 X5 | 1.712 → 1.617 | +0.095 | 12.35 → 14.23 | +1.88 |
| SA609 X5 | SA609 X6 | 1.686 → 1.683 | +0.003 | 14.98 → 17.85 | +2.87 |

**Five daughters scored worse on both metrics.** The median distance increase was **0.177** (range 0.098–0.279), and the median overlap loss was **2.91 pp** (0.32–4.13).

| Parent | Daughter | $W_1$ baseline → KARST | $\Delta W_1$ | Overlap (%) baseline → KARST | $\Delta O$ (pp) |
|---|---|---:|---:|---:|---:|
| SA000 X6 | SA001 X7 | 2.030 → 2.128 | -0.098 | 15.92 → 11.79 | -4.13 |
| SA609 X8 | SA609 X9 | 1.238 → 1.354 | -0.116 | 35.56 → 35.24 | -0.32 |
| SA000 X6 | SA000 X7 | 2.304 → 2.481 | -0.177 | 11.25 → 8.34 | -2.91 |
| SA000 X6 | SA609H1 X7 | 2.049 → 2.249 | -0.199 | 13.91 → 10.48 | -3.43 |
| SA609R2 X5 | SA609R2 X6 | 3.168 → 3.447 | -0.279 | 2.01 → 1.39 | -0.62 |

The three largest distance reductions all came from **SA609R2 X6** (1.331, 0.916, and 0.677). Three of the five losses came from **SA000 X6**. For SA609 X5 → X6, the distance fell by only 0.003 while overlap rose by 2.87 pp. The counts and medians describe the selected daughters, with multiple branches from the same parent contributing to the totals. [Full-precision differences for all 18 daughters][deltas]

## Reproducing the comparisons

Run the commands below from the KARST repository root with the R package `transport`, the local results archive, and the original ALFA-K Salehi data. `PAPER_DATA` defaults to `../ALFA-K/data` in a sibling checkout. The first script reads the saved five-replicate mean forecasts and cross-checks the archived scores. The second computes the baseline from the parent and daughter cell counts. Git ignores the local results directory. The result tables are included above for readers of the standalone document. [Mean forecasts][means] · [Archived scores][scores] · [Scoring script][recalc] · [Baseline script][basecode]

```sh
audit=benchmark/paper/results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints
Rscript "$audit/recalculate.R"
Rscript "$audit/compare_no_change.R"
```

This R code reconstructs the 18 daughter-level differences and the two sets of medians from the 45-row baseline table.

```r
audit <- "benchmark/paper/results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints"
points <- read.delim(file.path(audit, "no_change_45_points.tsv"))
sa609 <- subset(points, in_sa609_18)
sa609$w1_reduction <- with(sa609, no_change_w1 - karst_w1)
sa609$overlap_gain_pp <- with(sa609, 100 * (karst_overlap - no_change_overlap))
sa609$both_better <- with(sa609, w1_reduction > 0 & overlap_gain_pp > 0)
stopifnot(nrow(sa609) == 18, sum(sa609$both_better) == 13,
          with(sa609, all(sign(w1_reduction) == sign(overlap_gain_pp))))
sa609[order(-sa609$w1_reduction),
      c("parent_sample", "daughter_sample", "w1_reduction", "overlap_gain_pp")]
aggregate(cbind(w1_reduction, overlap_gain_pp) ~ both_better, sa609, median)
```

[grid]: beta_grid_salehi_6ae2a24.md
[versions]: ../README.md
[cohort]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paper_46_observations.tsv
[scores]: ../results_6ae2a24/tables/salehi_scores.Rds
[export]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/export_red_forecasts.R
[means]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/forecast_means.rds
[paired45]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/paired_45_published_alfak.tsv
[paircounts]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/paired_win_counts.tsv
[base45]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/no_change_45_points.tsv
[basecounts]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/no_change_win_counts.tsv
[deltas]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/sa609_18_karst_vs_no_change_deltas.tsv
[recalc]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/recalculate.R
[basecode]: ../results_6ae2a24/beta_grid/analysis/fig2d_original_audit/paired_endpoints/compare_no_change.R
