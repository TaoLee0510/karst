# History choice study

Results: `/share/lab_crd/taoli/Project/karst/validation/sparse_joint/history_choice/results`, analysed 2026-11-19 13:54 EST.

## Development benchmark (two time points)

168 datasets read, 168 fitted, 0 failed; kept candidate identical to the fit returned: 168 of 168.

Candidates: 564; posterior-mode objective differs from nlminb's by more than 1e-3 at 127 (22.5%), found by the second search at 4.
Current rule reproduces the fit's own choice at 168 of 168 datasets.

### Rule A against the current rule

Choice changed at 2 of 168 datasets (core 2, stress 0).

| distance | datasets | median_rmse_ratio | q25 | q75 | better | worse | sign_p_worse | median_spearman_diff | median_coverage_diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d0 |    2 | 0.978 | 0.96 | 0.995 |    1 |    1 | 0.75 | 0.0282 | 0.0233 |
| d1 |    2 | 0.977 | 0.963 | 0.991 |    1 |    1 | 0.75 | 0.0204 | 0.0117 |
| d2 |    2 | 0.966 | 0.956 | 0.976 |    2 | 0 |    1 | 0.0248 | 0.0247 |

### Rule B against the current rule

Choice changed at 5 of 168 datasets (core 3, stress 2).

| distance | datasets | median_rmse_ratio | q25 | q75 | better | worse | sign_p_worse | median_spearman_diff | median_coverage_diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d0 |    5 |    1 | 0.991 |    1 |    4 |    1 | 0.969 | 0.00582 | 0 |
| d1 |    5 |    1 | 0.992 |    1 |    4 |    1 | 0.969 | 0.00816 | 0 |
| d2 |    5 | 0.997 | 0.987 |    1 |    4 |    1 | 0.969 | 0.00937 | 0.0005 |

Median RMSE ratio over distances 0-2 together: 1.000.

### Rule oracle against the current rule

Choice changed at 115 of 168 datasets (core 65, stress 50).

| distance | datasets | median_rmse_ratio | q25 | q75 | better | worse | sign_p_worse | median_spearman_diff | median_coverage_diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d0 |  115 | 0.994 | 0.968 |    1 |   93 |   22 |    1 | 0 | 0 |
| d1 |  115 | 0.989 | 0.965 |    1 |  101 |   14 |    1 | 0.00000401 | 0 |
| d2 |  115 | 0.984 | 0.933 |    1 |  105 |   10 |    1 | 0.000000769 | 0.0005 |

### All datasets, by rule

| rule | distance | set | datasets | mean_rmse | mean_rmse_zero | median_spearman | mean_coverage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| current | d0 | core |   96 | 0.0262 | 0.0463 | 0.495 | 0.98 |
| A | d0 | core |   96 | 0.0262 | 0.0463 | 0.495 | 0.981 |
| B | d0 | core |   96 | 0.0262 | 0.0463 | 0.495 | 0.981 |
| oracle | d0 | core |   96 | 0.0261 | 0.0463 | 0.501 | 0.982 |
| current | d1 | core |   96 | 0.0445 | 0.0918 | 0.422 | 0.991 |
| A | d1 | core |   96 | 0.0445 | 0.0918 | 0.423 | 0.991 |
| B | d1 | core |   96 | 0.0445 | 0.0918 | 0.423 | 0.991 |
| oracle | d1 | core |   96 | 0.044 | 0.0918 | 0.412 | 0.992 |
| current | d2 | core |   96 | 0.0662 | 0.137 | 0.34 | 0.984 |
| A | d2 | core |   96 | 0.0661 | 0.137 | 0.342 | 0.984 |
| B | d2 | core |   96 | 0.0661 | 0.137 | 0.342 | 0.984 |
| oracle | d2 | core |   96 | 0.0652 | 0.137 | 0.339 | 0.986 |
| current | d0 | stress |   72 | 0.0471 | 0.0696 | 0.489 | 0.901 |
| A | d0 | stress |   72 | 0.0471 | 0.0696 | 0.489 | 0.901 |
| B | d0 | stress |   72 | 0.0471 | 0.0696 | 0.489 | 0.901 |
| oracle | d0 | stress |   72 | 0.0465 | 0.0696 | 0.51 | 0.904 |
| current | d1 | stress |   72 | 0.0639 | 0.113 | 0.524 | 0.959 |
| A | d1 | stress |   72 | 0.0639 | 0.113 | 0.524 | 0.959 |
| B | d1 | stress |   72 | 0.0639 | 0.113 | 0.524 | 0.959 |
| oracle | d1 | stress |   72 | 0.0625 | 0.113 | 0.536 | 0.963 |
| current | d2 | stress |   72 | 0.0853 | 0.156 | 0.443 | 0.964 |
| A | d2 | stress |   72 | 0.0853 | 0.156 | 0.443 | 0.964 |
| B | d2 | stress |   72 | 0.0853 | 0.156 | 0.443 | 0.964 |
| oracle | d2 | stress |   72 | 0.0827 | 0.156 | 0.451 | 0.971 |

### Datasets whose choice A or B changes

History chosen by each rule; mean RMSE over distances 0-2 relative to the current rule's choice.

| dataset | set | current | A | B | rmse_ratio_A | rmse_ratio_B |
| --- | --- | --- | --- | --- | --- | --- |
| dev_L02_P01_shift_3000_100_360d | stress |  3000 |  3000 |   360 |     1 |     1 |
| dev_L02_P01_typical_500_700_360d | core |  1000 |  1000 |   360 |     1 |     1 |
| dev_L02_P01_unbalanced_100_3000_360d | core |   360 |  1000 |  1000 | 0.9466 | 0.9466 |
| dev_L07_P01_typical_500_700_360d | core |   360 |  1000 |  1000 | 0.9992 | 0.9992 |
| dev_L11_P01_shift_3000_100_360d | stress |  3000 |  3000 |  1000 |     1 | 0.9935 |

## ALFA-K paper benchmark: ABM and missegregation simulations

377 tasks read, 377 fitted, 0 failed; kept candidate identical to the fit returned: 377 of 377.

Group left out, no state in any fit: shell2 (a fit's domain is the majors and their one-step neighbours).

Candidates: 2071; posterior-mode objective differs from nlminb's by more than 1e-3 at 402 (19.4%), found by the second search at 67; candidates whose mode was not found at all: 0.
Current rule reproduces the fit's own choice at 377 of 377 tasks.

Objective at the posterior mode minus nlminb's (ABM candidates):

| candidates | n | median_gap | share_above_1e3 | share_below_1e3 |
| --- | --- | --- | --- | --- |
| all | 2071 | 0 | 0.163 | 0.0309 |
| second search |   67 | 8.95 |    1 | 0 |
| first search | 2004 | 0 | 0.135 | 0.0319 |

### Set random (306 tasks)

Tasks with a candidate rescued by the second search: 0; with a candidate whose mode was not found: 0.
Rescued candidate chosen: current 0, A 0, B 0.

Rule A against the current rule: choice changed at 0 of 306 tasks.

Rule B against the current rule: choice changed at 3 of 306 tasks.

| group | tasks | median_spearman_diff | mean_spearman_diff | better | worse | sign_p_worse | median_pearson_diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all |    3 | 0.00501 | 0.00553 |    3 | 0 |    1 | 0.00495 |
| major |    3 | 0 | 0 | 0 | 0 | NA | -0.0148 |
| shell1 |    3 | 0.0055 | 0.00589 |    3 | 0 |    1 | 0.00574 |

Changed by design and method:

| design / method | tasks changed |
| --- | --- |
| w2_ends karst |    2 |
| w8_all karst |    1 |

Rule oracle against the current rule: choice changed at 206 of 306 tasks.

| group | tasks | median_spearman_diff | mean_spearman_diff | better | worse | sign_p_worse | median_pearson_diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all |  206 | 0.00768 | 0.0277 |  206 | 0 |    1 | 0.0081 |
| major |  206 | 0 | 0.0144 |   70 |   50 | 0.973 | 0.000408 |
| shell1 |  206 | 0.011 | 0.0297 |  203 |    3 |    1 | 0.0115 |

All tasks of the set, by rule:

| rule | group | tasks | median_spearman | mean_spearman | mean_pearson |
| --- | --- | --- | --- | --- | --- |
| current | all |  306 | 0.782 | 0.699 | 0.707 |
| A | all |  306 | 0.782 | 0.699 | 0.707 |
| B | all |  306 | 0.782 | 0.699 | 0.707 |
| oracle | all |  306 | 0.803 | 0.718 | 0.725 |
| current | major |  306 | 0.937 | 0.867 | 0.908 |
| A | major |  306 | 0.937 | 0.867 | 0.908 |
| B | major |  306 | 0.937 | 0.867 | 0.908 |
| oracle | major |  306 | 0.942 | 0.877 | 0.914 |
| current | shell1 |  306 | 0.786 | 0.689 | 0.693 |
| A | shell1 |  306 | 0.786 | 0.689 | 0.693 |
| B | shell1 |  306 | 0.786 | 0.689 | 0.693 |
| oracle | shell1 |  306 | 0.802 | 0.709 | 0.715 |

Tasks whose choice A or B changes (history chosen by each rule; Spearman correlation over all domain states, and each rule's difference from the current rule):

| task | current | A | B | oracle | current_rescued | spearman_current | diff_A | diff_B | diff_oracle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| abm w_0p8_m_5e-05_rep_42 w8_all karst | 60.5 | 60.5 |  100 | 3000 | FALSE | 0.874 | 0 | 0.00501 | 0.0129 |
| misseg w_0p8_m_2e-04_rep_10 w2_ends karst | 3000 | 3000 |  100 |   30 | FALSE | 0.398 | 0 | 0.00182 | 0.0801 |
| misseg w_0p8_m_8e-04_rep_15 w2_ends karst | 1000 | 1000 | 3000 |   30 | FALSE | 0.24 | 0 | 0.00976 | 0.171 |

### Set failed (71 tasks)

Tasks with a candidate rescued by the second search: 66; with a candidate whose mode was not found: 0.
Rescued candidate chosen: current 3, A 0, B 0.

Rule A against the current rule: choice changed at 3 of 71 tasks.

| group | tasks | median_spearman_diff | mean_spearman_diff | better | worse | sign_p_worse | median_pearson_diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all |    3 | -0.015 | -0.0063 |    1 |    2 |  0.5 | -0.0162 |
| major |    3 | 0 | 0 | 0 | 0 | NA | -0.00116 |
| shell1 |    3 | -0.0155 | -0.00694 |    1 |    2 |  0.5 | -0.0169 |

Rule B against the current rule: choice changed at 6 of 71 tasks.

| group | tasks | median_spearman_diff | mean_spearman_diff | better | worse | sign_p_worse | median_pearson_diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all |    6 | -0.0117 | -0.00737 |    1 |    5 | 0.109 | -0.013 |
| major |    6 | 0 | -0.00303 | 0 |    1 |  0.5 | 0.000394 |
| shell1 |    6 | -0.0132 | -0.00821 |    1 |    5 | 0.109 | -0.0141 |

Changed by design and method:

| design / method | tasks changed |
| --- | --- |
| w2_ends karst |    1 |
| w8_all karst |    2 |
| w8_all karst_state_space |    1 |
| w8_ends karst |    2 |

Rule oracle against the current rule: choice changed at 50 of 71 tasks.

| group | tasks | median_spearman_diff | mean_spearman_diff | better | worse | sign_p_worse | median_pearson_diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all |   50 | 0.0141 | 0.021 |   50 | 0 |    1 | 0.0119 |
| major |   50 | 0 | 0.0117 |   21 |   12 | 0.96 | 0.00145 |
| shell1 |   50 | 0.0147 | 0.0242 |   50 | 0 |    1 | 0.0171 |

All tasks of the set, by rule:

| rule | group | tasks | median_spearman | mean_spearman | mean_pearson |
| --- | --- | --- | --- | --- | --- |
| current | all |   71 | 0.818 | 0.788 | 0.787 |
| A | all |   71 | 0.818 | 0.788 | 0.787 |
| B | all |   71 | 0.818 | 0.787 | 0.786 |
| oracle | all |   71 | 0.838 | 0.803 | 0.804 |
| current | major |   71 | 0.939 | 0.891 | 0.92 |
| A | major |   71 | 0.939 | 0.891 | 0.92 |
| B | major |   71 | 0.939 | 0.89 | 0.92 |
| oracle | major |   71 | 0.93 | 0.899 | 0.939 |
| current | shell1 |   71 | 0.82 | 0.785 | 0.785 |
| A | shell1 |   71 | 0.82 | 0.785 | 0.784 |
| B | shell1 |   71 | 0.82 | 0.784 | 0.784 |
| oracle | shell1 |   71 | 0.848 | 0.802 | 0.805 |

Tasks whose choice A or B changes (history chosen by each rule; Spearman correlation over all domain states, and each rule's difference from the current rule):

| task | current | A | B | oracle | current_rescued | spearman_current | diff_A | diff_B | diff_oracle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| abm w_0p2_m_5e-05_rep_97 w8_all karst | 1000 | 1000 | 62.8 |  100 | FALSE | 0.406 | 0 | -0.000596 | 0.000243 |
| abm w_0p4_m_5e-05_rep_15 w8_all karst | 1000 | 1000 |  300 | 1000 | FALSE | 0.664 | 0 | -0.0165 | 0 |
| abm w_0p4_m_5e-05_rep_82 w8_all karst_state_space |  100 | 49.5 | 49.5 | 1000 | TRUE | 0.73 | -0.015 | -0.015 | 0.00358 |
| abm w_0p8_m_5e-05_rep_37 w2_ends karst | 1000 | 1000 |  100 | 1000 | FALSE | 0.895 | 0 | -0.00829 | 0 |
| abm w_1p6_m_5e-05_rep_32 w8_ends karst |  100 | 68.7 | 68.7 | 68.7 | TRUE | 0.881 | 0.0123 | 0.0123 | 0.0123 |
| abm w_1p6_m_5e-05_rep_36 w8_ends karst |  100 | 64.2 | 64.2 |  300 | TRUE | 0.849 | -0.0162 | -0.0162 | 0.00715 |

## Criteria for rule B (agreed 2026-11-19)

- Development benchmark, median RMSE ratio at most 1.00 over distances 0-2 and at most 1.02 at each: met (overall 1.000; d0 1.000, d1 1.000, d2 0.997).
- Development benchmark, median coverage at most 0.01 lower: met.
- ABM (random set), median Spearman difference at least 0 (all domain states): met (0.0050).
- ABM (random set), not significantly worse by a sign test: met (better 3, worse 0, one-sided p 1).
