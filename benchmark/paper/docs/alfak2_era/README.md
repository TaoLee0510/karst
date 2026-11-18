# Reports of the alfak2 era

The paper benchmark's reports from before KARST was renamed (2026-11-11 to
2026-11-17), kept as they were written. They use the names of that time:

- **alfak2** is KARST (package `alfak2` up to version 0.2.0); its method ids
  `alfak2`, `alfak2_pairwise`, `alfak2_state_space` are `karst`,
  `karst_pairwise`, `karst_state_space` now.
- **alfakR** is alfak: "alfakR dev2 f93b2d2" is alfak dev2@f93b2d2, the
  reference fits; the paper's own fits are alfak main@4ad82b0.
- **Commit ids** are alfak2's ids from before its history was cleaned of
  attribution lines (trees unchanged). `doc/commit_map.tsv` at the top of the
  repository maps each to its current id, e.g. 05a8855 → 2ac72b2 and
  348addd → 64814c7.

Their results are on the HPC in
`/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_<id>` (8327282:
`results/`), named by the old ids; alfak's fits of the KARST rounds are
linked from there.

| alfak2 commit (current id) | Report |
|---|---|
| 8327282 (0be61fa) | `report.md` |
| ff13cf4 (ef84165) | `report_ff13cf4.md` |
| 8512694 (8e4e209; multi-timepoint) | `report_8512694.md` |
| 05a8855 (2ac72b2; series methods) | `report_05a8855.md`, `fig2d_check_05a8855.md` |
| 827c4b5 (49aed09; `initial_state = "auto"`; a regression, not of record) | `report_827c4b5.md` |
| 348addd (64814c7; round of record before the rename) | `report_348addd.md` |
