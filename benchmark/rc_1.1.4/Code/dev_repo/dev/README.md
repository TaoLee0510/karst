# dev: automated run on overcoming KARST's limitations

The automated run (phases 0–3) approved by the user on 2026-11-19. The criteria and decision rules are in [PLAN.md](PLAN.md) (fixed before scoring, not modified). Literature and candidate approaches are in `doc/literature/karst_limitations_review.html` (local, not in git).

**Summary of rules**: no commits during the run; the package's `R/`, `src/`, `tests/` and defaults are not changed; model changes are patches in `dev/patches/` applied on b319cfd and built into a separate library; HPC working directory `/share/lab_crd/taoli/Project/karst_dev/` (code synced by `dev/hpc/sync.sh`, not via git); no cap on HPC, qos xxlarge, every array with an add-hours monitor.

## Directory

| Path | Contents |
|---|---|
| `PLAN.md` | pre-registered criteria and rules |
| `R/common.R` | shared functions: paths, e-set configuration, ABM with passage records, replay windows |
| `abm/abm_dev.cpp` | copy of the e-set ABM: records every passage; also RMF and NK landscapes |
| `hpc/` | `env.sh`, `submit.sh` (array + monitor), `array_task.sh`, `build.sh`, `make_src.sh` (baseline + patches → source package), `sync.sh` |
| `phase0/` … `phase3/` | scripts, task tables and Chinese-language reports for each phase |
| `patches/` | model changes |
| `results/` | local results (not in git) |

## Status

The time is when the row was written (EST, checked against the session records). Rows before 11-19 08:24 were originally labelled with estimated times that were 1–4 hours late; they have been corrected.

| Time (EST) | Item |
|---|---|
| 11-19 02:44 | PLAN.md fixed |
| 11-19 02:44 | Infrastructure: the dev ABM replays the e1 windows bit-for-bit identical to the original data; HPC base build `karst_base` (b319cfd) |
| 11-19 02:44 | D6+D4 process noise: run locally (neutral 240, replay 264, replicates 2160 tasks) |
| 11-19 02:52 | D6: c = 1.777 (G0 passed). D4: median dispersion of rare karyotypes at bottleneck 200 is 19.3 (G1 passed, R1 enters phase 1); the process variance of rare karyotypes is 0.47 times the drift prediction (G2 not flagged, phase 1 does not use the R5 fallback) |
| 11-19 02:52 | HPC: base_x1 22065376, base_x4 22065378, r1_x1 22067517, r1_x4 22067519 (1056 tasks each; monitor +1); R1 patch `patches/R1_nb1_drift.patch` (posterior information switched to the NB1 expected information) |
| 11-19 02:52 | D1, D2 handed to forked sub-agents (`phase0/d1_*`, `phase0/d2_*`) |
| 11-19 03:25 | Phase 0 done (`phase0/report.md`): D3 the high β comes from the prior and the determinant; D2 the far-away ranking reaches only the ceiling 0.27–0.36 (R3 must assess ranking); D5 Laplace/NUTS sd ratio 0.85; D1 RMF, NK passed |
| 11-19 03:25 | R1 scoring: all five failed, the R1 line is ended (`phase1/report.md`). Also: 1.1.0 reports `karst_inner_failure` on 2–6% of e-set datasets |
| 11-19 03:25 | In progress: R2(i) 22072764; R2(iii) ×1/×4/×16 22080995/22080997/22080999; paper ABM missegregation ×1/×4/×16 22084547/49/51; R4 misjudgement base/R4 22080370/22080374, clean 22080372/22080376, ABM w2 base/R4 22084543/22084545 |
| 11-19 03:25 | General tools: `R/adapters.R` (unified interface across evaluation sets), `phase1/adapter_fit.R`, `R/calibration.R`, paper ABM interface `R/paper_abm.R` (consistent with 2c60af1) |
| 11-19 03:42 | Phase 2: R3 patch `patches/R3_shell2.patch` (22086568 base / 22086570 R3), R9 patch `patches/R9_revert_link.patch` (stacked on R3, 22087219); `phase2/far_fit.R` assesses near and far at the same time. R3 and R9 ran out of memory in large numbers at 8 GB; to be rerun at 32 GB after the arrays finish (a resource problem, not a method failure) |
| 11-19 03:42 | Phase 3: R7 (history = recorded population age; 22087599/22087601/22087603/22087605), R6 β grid ×1/16, ×1/4, ×64 (22089769/22089871/22089873), missegregation simulations R2 (22088245/47/49, R2(i) 22088551); Lynch 2022 Table 1 saved as `phase3/lynch2022_table1.csv` |
| 11-19 03:42 | Forked sub-agents: R6 dosage prior, R5, R8(a), R10 (each spec written in `phase3/*_spec.md`) |
| 11-19 03:42 | R8(b) cross-family conformal: failed (the error / model sd scales of NK and GRF/RMF differ too much) |
| 11-19 04:39 | Phase 1 done: R1, R2(i), R2(iii), R4 all failed or not kept (`phase1/report.md`) |
| 11-19 04:39 | Phase 3: R6 β prior failed, R6 dosage prior failed, R7 failed, R8(b) failed, R10 failed (WGD cannot be assessed); R5, R8(a) in progress; the out-of-memory tasks of R3/R9 rerun at 32 GB (22096740, 22096742) |
| 11-19 04:39 | Infrastructure fixes: `sync.sh` no longer deletes task tables generated on HPC; `array_task.sh` retries start-up failures within 30 seconds |
| 11-19 08:24 | R3/R9 32 GB rerun: still 150/153 out of memory. Cause: the shell-2 aggregated counts link all shell-1 latent variables together, so the Laplace Hessian becomes dense (a cost of the method itself). Rerun at 256 GB: R3 22100757, R9 22100759 |
| 11-19 08:24 | Interim assessment on a subset (R3 162/312, R9 159/312): R3 and R9 both failed (far-away bias almost unchanged, d3–4 ranking −0.015 / −0.042, intervals excluding 0) |
| 11-19 08:24 | Infrastructure: in phase 2, 52 fits of base and R9 on RMF/NK failed because multiple tasks shared the Rcpp compile cache (same as R8(a)); `R/common.R` changed to one cache per process, failed outputs deleted and rerun: base 22100969, R9 22100971 |
| 11-19 08:24 | Phase 3 report filled in for R6 (β prior, dosage prior), R5, R10, all failed; R8(a) round 1 22099022, round 2 22099024 in progress |
| 11-19 12:09 | R8(a): failed (GRF far-away coverage 0.87, 11 cells <0.90; RMF too wide). Phase 3 fully decided (`phase3/report.md`); waiting for the 256 GB R3/R9 reruns and the base/R9 cache-failure reruns |
| 11-20 00:19 | The add-hours monitors for the R3/R9 256 GB reruns (22100758, 22100760) were placed on xxlarge and were cancelled at 11-19 20:20 by their own 12-hour limit. When this was found at 22:28, the 12 tasks still running were about 50 minutes from the 15-hour limit; they were manually extended to 18 hours with `add_hours_to_running_array_jobs.sh`, and the monitors were resubmitted on small QOS with a 3-day duration (22130270, 22130271). No task was lost to timeout |
| 11-20 00:19 | User decision: the out-of-memory R3/R9 tasks are not rerun again (HPC nodes have about 1 TB, so technically another try is possible). The cases out of memory at 256 GB were those on GRF with many cells or long intervals; very_sparse all completed |
| 11-20 00:19 | Phase 2 report first draft (`phase2/report.md`): R3, R9 both failed (on the datasets where all three methods completed). Overall report first draft `REPORT.md`: no candidate passed in phases 0–3. 5 R3 tasks and 1 R9 task still running; numbers to be updated when they finish |
| 11-20 05:08 | Placeholders in the phase 2 report filled in: R3 22100757 completed 74, out of memory 76; R9 22100759 completed 78, out of memory 75; after recomputing with all fits the verdicts are unchanged, and on λ 1.2 the ranking improvement of R3 narrows to 0.28 → 0.32. |
