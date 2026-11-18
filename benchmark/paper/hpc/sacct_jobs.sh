#!/usr/bin/env bash
# Peak memory and run time of every array task of the jobs in job registries
# (BENCH_RESULTS/jobs.tsv, written by hpc/submit.sh; hpc/jobs/<commit>.tsv for
# earlier rounds), from Slurm's accounting, for scripts/collect_resources.R.
# sacct runs on the host, not in the image.
#
#   hpc/sacct_jobs.sh OUT REGISTRY ...
#
# OUT: one row per array task and attempt: registry, job id, array index,
# state, peak resident memory of the batch step (KB), elapsed time, node and
# the task file (relative paths are relative to BENCH_DIR). The accounting
# samples memory every 30 s or so: a task shorter than that may show only the
# batch shell's few MB (run_task.R records its own peak, task_peak.Rds).
set -euo pipefail
cd "${BENCH_DIR:-$PWD}"
source hpc/env.sh
load_hpc_environment
out="$1"
shift
mkdir -p "$(dirname "${out}")"
tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT
printf 'registry\tjob_id\tarray_index\tstate\tmaxrss_kb\telapsed\tnode\ttask_file\n' >"${tmp}"
for registry in "$@"; do
  name="$(basename "${registry}" .tsv)"
  tail -n +2 "${registry}" | while IFS=$'\t' read -r job mode task_file _; do
    [[ "${mode}" == run ]] || continue
    # --duplicates: also the attempts before a requeue (e.g. after NODE_FAIL)
    sacct -j "${job}" --duplicates --format=JobID,State,MaxRSS,Elapsed,NodeList -n -P | awk -F'|' -v reg="${name}" \
      -v tf="${task_file}" 'BEGIN { OFS = "\t" }
      $1 ~ /\.batch$/ {
        id = $1; sub(/\.batch$/, "", id); n = split(id, p, "_")
        rss = $3; kb = rss + 0
        if (rss ~ /M$/) kb = kb * 1024; else if (rss ~ /G$/) kb = kb * 1024 * 1024
        print reg, p[1], (n > 1 ? p[2] : ""), $2, kb, $4, $5, tf
      }' >>"${tmp}"
  done
done
mv "${tmp}" "${out}"
trap - EXIT
echo "$(($(wc -l <"${out}") - 1)) task attempts written to ${out}"
