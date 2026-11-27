#!/bin/bash
# Submit as:
#   submit_hourly_add_hours.sh --time=12:00:00 ARRAY_ID ADD_HOURS
#
# Optional overrides:
#   submit_hourly_add_hours.sh --time=12:00:00 --threshold-hours=2 ARRAY_ID ADD_HOURS
#   sbatch --export=ALL,THRESHOLD_HOURS=4 submit_hourly_add_hours.sh ARRAY_ID ADD_HOURS
#   sbatch --export=ALL,DRY_RUN=TRUE submit_hourly_add_hours.sh ARRAY_ID ADD_HOURS

#SBATCH --job-name=hourly_add_hours
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --output=hourly_add_hours_%j.out
#SBATCH --error=hourly_add_hours_%j.err

set -euo pipefail

# rc freeze copy: the companion script is Code/monitor/add_hours_to_running_array_jobs.sh
ADD_HOURS_SCRIPT="${ADD_HOURS_SCRIPT:-${RC_HOME:?RC_HOME is not set}/Code/monitor/add_hours_to_running_array_jobs.sh}"
CHECK_INTERVAL_SECONDS="${CHECK_INTERVAL_SECONDS:-3600}"
THRESHOLD_HOURS="${THRESHOLD_HOURS:-6}"
DRY_RUN="${DRY_RUN:-FALSE}"
AUTO_CONFIRM="${AUTO_CONFIRM:-TRUE}"
WALLTIME="${WALLTIME:-7-00:00:00}"
MAX_MISSING_CHECKS="${MAX_MISSING_CHECKS:-3}"

usage() {
  cat <<'USAGE'
Usage:
  submit_hourly_add_hours.sh [options] ARRAY_ID ADD_HOURS
  sbatch [sbatch-options] submit_hourly_add_hours.sh ARRAY_ID ADD_HOURS

Arguments:
  ARRAY_ID   Parent SLURM array job id.
  ADD_HOURS  Number of hours passed to add_hours_to_running_array_jobs.sh.

Options:
  --time WALLTIME, --time=WALLTIME
                           Walltime passed to sbatch when the script is run
                           directly. Default: 7-00:00:00
  --threshold-hours N, --threshold-hours=N
                           Set THRESHOLD_HOURS. Default: 6
  --interval-seconds N, --interval-seconds=N
                           Set CHECK_INTERVAL_SECONDS. Default: 3600
  --dry-run TRUE|FALSE, --dry-run=TRUE|FALSE
                           Set DRY_RUN. Default: FALSE
  --auto-confirm TRUE|FALSE, --auto-confirm=TRUE|FALSE
                           Set AUTO_CONFIRM. Default: TRUE
  --max-missing-checks N, --max-missing-checks=N
                           Exit after this many consecutive checks cannot find
                           ARRAY_ID in squeue. Default: 3

Environment:
  ADD_HOURS_SCRIPT         Script to run each interval.
                           Default: $RC_HOME/Code/monitor/add_hours_to_running_array_jobs.sh
  CHECK_INTERVAL_SECONDS   Seconds between runs.
                           Default: 3600
  THRESHOLD_HOURS          Passed to add_hours_to_running_array_jobs.sh.
                           Default: 6
  DRY_RUN                  Passed to add_hours_to_running_array_jobs.sh.
                           Default: FALSE
  AUTO_CONFIRM             Passed to add_hours_to_running_array_jobs.sh.
                           Default: TRUE
  MAX_MISSING_CHECKS       Consecutive missing squeue checks before exit.
                           Default: 3

Resources:
  time=7-00:00:00, cpus-per-task=1, mem=1G, no qos specified.

Examples:
  Code/monitor/submit_hourly_add_hours.sh --time=12:00:00 --threshold-hours=2 17206553 2
  sbatch --time=12:00:00 --export=ALL,THRESHOLD_HOURS=2 Code/monitor/submit_hourly_add_hours.sh 17206553 2
USAGE
}

require_positive_integer() {
  local label="$1"
  local value="$2"

  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    echo "${label} must be a positive integer, got: ${value}" >&2
    exit 1
  fi
  if (( 10#${value} <= 0 )); then
    echo "${label} must be > 0, got: ${value}" >&2
    exit 1
  fi
}

require_walltime() {
  local value="$1"

  if ! [[ "${value}" =~ ^([0-9]+-)?[0-9]+(:[0-9]{1,2}){0,2}$ ]]; then
    echo "WALLTIME must be a Slurm time value like 12:00:00 or 7-00:00:00, got: ${value}" >&2
    exit 1
  fi
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

args=()
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --time=*)
      WALLTIME="${1#*=}"
      ;;
    --time)
      shift
      if [[ "$#" -eq 0 ]]; then
        echo "--time requires a value" >&2
        exit 1
      fi
      WALLTIME="$1"
      ;;
    --threshold-hours=*)
      THRESHOLD_HOURS="${1#*=}"
      ;;
    --threshold-hours)
      shift
      if [[ "$#" -eq 0 ]]; then
        echo "--threshold-hours requires a value" >&2
        exit 1
      fi
      THRESHOLD_HOURS="$1"
      ;;
    --interval-seconds=*)
      CHECK_INTERVAL_SECONDS="${1#*=}"
      ;;
    --interval-seconds)
      shift
      if [[ "$#" -eq 0 ]]; then
        echo "--interval-seconds requires a value" >&2
        exit 1
      fi
      CHECK_INTERVAL_SECONDS="$1"
      ;;
    --dry-run=*)
      DRY_RUN="${1#*=}"
      ;;
    --dry-run)
      shift
      if [[ "$#" -eq 0 ]]; then
        echo "--dry-run requires a value" >&2
        exit 1
      fi
      DRY_RUN="$1"
      ;;
    --auto-confirm=*)
      AUTO_CONFIRM="${1#*=}"
      ;;
    --auto-confirm)
      shift
      if [[ "$#" -eq 0 ]]; then
        echo "--auto-confirm requires a value" >&2
        exit 1
      fi
      AUTO_CONFIRM="$1"
      ;;
    --max-missing-checks=*)
      MAX_MISSING_CHECKS="${1#*=}"
      ;;
    --max-missing-checks)
      shift
      if [[ "$#" -eq 0 ]]; then
        echo "--max-missing-checks requires a value" >&2
        exit 1
      fi
      MAX_MISSING_CHECKS="$1"
      ;;
    --)
      shift
      args+=("$@")
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      args+=("$1")
      ;;
  esac
  shift
done

set -- "${args[@]}"

if [[ "$#" -ne 2 ]]; then
  usage >&2
  exit 1
fi

ARRAY_ID="$1"
ADD_HOURS="$2"

require_positive_integer "ARRAY_ID" "${ARRAY_ID}"
require_positive_integer "ADD_HOURS" "${ADD_HOURS}"
require_positive_integer "CHECK_INTERVAL_SECONDS" "${CHECK_INTERVAL_SECONDS}"
require_positive_integer "THRESHOLD_HOURS" "${THRESHOLD_HOURS}"
require_positive_integer "MAX_MISSING_CHECKS" "${MAX_MISSING_CHECKS}"
require_walltime "${WALLTIME}"

export PATH="/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:${PATH}"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "Required command not found: sbatch" >&2
    exit 1
  fi

  echo "Submitting hourly add-hours monitor"
  echo "  array_id: ${ARRAY_ID}"
  echo "  add_hours: ${ADD_HOURS}"
  echo "  walltime: ${WALLTIME}"
  echo "  threshold_hours: ${THRESHOLD_HOURS}"
  echo "  check_interval_seconds: ${CHECK_INTERVAL_SECONDS}"
  echo "  max_missing_checks: ${MAX_MISSING_CHECKS}"
  echo "  dry_run: ${DRY_RUN}"
  echo "  auto_confirm: ${AUTO_CONFIRM}"

  exec sbatch \
    --parsable \
    --job-name=hourly_add_hours \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=1G \
    --time="${WALLTIME}" \
    --output="${RC_RESULTS:?}/logs/monitor_%j.out" \
    --error="${RC_RESULTS:?}/logs/monitor_%j.err" \
    --export="ALL,ADD_HOURS_SCRIPT=${ADD_HOURS_SCRIPT},CHECK_INTERVAL_SECONDS=${CHECK_INTERVAL_SECONDS},THRESHOLD_HOURS=${THRESHOLD_HOURS},DRY_RUN=${DRY_RUN},AUTO_CONFIRM=${AUTO_CONFIRM},WALLTIME=${WALLTIME},MAX_MISSING_CHECKS=${MAX_MISSING_CHECKS}" \
    "$0" "${ARRAY_ID}" "${ADD_HOURS}"
fi

for cmd in squeue bash date grep sleep; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Required command not found: ${cmd}" >&2
    exit 1
  fi
done

if [[ ! -f "${ADD_HOURS_SCRIPT}" ]]; then
  echo "Missing add-hours script: ${ADD_HOURS_SCRIPT}" >&2
  exit 1
fi

echo "Starting hourly add-hours monitor"
echo "  array_id: ${ARRAY_ID}"
echo "  add_hours: ${ADD_HOURS}"
echo "  add_hours_script: ${ADD_HOURS_SCRIPT}"
echo "  check_interval_seconds: ${CHECK_INTERVAL_SECONDS}"
echo "  max_missing_checks: ${MAX_MISSING_CHECKS}"
echo "  threshold_hours: ${THRESHOLD_HOURS}"
echo "  dry_run: ${DRY_RUN}"
echo "  auto_confirm: ${AUTO_CONFIRM}"
echo "  requested_walltime: ${WALLTIME}"
echo "  slurm_job_id: ${SLURM_JOB_ID:-NA}"

missing_checks=0
while true; do
  echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] Checking array ${ARRAY_ID}"

  queue_rows=""
  if queue_rows="$(squeue -h -r -j "${ARRAY_ID}" -o "%i" 2>/dev/null)" &&
      [[ -n "${queue_rows//[[:space:]]/}" ]]; then
    missing_checks=0
  else
    missing_checks=$((missing_checks + 1))
    echo "Array ${ARRAY_ID} not found in recursive squeue check ${missing_checks}/${MAX_MISSING_CHECKS}."
    if (( missing_checks >= MAX_MISSING_CHECKS )); then
      echo "Array ${ARRAY_ID} was missing for ${missing_checks} consecutive checks; exiting."
      exit 0
    fi
    echo "Sleeping ${CHECK_INTERVAL_SECONDS} seconds before rechecking."
    sleep "${CHECK_INTERVAL_SECONDS}"
    continue
  fi

  set +e
  THRESHOLD_HOURS="${THRESHOLD_HOURS}" \
    DRY_RUN="${DRY_RUN}" \
    AUTO_CONFIRM="${AUTO_CONFIRM}" \
    bash "${ADD_HOURS_SCRIPT}" "${ARRAY_ID}" "${ADD_HOURS}"
  status="$?"
  set -e

  if [[ "${status}" -ne 0 ]]; then
    echo "add-hours command exited with status ${status}; continuing after sleep." >&2
  fi

  echo "Sleeping ${CHECK_INTERVAL_SECONDS} seconds."
  sleep "${CHECK_INTERVAL_SECONDS}"
done
