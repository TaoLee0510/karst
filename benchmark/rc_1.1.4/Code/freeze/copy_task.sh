#!/usr/bin/env bash
# One copy task of the rc freeze (array task; TASK_FILE from make_copy_tasks_*.sh).
# Copies the list's entries SRC -> DST with rsync -aL (real files, symlinks
# followed), then writes, under $RC_HOME/Code/freeze/manifest/:
#   stats/<label>.tsv    label, rsync exit, source files, source bytes, copied files, copied bytes
#   sha256/<label>.sha256  sha256 of every copied file, paths relative to RC_HOME
#   errors/<label>.txt   rsync's error lines (only when rsync failed)
#SBATCH --cpus-per-task=4
set -uo pipefail
RC=${RC_HOME:?}; TF=${TASK_FILE:?}; id=${SLURM_ARRAY_TASK_ID:?}
line=$(awk -F'\t' -v i="$id" 'NR == i + 1' "$TF")
IFS=$'\t' read -r label src dst list <<< "$line"
M=$RC/Code/freeze/manifest; mkdir -p "$M/stats" "$M/sha256" "$M/errors" "$dst"
echo "Host $(hostname)  task ${id} ${label}  $(date --iso-8601=seconds)"
echo "  ${src} -> ${dst} ($(wc -l < "$list") entries)"
rsync -aL -r --files-from="$list" "$src/" "$dst/" 2> "$M/errors/${label}.txt"; status=$?
[[ ${status} -eq 0 ]] && rm -f "$M/errors/${label}.txt"
count() { (cd "$1" && tr '\n' '\0' < "$list" | xargs -0 -I{} find $2 {} -type f -printf '%s\n' 2>/dev/null) | awk '{n++; s += $1} END {printf "%d\t%.0f", n, s}'; }
printf '%s\t%s\t%s\t%s\n' "$label" "$status" "$(count "$src" -L)" "$(count "$dst" "")" > "$M/stats/${label}.tsv"
rel=${dst#"$RC"/}
(cd "$RC" && while IFS= read -r e; do find "$rel/$e" -type f -print0; done < "$list" \
   | xargs -0 -r -P "${SLURM_CPUS_PER_TASK:-4}" -n 64 sha256sum) > "$M/sha256/${label}.sha256"
cat "$M/stats/${label}.tsv"
echo "Finished $(date --iso-8601=seconds)"
