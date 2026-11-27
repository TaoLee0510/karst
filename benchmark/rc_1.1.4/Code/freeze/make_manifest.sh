#!/usr/bin/env bash
# MANIFEST.tsv of the rc freeze: path (relative to RC_HOME), bytes, sha256 and
# whether the file is in the local copy, for every file under Code/, data/ and
# results_rc/. The sha256 of a copied file comes from its copy task
# (Code/freeze/manifest/sha256); files written or changed after the copies
# (all of Code/, results merged from the local machine, reports) are hashed
# here (sizes were checked by the copy tasks, Code/freeze/manifest/stats).
# Stops on a symbolic link or a file without a hash.
#   sbatch --cpus-per-task=8 ... --export=ALL,RC_HOME=... Code/freeze/make_manifest.sh
set -euo pipefail
RC=${RC_HOME:?}; cd "${RC}"
M=Code/freeze/manifest; W=$(mktemp -d "${TMPDIR:-/tmp}/manifest.XXXXXX"); trap 'rm -rf "${W}"' EXIT
find Code data results_rc -type l > "${W}/links"
[[ ! -s "${W}/links" ]] || { echo "symbolic links found:" >&2; cat "${W}/links" >&2; exit 1; }
find Code data results_rc -type f -printf '%p\t%s\n' | LC_ALL=C sort -t $'\t' -k1,1 > "${W}/files"
cat "${M}"/sha256/*.sha256 | awk '{h = $1; $1 = ""; sub(/^ +/, ""); print $0 "\t" h}' | LC_ALL=C sort -t $'\t' -k1,1 -u > "${W}/known"
# hash again: Code/, anything without a recorded hash, and the files replaced after the copy
{ awk -F'\t' '$1 ~ /^Code\//' "${W}/files" | cut -f1
  LC_ALL=C join -t $'\t' -v 1 "${W}/files" "${W}/known" | cut -f1
  printf '%s\n' results_rc/part1_paper/round_6ae2a24/analysis/fig2c_two_timepoints/oracle_both/manifest.R
} | LC_ALL=C sort -u > "${W}/rehash"
tr '\n' '\0' < "${W}/rehash" | xargs -0 -r -P "${SLURM_CPUS_PER_TASK:-4}" -n 64 sha256sum \
  | awk '{h = $1; $1 = ""; sub(/^ +/, ""); print $0 "\t" h}' | LC_ALL=C sort -t $'\t' -k1,1 > "${W}/fresh"
LC_ALL=C join -t $'\t' -v 1 "${W}/known" "${W}/fresh" > "${W}/known_kept"
LC_ALL=C sort -t $'\t' -k1,1 -m "${W}/known_kept" "${W}/fresh" > "${W}/hashes"
LC_ALL=C join -t $'\t' -a 1 "${W}/files" "${W}/hashes" > "${W}/joined"
awk -F'\t' 'NF < 3 {print "no hash: " $1; bad = 1} END {exit bad}' "${W}/joined"
# local copy (option C): no fit folders or fit objects, logs, image, built libraries, reference data
awk -F'\t' -v OFS='\t' 'BEGIN {print "path", "bytes", "sha256", "local"}
  {p = $1; loc = "yes"
   if (p ~ /\/fits\// || p ~ /\/objects\// || p ~ /\/logs\// || p ~ /^Code\/env\/lib\// || p ~ /\.sif$/ || p ~ /^data\/reference\//) loc = "no"
   print p, $2, $3, loc}' "${W}/joined" > MANIFEST.tsv
awk -F'\t' 'NR > 1 {n[$4]++; b[$4] += $2} END {for (k in n) printf "%s: %d files, %.2f GB\n", k, n[k], b[k] / 1e9}' MANIFEST.tsv
echo "MANIFEST.tsv: $(( $(wc -l < MANIFEST.tsv) - 1 )) files ($(wc -l < "${W}/rehash") hashed here)"
