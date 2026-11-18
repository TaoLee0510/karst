#!/usr/bin/env bash
# Link alfak's fits of an earlier round into this round's results, so that a
# new karst commit is compared with the same alfak fits without refitting
# them (alfak, its settings and the inputs are unchanged). Each
# fits/<source>/<dataset>/<design>/alfak folder of REFERENCE becomes a
# symbolic link in BENCH_RESULTS; the karst folders are not touched. The
# linked folders are only read (scoring, collection); run_task.R is run for
# the karst tasks only.
#
#   hpc/link_alfak.sh REFERENCE_RESULTS [METHOD]   e.g. "${ALFAK_REFERENCE}"
#
# METHOD (default alfak) links another method's folders instead, e.g.
# `hpc/link_alfak.sh results_05a8855 karst` for a round that refits only
# karst's candidate missegregation rates (ac36fd1) and keeps 05a8855's karst
# fits (5e-5, unchanged code path) as its fixed-rate arm. Only real folders
# are linked (find -type d): link from the round that fitted them.
#
# Run from the repository root. Refuses to replace anything that exists;
# with SKIP_EXISTING=1 it leaves existing folders as they are and links the
# rest (e.g. 827c4b5: a refitted check subset of the pairwise and state-space
# fits, the rest linked from 05a8855 once the subset was found identical).
set -euo pipefail
cd "${BENCH_DIR:-$PWD}"
source hpc/env.sh
reference="$(realpath "$1")"
method="${2:-alfak}"
# Name of the method's folders in REFERENCE: the alfak2-era rounds named
# alfak's folders after the package, alfakR.
source_name="${SOURCE_NAME:-$([[ "${method}" == alfak ]] && echo "alfak""R" || echo "${method}")}"
target="${BENCH_RESULTS}"
[[ "${reference}" != "$(realpath -m "${target}")" ]] || { echo "REFERENCE is BENCH_RESULTS" >&2; exit 1; }
linked=0
skipped=0
while IFS= read -r dir; do
  rel="${dir#"${reference}"/}"
  dest="${target}/$(dirname "${rel}")/${method}"
  if [[ -e "${dest}" || -L "${dest}" ]]; then
    [[ "${SKIP_EXISTING:-0}" == 1 ]] || { echo "Exists: ${dest}" >&2; exit 1; }
    skipped=$((skipped + 1))
    continue
  fi
  mkdir -p "$(dirname "${dest}")"
  ln -s "${dir}" "${dest}"
  linked=$((linked + 1))
done < <(find "${reference}/fits" -mindepth 4 -maxdepth 4 -type d -name "${source_name}" | sort)
echo "Linked ${linked} ${method} fit folders from ${reference} into ${target} (existing, left as they are: ${skipped})"
