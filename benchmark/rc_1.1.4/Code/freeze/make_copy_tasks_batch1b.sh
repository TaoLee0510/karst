#!/usr/bin/env bash
# Batch 1b of the rc freeze: what make_copy_tasks_batch1.sh left out. The built
# karst and extra R libraries (Code/env/lib), the staged 1.1.4 work tree that
# A3 / A4 ran with (karst_wip: library and source), and the dispersion table
# the rate-grid comparison reads (from the fig2d_misseg study).
#   bash make_copy_tasks_batch1b.sh            (on the HPC)
set -euo pipefail
P=/share/lab_crd/taoli/Project
RC=${RC_HOME:-$P/karst/benchmark/rc_1.1.4}
R6=$P/karst/benchmark/paper/results_6ae2a24
PAPER=$P/karst/benchmark/paper
T=$RC/Code/freeze/tasks; L=$T/lists_batch1b
mkdir -p "$L"
tf=$T/copy_batch1b.tsv
printf 'label\tsrc\tdst\tlist\n' > "$tf"
add() { local label=$1 src=$2 dst=$3; shift 3; printf '%s\n' "$@" > "$L/$label.txt"
        printf '%s\t%s\t%s\t%s\n' "$label" "$src" "$dst" "$L/$label.txt" >> "$tf"; }
add env_libs $PAPER/lib $RC/Code/env/lib karst_50e5889 karst_6ae2a24 extra_r
add env_wip_lib $R6/xv_target_code/lib $RC/Code/env/lib karst_wip
add env_wip_src $R6/xv_target_code $RC/Code/env/src karst_src
add a2_dispersion $R6/fig2d_misseg $RC/results_rc/part1_paper/salehi_rate_grid dispersion.tsv
echo "$(( $(wc -l < "$tf") - 1 )) tasks in $tf"
