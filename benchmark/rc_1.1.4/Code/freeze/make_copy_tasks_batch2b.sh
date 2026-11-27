#!/usr/bin/env bash
# Batch 2b of the rc freeze: the code, libraries and inputs of the batch 2
# studies (C1-C7, C8a) that are not in their result folders.
# - Code/dev_repo: the karst_dev repository copy the dev studies ran from
#   (dev/, validation/, benchmark/, the package R/ and src/), as run;
# - Code/env: karst_cd342db (1.1.2) and karst_v4 (1.1.3 + dev patches R1, V4)
#   libraries and sources;
# - data/dev: the alfak2-era e1-e4 simulated inputs (DEV_BACKUP), the CE c2
#   grid and the dev/realdata/A Salehi table the dev studies read;
# - data/pancan: the NUMBAT-derived karyotype counts (ALFA-K inputs);
# - data/own_benchmark_dev: the development split's datasets and nodes used by C7.
#   bash make_copy_tasks_batch2b.sh            (on the HPC)
set -euo pipefail
P=/share/lab_crd/taoli/Project
RC=${RC_HOME:-$P/karst/benchmark/rc_1.1.4}
KD=$P/karst_dev; KDR=$KD/results
BK=$P/alfak2_results_backup_2026-11-17/validation/sparse_joint
PAN=$P/PANcanKFLs/Results/data
OWNDEV=$P/karst/benchmark/results/sparse_two_timepoint
T=$RC/Code/freeze/tasks; L=$T/lists_batch2b
mkdir -p "$L"
tf=$T/copy_batch2b.tsv
printf 'label\tsrc\tdst\tlist\n' > "$tf"
add() { local label=$1 src=$2 dst=$3; shift 3; printf '%s\n' "$@" > "$L/$label.txt"
        printf '%s\t%s\t%s\t%s\n' "$label" "$src" "$dst" "$L/$label.txt" >> "$tf"; }
add code_dev_repo $KD/repo $RC/Code/dev_repo dev validation benchmark R src DESCRIPTION NAMESPACE
add env_dev_libs $KD/lib $RC/Code/env/lib karst_cd342db karst_v4
add env_dev_src $KD/src $RC/Code/env/src karst_cd342db.tar.gz karst_v4.tar.gz
add data_dev_backup $BK $RC/data/dev/backup multi_timepoint/results/series/datasets multi_timepoint/results/series/nodes \
  series_methods/results/inputs/e2 series_methods/results/inputs/e3 series_methods/results/inputs/e4
add data_dev_ce $KDR/realdata/CE/c2 $RC/data/dev/realdata/CE/c2 grid
add data_dev_a $KDR/realdata/A $RC/data/dev/realdata/A salehi.rds
add data_pancan $PAN $RC/data/pancan $(cd $PAN && ls -d processed/ALFA_K/{GBM,HNSC,Ovarian,Melanoma}/{ALFA-K_inputs,karyotype_count})
ids=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "dataset_id") c = i; next} {print $c}' $P/karst_patch/marginal_pl/tasks/own_dev.tsv | sort -u)
add data_own_dev $OWNDEV $RC/data/own_benchmark_dev dataset_index_dev.csv $(for i in $ids; do echo datasets/$i.rds; echo nodes/$i.rds; done)
echo "$(( $(wc -l < "$tf") - 1 )) tasks in $tf; own dev datasets: $(echo "$ids" | wc -l)"
