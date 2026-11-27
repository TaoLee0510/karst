#!/usr/bin/env bash
# Install CRAN packages the image lacks into lib/extra_r, with the R and
# compilers of the image (as hpc/build_karst.sh). Run as a Slurm job from
# the repository root:  sbatch [resources] hpc/build_extra_r.sh
# Installs every source tarball in src/extra/ (downloaded on the login node,
# e.g. glmmTMB for the paper's Fig. 4 models) in the order given by
# src/extra/ORDER, and refuses to overwrite an installed package.
set -euo pipefail
cd "${BENCH_DIR:-$PWD}"
source hpc/env.sh
load_hpc_environment
# rc freeze: the source tarballs are no longer kept (only the built library
# Code/env/lib/extra_r); put them in Code/env/src/extra to rebuild.
SRC="${RC_CODE}/env/src/extra"
[[ -s "${SRC}/ORDER" ]] || { echo "Missing ${SRC}/ORDER" >&2; exit 1; }
mkdir -p "${EXTRA_LIB}"
while read -r tarball; do
  pkg="${tarball%%_*}"
  [[ ! -e "${EXTRA_LIB}/${pkg}" ]] || { echo "Refusing to overwrite ${EXTRA_LIB}/${pkg}" >&2; exit 1; }
done <"${SRC}/ORDER"
# Inside the image R's file.access() reports the NFS share as not writable:
# install on the node's local disk, then copy.
LOCAL="$(mktemp -d "${TMPDIR:-/tmp}/extra_r.XXXXXX")"
trap 'rm -rf "${LOCAL}"' EXIT
apptainer exec --cleanenv --bind "${RC_CODE}/env:${RC_CODE}/env" \
  --env "EXTRA_LIB=${EXTRA_LIB}" --env "SRC=${SRC}" --env "LOCAL=${LOCAL}" "${BENCH_SIF}" /bin/bash -c '
  set -euo pipefail
  source /usr/local/lib/pancankfls-runtime-env.sh
  export TMPDIR="${LOCAL}"
  mkdir -p "${LOCAL}/library"
  export R_LIBS_USER="${R_LIBS_USER:-}:${LOCAL}/library"
  while read -r tarball; do
    /app/eb/software/R/4.4.2-gfbf-2024a/bin/R CMD INSTALL --library="${LOCAL}/library" "${SRC}/${tarball}"
  done <"${SRC}/ORDER"
  cp -a "${LOCAL}/library/." "${EXTRA_LIB}/"
  export R_LIBS_USER="${R_LIBS_USER}:${EXTRA_LIB}"
  /app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript --vanilla -e "
    for (t in readLines(file.path(Sys.getenv(\"SRC\"), \"ORDER\"))) {
      p <- sub(\"_.*$\", \"\", t)
      suppressPackageStartupMessages(library(p, character.only = TRUE))
      cat(\"EXTRA_R\t\", p, \"\t\", as.character(packageVersion(p)), \"\t\", find.package(p), \"\n\", sep = \"\")
    }"
'
echo "EXTRA_R_INSTALL_STATUS	PASS	${EXTRA_LIB}"
