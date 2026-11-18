#!/usr/bin/env bash
# Install karst into a side library with the R and compilers of the image,
# as PANcanKFLs installs alfak (code/container/scripts/
# build_alfakr_f93b2d2_package.sh). Run as a Slurm job from the repository
# root:  sbatch [resources] hpc/build_karst.sh
# Needs src/karst_<commit>.tar.gz (git archive of the commit in
# hpc/karst.lock) and refuses to overwrite an existing install.
set -euo pipefail
cd "${BENCH_DIR:-$PWD}"
source hpc/env.sh
load_hpc_environment

SRC_TAR="${BENCH_DIR}/src/karst_${KARST_COMMIT}.tar.gz"
BUILD="${BENCH_DIR}/build/karst_${KARST_COMMIT}"
[[ -s "${SRC_TAR}" ]] || { echo "Missing ${SRC_TAR}" >&2; exit 1; }
printf '%s  %s\n' "${KARST_SHA256}" "${SRC_TAR}" | sha256sum --check --strict
[[ ! -e "${KARST_LIB}/karst" ]] || { echo "Refusing to overwrite ${KARST_LIB}/karst" >&2; exit 1; }
mkdir -p "${KARST_LIB}" "${BUILD}/source" "${BUILD}/tmp"
tar -xzf "${SRC_TAR}" -C "${BUILD}/source"
grep -qx 'Package: karst' "${BUILD}/source/DESCRIPTION"
printf '%s\n' \
  'CXXFLAGS = -O2 -ftree-vectorize -march=x86-64 -fno-math-errno' \
  'CXX17FLAGS = -O2 -ftree-vectorize -march=x86-64 -fno-math-errno' >"${BUILD}/Makevars"

# Inside the image R's file.access() reports the NFS share as not writable,
# so R CMD INSTALL refuses it: install on the node's local disk, then copy.
LOCAL="$(mktemp -d "${TMPDIR:-/tmp}/karst_lib.XXXXXX")"
trap 'rm -rf "${LOCAL}"' EXIT
apptainer exec --cleanenv --bind "${BENCH_DIR}:${BENCH_DIR}" \
  --env "KARST_LIB=${KARST_LIB}" --env "BUILD=${BUILD}" --env "LOCAL=${LOCAL}" "${BENCH_SIF}" /bin/bash -c '
  set -euo pipefail
  source /usr/local/lib/pancankfls-runtime-env.sh
  export TMPDIR="${LOCAL}"
  export R_MAKEVARS_USER="${BUILD}/Makevars"
  mkdir -p "${LOCAL}/library"
  /app/eb/software/R/4.4.2-gfbf-2024a/bin/R CMD INSTALL --library="${LOCAL}/library" "${BUILD}/source"
  cp -a "${LOCAL}/library/karst" "${KARST_LIB}/"
  export R_LIBS_USER="${KARST_LIB}:${R_LIBS_USER:-}"
  /app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript --vanilla -e "
    stopifnot(normalizePath(find.package(\"karst\")) == normalizePath(file.path(Sys.getenv(\"KARST_LIB\"), \"karst\")))
    cat(\"KARST_VERSION\t\", as.character(packageVersion(\"karst\")), \"\n\", sep = \"\")
    cat(\"ALFAKR\t\", as.character(packageVersion(\"alfakR\")), \"\t\", find.package(\"alfakR\"), \"\n\", sep = \"\")
  "
  ldd "${KARST_LIB}/karst/libs/karst.so"
  if ldd "${KARST_LIB}/karst/libs/karst.so" | grep -q "not found"; then
    echo "Unresolved karst shared-library dependency" >&2
    exit 1
  fi
'
echo "KARST_INSTALL_STATUS	PASS	${KARST_LIB}"
