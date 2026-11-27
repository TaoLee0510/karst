#!/usr/bin/env bash
# Install a karst source archive into lib/NAME with the image's R:
#   sbatch [resources] dev/hpc/build.sh NAME     (needs src/NAME.tar.gz)
# As benchmark/paper/hpc/build_karst.sh: install on node-local disk, then copy.
set -euo pipefail
source "${RC_HOME:?}/Code/dev_repo/dev/hpc/env.sh"   # rc freeze: RC_HOME from Code/rc_env.sh (exported)
load_hpc_environment
name="$1"
src="${RC_CODE}/env/src/${name}.tar.gz"
lib="${RC_CODE}/env/lib/${name}"
[[ -s "${src}" ]] || { echo "Missing ${src}" >&2; exit 1; }
[[ ! -e "${lib}/karst" ]] || { echo "Refusing to overwrite ${lib}/karst" >&2; exit 1; }
build="$(mktemp -d "${TMPDIR:-/tmp}/karst_dev_build.XXXXXX")"
mkdir -p "${lib}" "${build}/source"
tar -xzf "${src}" -C "${build}/source"
printf '%s\n' 'CXXFLAGS = -O2 -ftree-vectorize -march=x86-64 -fno-math-errno' \
  'CXX17FLAGS = -O2 -ftree-vectorize -march=x86-64 -fno-math-errno' >"${build}/Makevars"
local_dir="$(mktemp -d "${TMPDIR:-/tmp}/karst_dev_lib.XXXXXX")"
trap 'rm -rf "${local_dir}" "${build}"' EXIT
apptainer exec --cleanenv --bind "${RC_CODE}/env:${RC_CODE}/env" --bind "${build}:${build}" --env "LIB=${lib}" --env "BUILD=${build}" \
  --env "LOCAL=${local_dir}" "${DEV_SIF}" /bin/bash -c '
  set -euo pipefail
  source /usr/local/lib/pancankfls-runtime-env.sh
  export TMPDIR="${LOCAL}" R_MAKEVARS_USER="${BUILD}/Makevars"
  mkdir -p "${LOCAL}/library"
  /app/eb/software/R/4.4.2-gfbf-2024a/bin/R CMD INSTALL --library="${LOCAL}/library" "${BUILD}/source"
  cp -a "${LOCAL}/library/karst" "${LIB}/"
  export R_LIBS_USER="${LIB}:${R_LIBS_USER:-}"
  /app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript --vanilla -e "stopifnot(normalizePath(find.package(\"karst\")) == normalizePath(file.path(Sys.getenv(\"LIB\"), \"karst\"))); cat(\"karst\", as.character(packageVersion(\"karst\")), \"\n\")"
'
echo "BUILD_STATUS PASS ${lib}"
