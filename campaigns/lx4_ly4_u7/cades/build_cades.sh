#!/usr/bin/env bash
set -eo pipefail

source /etc/profile
set -u

SOURCE_DIR=${1:?usage: build_cades.sh SOURCE_DIR BIN_DIR [intel|gcc]}
BIN_DIR=${2:?usage: build_cades.sh SOURCE_DIR BIN_DIR [intel|gcc]}
TOOLCHAIN=${3:-intel}
BUILD_ROOT=${BUILD_ROOT:-"$SOURCE_DIR/build-cades-$TOOLCHAIN"}

module purge
case "$TOOLCHAIN" in
  intel)
    module load intel/2024.1.0 mkl/2024.1.0 cmake/3.24.4
    CXX_COMPILER=$(command -v icpx)
    BLA_VENDOR_VALUE=Intel10_64lp_seq
    ;;
  gcc)
    module load gcc/12.2.0 mkl/2024.1.0 cmake/3.24.4
    CXX_COMPILER=$(command -v g++)
    BLA_VENDOR_VALUE=Intel10_64lp_seq
    ;;
  *)
    echo "Unknown toolchain: $TOOLCHAIN" >&2
    exit 2
    ;;
esac

cmake -S "$SOURCE_DIR" -B "$BUILD_ROOT" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$CXX_COMPILER" \
  -DBLA_VENDOR="$BLA_VENDOR_VALUE" \
  -DBUILD_TESTING=ON
cmake --build "$BUILD_ROOT" --parallel 8
ctest --test-dir "$BUILD_ROOT" --output-on-failure

mkdir -p "$BIN_DIR"
for executable in \
  ftlm_n_vs_mu \
  ftlm_reduce_checkpoint \
  ftlm_ed_n_vs_mu \
  ftlm_conductivity \
  ftlm_ed_conductivity; do
  cp -p "$BUILD_ROOT/$executable" "$BIN_DIR/$executable"
done
{
  echo "toolchain=$TOOLCHAIN"
  echo "compiler=$($CXX_COMPILER --version | head -1)"
  echo "mklroot=${MKLROOT:-}"
  echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  sha256sum "$BIN_DIR"/ftlm_*
} | tee "$BIN_DIR/BUILD_PROVENANCE.txt"
