#!/bin/bash
# Copy what the macOS packaging needs from one architecture's build directory,
# so it can be published as a pipeline artifact and combined with lipo later.
#
#   macos-stage.sh build-dir output-dir
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 build-dir output-dir" >&2
    exit 1
fi

BUILD_DIR=$(cd "$1" && pwd)
OUT_DIR=$2
SOLVERS=bonmin,couenne,ipopt
if [ -f "$BUILD_DIR/coin-solvers.txt" ]; then
    SOLVERS=$(cat "$BUILD_DIR/coin-solvers.txt")
fi

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/bin"
for solver in ${SOLVERS//,/ }; do
    cp -p "$BUILD_DIR/bin/$solver" "$OUT_DIR/bin/"
done

# The solvers must only depend on libraries every Mac has (not, for example,
# Homebrew's libgcc_s or libgfortran).
status=0
for solver in ${SOLVERS//,/ }; do
    extra=$(otool -L "$OUT_DIR/bin/$solver" | tail -n +2 | awk '{print $1}' \
        | grep -v -E '^(/usr/lib/|/System/)' || true)
    if [ -n "$extra" ]; then
        echo "$solver depends on non-system libraries:" >&2
        echo "$extra" >&2
        status=1
    fi
done
[ "$status" -eq 0 ] || exit 1
# Static libraries and the files package-solvers.sh reads from the build root.
cp -p "$BUILD_DIR"/lib*.a "$OUT_DIR/"
cp -p "$BUILD_DIR/coin-versions.txt" "$BUILD_DIR/coin-solvers.txt" "$OUT_DIR/"
