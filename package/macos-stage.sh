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
# Static libraries and the files package-solvers.sh reads from the build root.
cp -p "$BUILD_DIR"/lib*.a "$OUT_DIR/"
cp -p "$BUILD_DIR/coin-versions.txt" "$BUILD_DIR/coin-solvers.txt" "$OUT_DIR/"
