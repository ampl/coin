#!/bin/bash
# Combine the arm64 and x86_64 macOS builds (see macos-stage.sh) into one
# directory of universal binaries and static libraries, laid out like a build
# directory so package-solvers.sh can package and sign it.
#
#   macos-lipo.sh arm64-dir x86_64-dir output-dir
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 arm64-dir x86_64-dir output-dir" >&2
    exit 1
fi

ARM_DIR=$(cd "$1" && pwd)
INTEL_DIR=$(cd "$2" && pwd)
OUT_DIR=$3

merge() {
    # $1: path relative to both input directories (and the output directory)
    if [ ! -f "$INTEL_DIR/$1" ]; then
        echo "Missing $INTEL_DIR/$1" >&2
        exit 1
    fi
    lipo -create "$ARM_DIR/$1" "$INTEL_DIR/$1" -output "$OUT_DIR/$1"
    archs=$(lipo -archs "$OUT_DIR/$1")
    for arch in arm64 x86_64; do
        case " $archs " in
            *" $arch "*) ;;
            *) echo "$1 is not universal: $archs" >&2; exit 1 ;;
        esac
    done
    echo "$1: $archs"
}

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/bin"
for path in "$ARM_DIR"/bin/* "$ARM_DIR"/lib*.a; do
    merge "${path#"$ARM_DIR"/}"
done
# lipo keeps neither mode bits nor the build metadata; restore what's needed.
chmod 755 "$OUT_DIR"/bin/*
cp -p "$ARM_DIR/coin-versions.txt" "$ARM_DIR/coin-solvers.txt" "$OUT_DIR/"
# The Cbc libraries' generated config headers don't depend on the
# architecture; check that, and take one copy.
if ! diff -r "$ARM_DIR/cbc-libs-include" "$INTEL_DIR/cbc-libs-include"; then
    echo "The arm64 and x86_64 cbc-libs config headers differ" >&2
    exit 1
fi
cp -Rp "$ARM_DIR/cbc-libs-include" "$OUT_DIR/"
