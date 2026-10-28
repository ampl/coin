#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 build-dir output-dir platform [YYYYMMDD]" >&2
    exit 1
fi

SOURCE_DIR=$(cd "$(dirname "$0")" && pwd)
BUILD_DIR=$(cd "$1" && pwd)
mkdir -p "$2"
OUTPUT_DIR=$(cd "$2" && pwd)
PLATFORM=$3
PACKAGE_DATE=${4:-$(date +%Y%m%d)}

BIN_DIR=$BUILD_DIR/bin
EXE_SUFFIX=
case "$PLATFORM" in
    mswin*)
        BIN_DIR=$BIN_DIR/Release
        EXE_SUFFIX=.exe
        ;;
esac

STAGING_DIR=$(mktemp -d)
trap 'rm -rf "$STAGING_DIR"' EXIT

for solver in cbc bonmin couenne ipopt; do
    package_dir=$STAGING_DIR/$solver
    mkdir -p "$package_dir/docs"
    cp "$BIN_DIR/$solver$EXE_SUFFIX" "$package_dir/"
    cp "$SOURCE_DIR/Cbc/Cbc/LICENSE" "$package_dir/coin-license.txt"
    cp "$BUILD_DIR/coin-versions.txt" "$package_dir/"
    cp "$SOURCE_DIR/docs/CHANGES.$solver.md" "$package_dir/docs/"
    if [ -n "$EXE_SUFFIX" ] && [ "$solver" != cbc ]; then
        cp "$BIN_DIR/libipoptfort.dll" "$package_dir/"
    fi
    # Keep macOS resource-fork metadata out of the distributable archives.
    COPYFILE_DISABLE=1 tar -czf "$OUTPUT_DIR/solvers-$solver.$PLATFORM.$PACKAGE_DATE.tar.gz" \
        -C "$package_dir" .
done
