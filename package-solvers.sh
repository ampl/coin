#!/bin/bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 build-dir output-dir platform" >&2
    exit 1
fi

SOURCE_DIR=$(cd "$(dirname "$0")" && pwd)
BUILD_DIR=$(cd "$1" && pwd)
PLATFORM=$3
BIN_DIR=$BUILD_DIR/bin
case "$PLATFORM" in
    mswin*) BIN_DIR=$BIN_DIR/Release ;;
esac

cmake "-DPLATFORM=$PLATFORM" "-DBINDIR=$BIN_DIR" "-DOUTDIR=$2" \
    "-DSOURCE_DIR=$SOURCE_DIR" "-DVERSIONS_FILE=$BUILD_DIR/coin-versions.txt" \
    -P "$SOURCE_DIR/package/package-ampl.cmake"
