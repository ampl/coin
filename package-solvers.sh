#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 build-dir output-dir platform [YYYYMMDD]" >&2
    exit 1
fi

SOURCE_DIR=$(cd "$(dirname "$0")" && pwd)
BUILD_DIR=$(cd "$1" && pwd)
PLATFORM=$3
PACKAGE_DATE=${4:-$(date +%Y%m%d)}
BIN_DIR=$BUILD_DIR/bin
case "$PLATFORM" in
    mswin*) BIN_DIR=$BIN_DIR/Release ;;
esac

cmake "-DPLATFORM=$PLATFORM" "-DBINDIR=$BIN_DIR" "-DOUTDIR=$2" \
    "-DSOURCE_DIR=$SOURCE_DIR" "-DVERSIONS_FILE=$BUILD_DIR/coin-versions.txt" \
    "-DPACKAGE_DATE=$PACKAGE_DATE" -P "$SOURCE_DIR/package/package-ampl.cmake"
