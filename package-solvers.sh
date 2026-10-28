#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 build-dir output-dir platform [--signed]" >&2
    exit 1
fi

SIGNING_IDENTITY=
if [ "$#" -eq 4 ]; then
    if [ "$4" != --signed ]; then
        echo "Unknown option: $4" >&2
        exit 1
    fi
    SIGNING_IDENTITY=${MACOS_SIGNING_IDENTITY:-Developer ID Application: AMPL Optimization Inc. (ZNNBG5892S)}
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
    "-DSIGNING_IDENTITY=$SIGNING_IDENTITY" "-DSIGNING_KEYCHAIN=${MACOS_SIGNING_KEYCHAIN:-}" \
    -P "$SOURCE_DIR/package/package-ampl.cmake"
