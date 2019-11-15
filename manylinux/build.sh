#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 [32|64]"
elif [ -d "/base/manylinux/" ]; then
    BUILD_DIR=/base/manylinux/linux$1/
    rm -rf $BUILD_DIR
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR
    cmake /base/
    make all coin-versions -j2
    make test || true
else
    echo "Must be run inside a docker container."
fi
