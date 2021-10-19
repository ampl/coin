#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 nproc arch [32|64]"
elif [ -d "/base/manylinux/" ]; then
    NRPOC=$1
    ARCH=$2
    NBITS=$3
    BUILD_DIR=/base/manylinux/linux-$ARCH/
    rm -rf $BUILD_DIR
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR
    cp -r /base /tmp/coin
    cmake /tmp/coin/ -DARCH=$NBITS # -DNO_AVX_HARDWARE=1
    make all coin-versions -j$NPROC
    make test || true
    make package
else
    echo "Must be run inside a docker container."
fi
