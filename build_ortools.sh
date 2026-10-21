#!/bin/bash
# Linux builds in the manylinux containers (old glibc, for compatibility)
# Usage: build_ortools.sh [x64|arm64]
set -ex
cd `dirname $0`

toexecute=do_build_ortools.sh
if [ "$1" = "x64" ]; then
    image="quay.io/pypa/manylinux_2_28_x86_64"
    builddir="build-x64"
    nthreads=2
    packagename="linux-x64"
    user=""
elif [ "$1" = "arm64" ]; then
    image="quay.io/pypa/manylinux_2_28_aarch64"
    builddir="build-arm64"
    nthreads=4
    packagename="linux-arm64"
    user="-u $(id -u):$(id -g)"
else
    echo "Usage: $0 [x64|arm64]"
    exit 1
fi

docker run --rm $image /bin/bash -c "ldd --version"
docker run $user -v `pwd`:/ortools-build $image /bin/bash -c "chmod +x /ortools-build/$toexecute; /ortools-build/$toexecute $builddir $nthreads $packagename"
