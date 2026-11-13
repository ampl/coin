#!/bin/bash
set -ex
cd `dirname $0`

if [ $# -eq 0 ]; then
    builddirname="build"
    nthreads=2
else
    builddirname="$1"
    nthreads=$2
fi

if [ $# -ge 3 ]; then
    preprocdef="${@:3}"
else
    preprocdef=""
fi

if [ -z "${HOME:-}" ] || [ ! -w "$HOME" ]; then
  export HOME=/tmp/buildhome
  mkdir -p "$HOME"
fi

export TAR_OPTIONS="--no-same-owner"

git config --global --add safe.directory '*'

cd HiGHS
mkdir -p $builddirname
cd $builddirname


export PATH=/usr/local/cuda/bin:$PATH
# The CUDA images also carry GCC 13, because CUDA 12.9 can't use the GCC 14
# of the base image (see docker/dockerfile).
if [ -f /opt/rh/gcc-toolset-13/enable ]; then
  source /opt/rh/gcc-toolset-13/enable
fi
cmake .. -DBUILD_EXAMPLES=off -DZLIB=off -DCPack_CMake_INCLUDED=0 -DBUILD_SHARED_LIBS=ON  $preprocdef
cmake --build . --config Release -j$nthreads
cpack -G ZIP