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


# Install gcc 13 because gcc 14 in incompatible with cuda 12.9
#dnf install -y gcc-toolset-13 gcc-toolset-13-gcc gcc-toolset-13-gcc-c++
#source /opt/rh/gcc-toolset-13/enable

cd HiGHS
export TAR_OPTIONS="--no-same-owner"
git config --global --add safe.directory '*'
export PATH=/usr/local/cuda/bin:$PATH

mkdir -p $builddirname
cd $builddirname



cmake .. -DBUILD_EXAMPLES=off -DZLIB=off -DCPack_CMake_INCLUDED=0 -DBUILD_SHARED_LIBS=ON  $preprocdef
cmake --build . --config Release -j$nthreads
cpack -G ZIP