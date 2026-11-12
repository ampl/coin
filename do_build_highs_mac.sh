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

if [ $# -eq 3 ]; then
    preprocdef=$3
else
    preprocdef=""
fi

if [ -z "${HOME:-}" ] || [ ! -w "$HOME" ]; then
  export HOME=/tmp/buildhome
  mkdir -p "$HOME"
fi

export TAR_OPTIONS="--no-same-owner"

git config --global --add safe.directory `dirname $0`/HiGHS

cd HiGHS

mkdir -p $builddirname
cd $builddirname

cmake .. -DBUILD_EXAMPLES=off -DHIPO=on  -DZLIB=off -DCMAKE_RELEASE_POSTFIX="-ampl" -DCPack_CMake_INCLUDED=0 $preprocdef
cmake --build . --config Release -j$nthreads
cpack -G ZIP