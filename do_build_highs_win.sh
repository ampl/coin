#!/bin/bash
cd `dirname $0`

if [ $# -eq 0 ]; then
    builddirname="build"
    nthreads=2
else
    builddirname="$1"
    nthreads=$2
fi

if [ $# -eq 4 ]; then
    preprocdef="$3"
else
    preprocdef=""
fi

git config --global --add safe.directory `dirname $0`/HiGHS
cd HiGHS
vcpkg new --application
vcpkg add port metis
vcpkg add port openblas
mkdir $builddirname
cd $builddirname

cmake .. -DBUILD_EXAMPLES=off  -DZLIB=off $preprocdef -DBUILD_OPENBLAS=on -DHIPO=on -DCPack_CMake_INCLUDED=0
cmake --build . --config Release -j$nthreads
cpack -G ZIP