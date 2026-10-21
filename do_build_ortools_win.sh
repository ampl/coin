#!/bin/bash
# Build, install and package OR-Tools (CP-SAT) as static libraries on Windows
# (Visual Studio 2022, dynamic runtime /MD; run from Git Bash).
#
# Usage: do_build_ortools_win.sh [builddir [nthreads [packagename]]]
# The OR-Tools sources are expected in ./or-tools; the result is
# or-tools/<builddir>/ortools-<packagename>.zip
set -ex
cd `dirname $0`

source ./ortools_cmake_options.sh

builddirname="${1:-build}"
nthreads="${2:-2}"
packagename="${3:-win64}"

git config --global --add safe.directory '*'

cd or-tools

# OR-Tools builds zlib and bzip2 as DLLs even for a static build; make them
# static so that there are no DLLs for them to distribute
patchfile=../patches/ortools-static-zlib-bzip2.patch
if git apply --reverse --check --ignore-whitespace $patchfile 2>/dev/null; then
    echo "$patchfile already applied"
else
    git apply --ignore-whitespace $patchfile
fi

installdir="$PWD/$builddirname/install"

cmake -S . -B $builddirname -G "Visual Studio 17 2022" -A x64 \
      -DCMAKE_INSTALL_PREFIX="$installdir" \
      $ORTOOLS_CMAKE_OPTIONS
cmake --build $builddirname --config Release --target install --parallel $nthreads

# Package (cmake can create zip files on every platform)
cd "$installdir"
cmake -E tar cf "../ortools-$packagename.zip" --format=zip -- .
