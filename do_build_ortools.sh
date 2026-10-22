#!/bin/bash
# Build, install and package OR-Tools (CP-SAT) as static libraries.
# Runs on Linux (inside the manylinux containers, see build_ortools.sh),
# on macOS, and on Windows (Visual Studio, dynamic runtime /MD, from Git Bash).
#
# Usage: do_build_ortools.sh [builddir [nthreads [packagename [cmake options...]]]]
#   example: ./do_build_ortools.sh build-x64 2 linux-x64
#   example (macOS universal): ./do_build_ortools.sh build 2 osx64 \
#                    -DCMAKE_OSX_ARCHITECTURES='arm64;x86_64'
#   example (Windows): ./do_build_ortools.sh build 2 win64
# The OR-Tools sources are expected in ./or-tools; the result is
# or-tools/<builddir>/ortools-<packagename>.zip
set -ex
cd `dirname $0`

source ./ortools_cmake_options.sh

builddirname="${1:-build}"
nthreads="${2:-2}"
packagename="${3:-package}"
if [ $# -ge 4 ]; then
    extraoptions="${@:4}"
else
    extraoptions=""
fi

if [ -z "${HOME:-}" ] || [ ! -w "$HOME" ]; then
  export HOME=/tmp/buildhome
  mkdir -p "$HOME"
fi

export TAR_OPTIONS="--no-same-owner"

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

mkdir -p $builddirname
installdir="$PWD/$builddirname/install"

# CMAKE_BUILD_TYPE is ignored by multi-config generators (the default one
# on Windows, Visual Studio): there the configuration comes from --config
cmake -S . -B $builddirname -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$installdir" \
      $ORTOOLS_CMAKE_OPTIONS $extraoptions
cmake --build $builddirname --config Release --parallel $nthreads
cmake --install $builddirname --config Release

# Package (cmake can create zip files on every platform)
cd "$installdir"
cmake -E tar cf "../ortools-$packagename.zip" --format=zip -- .
