#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 nproc arch [32|64]"
elif [ -d "/base/manylinux/" ]; then
    NPROC=$1
    ARCH=$2
    NBITS=$3
    # The manylinux images only ship the C/C++ toolchain needed to build
    # Python extensions; gfortran (needed for Ipopt/Mumps/Lapack) and cmake
    # aren't preinstalled. Mumps 5.x's Fortran needs a reasonably modern
    # gfortran (older ones mishandle c_loc() on derived-type array
    # components), so prefer the newest available gcc-toolset (AlmaLinux's
    # devtoolset-style newer-GCC packages) over the base distro gcc-gfortran,
    # which on EL8-based images is a much older GCC 8.
    dnf install -y cmake
    newest_toolset=$(dnf list --available 'gcc-toolset-*-gcc-gfortran' 2>/dev/null \
      | awk '/^gcc-toolset-/{print $1}' | sort -t- -k3 -n -r | head -1)
    if [ -n "$newest_toolset" ]; then
        toolset_name=${newest_toolset%%-gcc-gfortran*}
        # Install the matching C and C++ compiler subpackages too, not just
        # -gcc-gfortran: otherwise plain-C code (e.g. ThirdParty/asl) can
        # silently fall back to the image's own newer default gcc, producing
        # object files whose LTO bytecode version doesn't match what this
        # toolset's linker expects ("bytecode stream ... generated with LTO
        # version X instead of the expected Y").
        dnf install -y "$toolset_name-gcc" "$toolset_name-gcc-c++" "$newest_toolset"
        source /opt/rh/$toolset_name/enable
    else
        dnf install -y gcc gcc-c++ gcc-gfortran
    fi
    gfortran --version
    BUILD_DIR=/base/manylinux/linux-$ARCH/
    rm -rf $BUILD_DIR
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR
    cp -r /base /tmp/coin
    # ThirdParty/asl enables LTO by default (USE_LTO option, ON for Release);
    # turn it off - we don't need the optimization, and it's exactly this
    # kind of cross-toolchain fragility (LTO bytecode is only compatible
    # between matching compiler versions) that's easiest to just not depend on
    # in a build that mixes gfortran/gcc/g++ across separately-installed
    # packages.
    cmake /tmp/coin/ -DARCH=$NBITS -DUSE_LTO=OFF # -DNO_AVX_HARDWARE=1
    make all coin-versions -j$NPROC
    make test || true
    make package
else
    echo "Must be run inside a docker container."
fi
