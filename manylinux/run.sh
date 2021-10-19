#!/bin/bash
cd `dirname $0`

NPROC=`nproc`

rm -rf linux-/*

docker run -v `pwd`/../:/base --rm fdabrandao/manylinux32:coin_cmake-3.12.4 /base/manylinux/build.sh $NPROC intel32 32
docker run -v `pwd`/../:/base --rm fdabrandao/manylinux64:coin_cmake-3.12.4 /base/manylinux/build.sh $NPROC intel64 64

# docker run -v `pwd`/../:/base --rm dockcross/linux-armv7 /base/manylinux/build.sh $NPROC armv7 32
# docker run -v `pwd`/../:/base --rm fdabrandao/manylinux2014_aarch64:cmake_3.12.4 /base/manylinux/build.sh $NPROC arm64 arm64
# docker run -v `pwd`/../:/base --rm dockcross/manylinux2014-aarch64 /base/manylinux/build.sh $NPROC arm64 arm64
# docker run -v `pwd`/../:/base --rm messense/manylinux2014-cross:armv7l /base/manylinux/build.sh $NPROC arm64 arm64
# docker run -v `pwd`/../:/base --rm fdabrandao/manylinux2014_aarch64:cmake_3.12.4  /base/manylinux/build.sh $NPROC arm64 arm64


