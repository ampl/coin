#!/bin/bash
cd `dirname $0`

rm -rf linux32/*
rm -rf linux64/*
mkdir linux32 linux64
docker run -v `pwd`/../:/base --rm fdabrandao/manylinux32:coin.v1 /base/manylinux/build.sh 32
docker run -v `pwd`/../:/base --rm fdabrandao/manylinux64:coin.v1 /base/manylinux/build.sh 64
