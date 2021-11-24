#!/bin/bash
cd `dirname $0`
mkdir -p tmp
cd tmp
rm -rf *

CBC_VERSION=2.10.5
IPOPT_VERSION=3.12.13
BONMIN_VERSION=1.8.8
COUENNE_VERSION=0.5.8

curl -O https://www.coin-or.org/download/source/Cbc/Cbc-$CBC_VERSION.tgz
tar xzvf Cbc-$CBC_VERSION.tgz
find ../Cbc/* ! -iname "*cmake*" -delete
cp -r Cbc-$CBC_VERSION/* ../Cbc/

curl -O https://www.coin-or.org/download/source/Ipopt/Ipopt-$IPOPT_VERSION.tar.gz
tar xzvf Ipopt-$IPOPT_VERSION.tar.gz
find ../Ipopt/* ! -iname "*cmake*" -delete
cp -r Ipopt-releases-$IPOPT_VERSION/* ../Ipopt/

curl -O https://www.coin-or.org/download/source/Bonmin/Bonmin-$BONMIN_VERSION.tgz
tar xzvf Bonmin-$BONMIN_VERSION.tgz
find ../Bonmin/* ! -iname "*cmake*" -delete
cp -r Bonmin-$BONMIN_VERSION/* ../Bonmin/

curl -O https://www.coin-or.org/download/source/Couenne/Couenne-$COUENNE_VERSION.tgz
tar xzvf Couenne-$COUENNE_VERSION.tgz
find ../Couenne/* ! -iname "*cmake*" -delete
cp -r Couenne-$COUENNE_VERSION/* ../Couenne/


