#!/bin/bash
cd `dirname $0`
mkdir -p tmp
cd tmp
rm -rf *

Cbc_VERSION=2.10.8
Clp_VERSION=1.17.7
CoinUtils_VERSION=2.11.6
Osi_VERSION=0.108.7
Cgl_VERSION=0.60.6

IPOPT_VERSION=3.12.13

Bonmin_VERSION=1.8.9

Couenne_VERSION=0.5.8

get_rel () {
    declare version=$1_VERSION
    fname=$1-${!version}.tar.gz
    curl -L https://github.com/coin-or/$1/archive/refs/tags/releases/${!version}.tar.gz > $fname
    tar xzvf $fname
}

copy_sublib () {
    get_rel $1
    declare version=$1_VERSION
    dirname=$1-releases-${!version}
    cp -r $dirname/$1 ../Cbc/
}
get_rel Cbc 
find ../Cbc/* ! -iname "*cmake*" -delete
cp -r Cbc-releases-$Cbc_VERSION/* ../Cbc/

copy_sublib Clp
copy_sublib CoinUtils
copy_sublib Osi
copy_sublib Cgl

# IPOPT needs substantial work
#curl -O https://www.coin-or.org/download/source/Ipopt/Ipopt-$IPOPT_VERSION.tar.gz
#tar xzvf Ipopt-$IPOPT_VERSION.tgz
#find ../Ipopt/* ! -iname "*cmake*" -delete
#cp -r Ipopt-$IPOPT_VERSION/* ../Ipopt/

get_rel Bonmin 
find ../Bonmin/* ! -iname "*cmake*" -delete
cp -r Bonmin-releases-$Bonmin_VERSION/* ../Bonmin/

# After updating couenne, in \Couenne\src\cut\sdpcuts\CouenneMatrix.hpp
# make function CouenneSparaseVector::compare_scalars::() const
#get_rel Couenne 
#find ../Couenne/* ! -iname "*cmake*" -delete
#cp -r Couenne-releases-$Couenne_VERSION/* ../Couenne/




