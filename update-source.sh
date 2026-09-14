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

IPOPT_VERSION=3.14.20
Mumps_VERSION=5.9.1

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

get_rel Ipopt
find ../Ipopt/Ipopt/* ! -iname "*cmake*" -delete
cp -r Ipopt-releases-$IPOPT_VERSION/* ../Ipopt/Ipopt/

# ../Ipopt/configure.ac (note: NOT ../Ipopt/Ipopt/configure.ac, which was
# just replaced above) is a separate, outer wrapper-project file that our
# CMakeLists.txt's configure() helper reads purely for the IPOPT_VERSION
# string it bakes into config_ipopt.h - keep its version in sync so that
# doesn't silently go stale (it did once already: an earlier version bump
# updated Ipopt/Ipopt/ but left this file at 3.12.13, so `ipopt -v` kept
# reporting the old version despite running genuinely newer code).
sed -i "s/AC_INIT(\[Ipopt\],\[[0-9.]*\]/AC_INIT([Ipopt],[$IPOPT_VERSION]/" ../Ipopt/configure.ac

# Ipopt's Mumps linear solver interface needs headers/macros (e.g. the
# MUMPS_VERSION macro, mumps_int_def.h) only present in newer MUMPS releases
# than the one previously vendored here (4.10.0), so it's fetched separately
# from the ThirdParty-Mumps mirror rather than via get_rel (which only knows
# the coin-or/* GitHub layout). PORD/METIS/SCOTCH ordering sources are
# intentionally not built (see Ipopt/CMakeLists.txt), so this is only used
# for its own AMD/QAMD/AMF orderings plus the numerical kernels.
curl -L http://coin-or-tools.github.io/ThirdParty-Mumps/MUMPS_$Mumps_VERSION.tar.gz > MUMPS_$Mumps_VERSION.tar.gz
tar xzvf MUMPS_$Mumps_VERSION.tar.gz
find ../ThirdParty/Mumps/MUMPS/* -delete
cp -r MUMPS_$Mumps_VERSION/* ../ThirdParty/Mumps/MUMPS/

get_rel Bonmin 
find ../Bonmin/* ! -iname "*cmake*" -delete
cp -r Bonmin-releases-$Bonmin_VERSION/* ../Bonmin/

# After updating couenne, in \Couenne\src\cut\sdpcuts\CouenneMatrix.hpp
# make function CouenneSparaseVector::compare_scalars::() const
#get_rel Couenne 
#find ../Couenne/* ! -iname "*cmake*" -delete
#cp -r Couenne-releases-$Couenne_VERSION/* ../Couenne/




