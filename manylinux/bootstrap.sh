#!/bin/bash

cd /tmp/
curl -O https://cmake.org/files/v3.12/cmake-3.12.4.tar.gz
tar xzvf cmake-3.12.4.tar.gz
cd cmake-3.12.4
./bootstrap
gmake
gmake install
rm -rf /tmp/cmake*