#!/bin/bash

cd /tmp/
curl -O https://cmake.org/files/v3.5/cmake-3.5.0.tar.gz
tar xzvf cmake-3.5.0.tar.gz
cd cmake-3.5.0
./bootstrap
gmake
gmake install
