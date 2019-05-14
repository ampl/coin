#!/bin/bash
# Usage example: ./try-local.sh -b coin-trusty64 -b coin-win64
buildbot try --connect=pb --master=127.0.0.1:5555 --username=coin --passwd=coin --vc=git $@

