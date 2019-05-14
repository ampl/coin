#!/bin/bash
# Usage example: ./try.sh -b coin-trusty64 -b coin-win64
ssh -M -S my-ctrl-socket -fnNT -L 4444:localhost:5555 ampl.com
ssh -S my-ctrl-socket -O check ampl.com
buildbot try --connect=pb --master=127.0.0.1:4444 --username=coin --passwd=coin --vc=git $@
ssh -S my-ctrl-socket -O exit ampl.com
