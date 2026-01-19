#!/bin/bash

rm -rf ./solver/build/
mkdir -p ./solver/build && cd ./solver/build
cmake ..
make install
mpirun ./numsim_parallel ../../cfg/scenarios/lid_driven_cavity.txt
