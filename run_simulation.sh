#!/bin/bash

rm -rf build
mkdir -p build && cd build
cmake ..
make install
mpirun ./numsim_parallel ..cfg/scenarios/lid_driven_cavity.txt
