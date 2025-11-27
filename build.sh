#!/bin/bash

rm -rf build
mkdir -p build && cd build
cmake ..
make install
mpirun ./numsim_parallel ../scenarios/lid_driven_cavity.txt
