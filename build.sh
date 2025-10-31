#!/bin/bash

rm -rf build
mkdir -p build && cd build
cmake ..
make install
./numsim ../scenarios/lid_driven_cavity.txt
