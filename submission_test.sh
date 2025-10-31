#!/bin/bash

# unzip submission.zip                    # this unzips the uploaded archive, filename may be different
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make install
./numsim ../scenarios/lid_driven_cavity.txt
# instead of lid_driven_cavity.txt, other files will be tested and should produce the respective results
