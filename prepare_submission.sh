#!/bin/bash

echo "Removing old submission.zip"
rm -f submission.zip

echo "Creating new submission.zip with CMakeLists.txt and src"
zip -r submission.zip CMakeLists.txt src > /dev/null

echo "Testing created submission.zip"
tmp_dir=$(mktemp -d)
cp submission.zip "$tmp_dir"
cp scenarios/lid_driven_cavity.txt "$tmp_dir"
cd "$tmp_dir" || exit 1

echo "Unzipping in tmp dir"
unzip -q submission.zip

echo "Building..."
mkdir -p build && cd build || exit 1
cmake -DCMAKE_BUILD_TYPE=Release .. > /dev/null 2>&1
make install > /dev/null 2>&1

echo "Running with lid_driven_cavity scenario"
if ./numsim ../lid_driven_cavity.txt > /dev/null 2>&1; then
    echo "submission.zip is valid!"
else
    echo "submission.zip failed the test."
fi
