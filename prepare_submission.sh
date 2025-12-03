#!/bin/bash

echo "Removing old submission.zip"
rm -f submission.zip

echo "Creating new submission.zip with CMakeLists.txt and src"
zip -r submission.zip CMakeLists.txt src > /dev/null

echo "Created submission.zip"
