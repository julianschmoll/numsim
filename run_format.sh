#!/bin/bash

find . -iname '*.h' -o -iname '*.cpp' | clang-format --style=file -i --files=/dev/stdin