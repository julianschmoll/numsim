cmake_minimum_required(VERSION 3.18)
project(numsim_parallel)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -funroll-loops -flto -ffinite-math-only -fprefetch-loop-arrays -fno-trapping-math -ftree-vectorizer-verbose=2 -freciprocal-math -fstack-arrays -mavx512f -faggressive-loop-optimizations")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
endif()

add_subdirectory(src)
