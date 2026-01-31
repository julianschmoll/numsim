cmake_minimum_required(VERSION 3.18)
project(numsim_parallel)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS_RELEASE )
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
endif()

add_subdirectory(src)
