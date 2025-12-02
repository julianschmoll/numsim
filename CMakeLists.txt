cmake_minimum_required(VERSION 3.18)
project(numsim_parallel)
add_subdirectory(src)

# Define NDEBUG if not building in Debug mode
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE NDEBUG)
    set(CMAKE_CXX_FLAGS "-O3")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
endif()
