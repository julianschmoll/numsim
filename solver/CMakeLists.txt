cmake_minimum_required(VERSION 3.18)
project(numsim_parallel)
add_subdirectory(src)

# Define NDEBUG if not building in Debug mode
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE NDEBUG)
    #target_compile_definitions(${PROJECT_NAME} PRIVATE USE_CG)
    # Those flags are really aggressive but seem to work just fine (we might lose some accuracy though)
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -funroll-loops -flto -ffinite-math-only -fprefetch-loop-arrays -fno-trapping-math -ftree-vectorizer-verbose=2 -freciprocal-math -fstack-arrays -mavx512f -faggressive-loop-optimizations")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
endif()
