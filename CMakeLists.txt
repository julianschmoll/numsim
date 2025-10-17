cmake_minimum_required(VERSION 3.10)

project(numsim)


find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(Eigen3 3.4 REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_PATH})

# find_package(OpenMP REQUIRED)

add_executable(numsim src/main.cpp)
target_include_directories(numsim PRIVATE include)

set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
