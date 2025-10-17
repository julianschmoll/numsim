cmake_minimum_required(VERSION 3.8) # 3.8 as in https://numsim-exercises.readthedocs.io

project(numsim)

add_executable(numsim src/main.cpp)
target_include_directories(numsim PRIVATE include)

# required for second submission:
# find_package(MPI REQUIRED)
# include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(VTK 8 COMPONENTS CommonCore)
if (VTK_FOUND)
    message("Found VTK ${VTK_VERSION} in \"${VTK_DIR}\".")
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    add_definitions(-DVTK_AVAILABLE) # add_compile_definitions(VTK_AVAILABLE) # CMake >= 3.12
else()
    message("No VTK version >= 8 found.")
endif(VTK_FOUND)

add_compile_options(-Wall -Wextra)

set(CMAKE_CXX_STANDARD 14) # 14 as in https://numsim-exercises.readthedocs.io
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
