cmake_minimum_required(VERSION 3.8)
project(numsim)
add_subdirectory(src)

# Define NDEBUG if not building in Debug mode
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE NDEBUG)
endif ()
