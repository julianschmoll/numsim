cmake_minimum_required(VERSION 3.8)
project(numsim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME}
        main.cpp
        output_writer/write_paraview_output.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/output_writer)

find_package(VTK REQUIRED)

if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    message(STATUS "Found VTK ${VTK_VERSION} in \"${VTK_DIR}\"")
else()
    message(FATAL_ERROR "VTK not found. Please install VTK and set VTK_DIR if needed.")
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
