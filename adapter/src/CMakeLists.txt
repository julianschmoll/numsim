cmake_minimum_required(VERSION 3.18)
project(numsim_adapter)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../solver ${CMAKE_CURRENT_BINARY_DIR}/solver)

add_executable(${PROJECT_NAME}
        main.cpp
)

# Links to our solver library
target_link_libraries(${PROJECT_NAME} PRIVATE solver)

find_package(precice REQUIRED CONFIG)

target_link_libraries(${PROJECT_NAME} PRIVATE precice::precice)

# install adapter executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
