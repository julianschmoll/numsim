cmake_minimum_required(VERSION 3.8)
project(numsim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME}
        # main stuff
        settings.cpp
        main.cpp

        #grid stuff
        grid/array2d.cpp
        grid/staggeredGrid.cpp
        grid/dataField.cpp

        # solver stuff
        simulation/simulation.cpp
        simulation/pressureSolver.cpp
        simulation/discreteOperators.cpp

        # outputWriter stuff
        outputWriter/outputWriter.cpp
        outputWriter/outputWriterParaview.cpp
        outputWriter/outputWriterText.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)

if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    message(STATUS "Found VTK ${VTK_VERSION} in \"${VTK_DIR}\"")
else ()
    message(FATAL_ERROR "VTK not found. Please install VTK and set VTK_DIR if needed.")
endif ()

# install numsim executable in build directory and get rid of annoying warnings
# make_policy(SET CMP0177 NEW # TODO: does not compile!! Added in version 3.31.
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
