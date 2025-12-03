cmake_minimum_required(VERSION 3.18)
project(numsim_parallel)

set(CMAKE_CXX_STANDARD 17)
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
    simulation/pressureSolver/redBlackSolver.cpp
    simulation/pressureSolver/conjugateGradientSolver.cpp
    simulation/discreteOperators.cpp
    simulation/partitioning.cpp

    # outputWriter stuff
    outputWriter/outputWriter.cpp
    outputWriter/outputWriterParaview.cpp
    outputWriter/outputWriterText.cpp
    outputWriter/outputWriterTextParallel.cpp
    outputWriter/outputWriterParaviewParallel.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if (VTK_FOUND)
    include_directories(${VTK_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    message(STATUS "Found VTK ${VTK_VERSION} in \"${VTK_DIR}\"")
else()
    message(FATAL_ERROR "VTK not found. Please install VTK and set VTK_DIR if needed.")
endif()

# install numsim executable in build directory and get rid of annoying warnings
# make_policy(SET CMP0177 NEW # TODO: does not compile!! Added in version 3.31.
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)
