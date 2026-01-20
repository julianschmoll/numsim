#include "macros.h"
#include "settings.h"
#include "simulation/simulation.h"
#include "precice/precice.hpp"

#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sstream>

void runSimulation(const Settings &settings, const std::string &folderName) {
    Simulation simulation{settings, folderName};
    simulation.run();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int ownRankNo = 0;
    int nRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    // we need an input file being specified
    if (argc != 3) {
        std::cout << "The adapter was called with an incorrect number of arguments.\n";
        std::cout << "Usage: ./numsim_adapter [precice config] [scenario file]\n\n";
        std::cout << "Parameter description\n";
        std::cout << "    precice config:  Path and filename of preCICE configuration\n";
        std::cout << "    scenario file:   Path and filename of settings file for the solver\n";
        return EXIT_FAILURE;
    }
    const std::string preciceConfigPath(argv[1]);
    const std::string settingsPath(argv[2]);

    Settings settings;
    settings.loadFromFile(settingsPath);
    DEBUG(settings.printSettings());

    runSimulation(settings, "out");

    MPI_Finalize();
    return EXIT_SUCCESS;
}
