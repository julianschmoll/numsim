#include "macros.h"
#include "settings.h"
#include "simulation/simulation.h"
#include "precice/precice.hpp"

#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sstream>


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

    Simulation simulation(settings, "out");
    precice::Participant participant("Fluid", preciceConfigPath, ownRankNo, nRanks);

    participant.initialize();
    // prepare step (set boundaries once I guess)

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            DEBUG(std::cout << "Writing iteration checkpoint\n");
            simulation.saveState();
        }

        // ToDo: Handle solver time step being smaller than precice time step
        double dt = participant.getMaxTimeStepSize();

        // Read Boundary Data and apply conditions
        // participant.readData(...)

        participant.startProfilingSection("Fluid Solver Step");
        // Solve Fluid and Calculate Forces
        participant.stopLastProfilingSection();

        // Write with write data
        // participant.writeData(...);

        // Advance time step
        participant.advance(dt);

        // Reload Checkpoint if convergence failed
        if (participant.requiresReadingCheckpoint()) {
            DEBUG(std::cout << "Reading iteration checkpoint\n");
            simulation.reloadLastState();
        } else {
            DEBUG(std::cout << "Advancing in time\n");
            simulation.setVelocities();
            // simulation.writeOutput(...);
        }
    }

    participant.finalize();
    DEBUG(std::cout << "Finalizing\n");

    MPI_Finalize();
    return EXIT_SUCCESS;
}
