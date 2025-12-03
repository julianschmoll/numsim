#include "settings.h"
#include "simulation/simulation.h"
#include "macros.h"

#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    // we need an input file being specified
    if (argc == 1) {
        std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
        return EXIT_FAILURE;
    }
    const std::string filename = argv[1];

    Settings settings;
    settings.loadFromFile(filename);

    DEBUG(settings.printSettings());

    MPI_Init(&argc, &argv);

    int ownRankNo = 0;
    int nRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    {
        Simulation simulation{settings};
        simulation.run();
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
