#include "macros.h"
#include "settings.h"
#include "simulation/simulation.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <mpi.h>

void runSimulation(const Settings &settings, const std::string &folderName) {
    Simulation simulation{settings, folderName};
    simulation.run();
}

void generateTrainingData(Settings &settings) {
    std::cout << "Generating training data..." << std::endl;
    for (int i = 0; i < settings.nSamples; ++i) {
        const double interpolate = static_cast<double>(i) / (settings.nSamples - 1);
        settings.re = settings.reynoldsRange[0] + interpolate * (settings.reynoldsRange[1] - settings.reynoldsRange[0]);
        settings.dirichletBcTop[0] = settings.velocityRange[0] + interpolate * (settings.velocityRange[1] - settings.velocityRange[0]);
        // help, I hate handling strings in cpp
        std::ostringstream index;
        index << std::setw(4) << std::setfill('0') << i;
        std::string folder = "train/out_" + index.str();
        runSimulation(settings, folder);
    }
}

int main(int argc, char *argv[]) {
    // we need an input file being specified
    if (argc != 2) {
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

#ifdef USE_CG
    settings.pressureSolver = IterSolverType::CG;
#endif

    if (settings.generateTrainingData) {
        generateTrainingData(settings);
    } else {
        runSimulation(settings, "out");
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
