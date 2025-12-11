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

int main(int argc, char *argv[]) {
    // we need an input file being specified
    if (argc == 1) {
        std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
        return EXIT_FAILURE;
    }
    const std::string filename = argv[1];
    const bool generateTrainingData = std::stoi(argv[2]) != 0;

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

    if (!generateTrainingData) {
        runSimulation(settings, "out");
    }
    else {
        // ToDo: This is nasty and should probably be configurable somehow
        const int N = 101;
        const double Re_min = 500.0;
        const double Re_max = 1500.0;
        const double u_min = 0.5;
        const double u_max = 1.5;

        for (int i = 0; i < N; ++i) {
            double t = static_cast<double>(i) / (N - 1);
            settings.re = Re_min + t * (Re_max - Re_min);
            settings.dirichletBcTop[0] = u_min + t * (u_max - u_min);
            std::ostringstream folder;
            folder << "train/out_" << std::setw(3) << std::setfill('0') << i;
            runSimulation(settings, folder.str());
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
