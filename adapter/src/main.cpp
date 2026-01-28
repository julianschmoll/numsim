// cpp
#include "macros.h"
#include "settings.h"
#include "simulation/simulation.h"
#include "simulation/partitioning.h"
#include "precice/precice.hpp"

#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <algorithm>


void printVector1D(const std::string &label,
                   const std::vector<double> &v,
                   int maxElements = 10) {
    std::cout << "[adapter-debug] " << label << " size=" << v.size() << " {";
    int n = std::min<int>(static_cast<int>(v.size()), maxElements);
    for (int i = 0; i < n; ++i) {
        std::cout << v[i];
        if (i + 1 < n) std::cout << ", ";
    }
    if (v.size() > static_cast<size_t>(maxElements)) std::cout << ", ...";
    std::cout << "}\n";
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int ownRankNo = 0;
    int nRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    if (argc != 3) {
        std::cout << "The adapter was called with an incorrect number of arguments.\n";
        std::cout << "Usage: ./numsim_adapter [precice config] [scenario file]\n\n";
        return EXIT_FAILURE;
    }
    const std::string preciceConfigPath(argv[1]);
    const std::string settingsPath(argv[2]);

    Settings settings;
    settings.loadFromFile(settingsPath);
    DEBUG(settings.printSettings());

    precice::string_view fluidMeshNodes = "Fluid-Mesh-Nodes";
    precice::string_view fluidMeshFaces = "Fluid-Mesh-Faces";
    precice::string_view displacementDelta = "Displacement";
    precice::string_view force = "Force";

    precice::Participant participant("Fluid", preciceConfigPath, ownRankNo, nRanks);

    Simulation simulation{settings, "out"};
    auto partitioning = simulation.getPartitioning();
    auto diskOps = simulation.getDiscreteOperators();
    auto physicalSize = settings.physicalSize;
    int meshDim = participant.getMeshDimensions(fluidMeshNodes);

    DEBUG(std::cout << "Mesh dimension: " << meshDim << "\n");

    // +1 because of vertices, not faces
    auto verticesWidth = partitioning->nCellsLocal()[0] + 1;
    auto facesWidth = partitioning->nCellsLocal()[0] +1;
    std::cout << "[adapter-debug] " << "Vertices width: " << verticesWidth << "\n";
    std::cout << "[adapter-debug] " << "Faces width: " << facesWidth << "\n";
    auto vertexSize = 2 * verticesWidth; // number of vertices at wet surface
    auto facesSize = 2 * facesWidth;

    DEBUG(std::cout << "Vertex size: " << vertexSize << "\n");
    DEBUG(std::cout << "Mesh size: " << facesSize << "\n");

    std::vector<double> nodeCoords(vertexSize*meshDim);
    std::vector<precice::VertexID> nodeIDs(vertexSize);

    std::vector<double> faceCoords(facesSize*meshDim);
    std::vector<precice::VertexID> faceIDs(facesSize);

    participant.setMeshVertices(fluidMeshNodes, nodeCoords, nodeIDs);
    participant.setMeshVertices(fluidMeshFaces, faceCoords, faceIDs);

    int displacementsDim = participant.getDataDimensions(fluidMeshNodes, displacementDelta);
    std::vector<double> displacements(vertexSize*displacementsDim, 0.0);

    int forcesDim = participant.getDataDimensions(fluidMeshFaces, force);
    std::vector<double> forces(facesSize*forcesDim, 0.0);

    DEBUG(std::cout << forces.size() << " force entries\n");

    participant.initialize();

    double currentTime = 0.0;

    std::cout << "[adapter-debug] " << "Starting coupling loop\n";

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            simulation.saveState();
        }

        std::cout << "[adapter-debug] " << "Coupling loop\n";

        double preciceDt = participant.getMaxTimeStepSize();
        TimeSteppingInfo tsInfo = simulation.computeTimeStepWidth(currentTime);
        double dt = std::min(preciceDt, tsInfo.timeStepWidth);
        simulation.setTimeStepWidth(dt);

        participant.readData(fluidMeshNodes, displacementDelta, nodeIDs, dt, displacements);

        printVector1D("Displacements (read)", displacements, 44);

        // setDisplacements(topDisplacements, bottomDisplacements);

        participant.startProfilingSection("Fluid Solver Step");
        simulation.advanceFluidSolver(dt);
        // getForces();
        participant.stopLastProfilingSection();

        printVector1D("Forces (before write)", forces, 44);

        participant.writeData(fluidMeshFaces, force, faceIDs, forces);

        participant.advance(dt);

        if (participant.requiresReadingCheckpoint()) {
            simulation.reloadLastState();
        } else {
            simulation.updateSolid();
            currentTime += dt;
            // int currentSec = static_cast<int>(currentTime);
            // int lastSec = static_cast<int>(currentTime - dt);
            // simulation.writeOutput(currentTime, currentSec, lastSec);
        }
    }

    participant.finalize();
    DEBUG(std::cout << "Finalizing\n");

    MPI_Finalize();
    return EXIT_SUCCESS;
}

void getForces() {

}

void setDisplacements() {

}
