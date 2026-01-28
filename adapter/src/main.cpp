#include "macros.h"
#include "precice/precice.hpp"
#include "settings.h"
#include "simulation/partitioning.h"
#include "simulation/simulation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>


void printVector1D(const std::string &label, const std::vector<double> &v, int maxElements = 10) {
    std::cout << "[adapter-debug] " << label << " size=" << v.size() << " {";
    int n = std::min<int>(static_cast<int>(v.size()), maxElements);
    for (int i = 0; i < n; ++i) {
        std::cout << v[i];
        if (i + 1 < n)
            std::cout << ", ";
    }
    if (v.size() > static_cast<size_t>(maxElements))
        std::cout << ", ...";
    std::cout << "}\n";
}

void setDisplacements(Simulation *simulation, std::vector<double> &displacements) {
    std::cout << "[adapter-debug] Setting Displacements\n";
    // TODO: Set Displacements in simulation from vector
}

void getForces(Simulation *simulation, std::vector<double> &forces) {
    std::cout << "[adapter-debug] Getting Forces\n";
    // TODO: Get Forces from simulation and write in vector
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
    precice::string_view displacementDelta = "DisplacementDelta";
    precice::string_view force = "Force";

    // Scope is here so MPI mem types are destroyed before MPI_Finalize
    {
        precice::Participant participant("Fluid", preciceConfigPath, ownRankNo, nRanks);
        Simulation simulation{settings, "/home/juschli/git/numsim_new/numsim/out"};
        auto partitioning = simulation.getPartitioning();
        auto diskOps = simulation.getDiscreteOperators();
        auto physicalSize = settings.physicalSize;
        int meshDim = participant.getMeshDimensions(fluidMeshNodes);

        DEBUG(std::cout << "Mesh dimension: " << meshDim << "\n");

        // +1 because of vertices, not faces
        auto verticesWidth = partitioning->nCellsLocal()[0] + 1;
        auto facesWidth = partitioning->nCellsLocal()[0] + 1;
        std::cout << "[adapter-debug] " << "Vertices width: " << verticesWidth << "\n";
        std::cout << "[adapter-debug] " << "Faces width: " << facesWidth << "\n";
        auto vertexSize = 2 * verticesWidth; // number of vertices at wet surface
        auto facesSize = 2 * facesWidth;

        DEBUG(std::cout << "Vertex size: " << vertexSize << "\n");
        DEBUG(std::cout << "Mesh size: " << facesSize << "\n");

        std::vector<double> nodeCoords(vertexSize * meshDim);
        std::vector<precice::VertexID> nodeIDs(vertexSize);

        std::vector<double> faceCoords(facesSize * meshDim);
        std::vector<precice::VertexID> faceIDs(facesSize);

        std::cout << "Phyiscal Size:" << physicalSize[0] << ", " << physicalSize[1] << "\n";

        int idx = 0;
        for (int row = 0; row < 2; ++row) {
            double y_pos = (row == 0) ? 0.0 : static_cast<double>(physicalSize[1]);
            for (int i = 0; i < verticesWidth; ++i) {
                double x_pos =
                    (verticesWidth == 1) ? 0.0 : (static_cast<double>(i) / static_cast<double>(verticesWidth - 1)) * static_cast<double>(physicalSize[0]);
                for (int d = 0; d < meshDim; ++d) {
                    double v = (d == 0) ? x_pos : (d == 1) ? y_pos : 0.0;
                    nodeCoords[idx++] = v;
                }
            }
        }

        idx = 0;
        for (int row = 0; row < 2; ++row) {
            double y_pos = (row == 0) ? 0.0 : static_cast<double>(physicalSize[1]);
            for (int i = 0; i < facesWidth; ++i) {
                double x_pos =
                    (facesWidth == 1) ? 0.0 : (static_cast<double>(i) / static_cast<double>(facesWidth - 1)) * static_cast<double>(physicalSize[0]);
                for (int d = 0; d < meshDim; ++d) {
                    double v = (d == 0) ? x_pos : (d == 1) ? y_pos : 0.0;
                    faceCoords[idx++] = v;
                }
            }
        }

        participant.setMeshVertices(fluidMeshNodes, nodeCoords, nodeIDs);
        participant.setMeshVertices(fluidMeshFaces, faceCoords, faceIDs);

        int displacementsDim = participant.getDataDimensions(fluidMeshNodes, displacementDelta);
        std::vector displacements(vertexSize * displacementsDim, 0.0);

        int forcesDim = participant.getDataDimensions(fluidMeshFaces, force);
        std::vector forces(facesSize * forcesDim, 0.0);

        DEBUG(std::cout << forces.size() << " force entries\n");

        participant.initialize();

        int meshSize = participant.getMeshVertexSize("Solid-Mesh");
        int dim = participant.getMeshDimensions("Solid-Mesh");

        std::vector<double> coords(meshSize * dim);
        std::vector<precice::VertexID> ids(meshSize);

        participant.getMeshVertexIDsAndCoordinates("Solid-Mesh", ids, coords);

        double currentTime = 0.0;

        std::cout << "[adapter-debug] " << "Starting coupling loop\n";

        printVector1D("nodeCoords", nodeCoords, 50);
        printVector1D("faceCoords", faceCoords, 50);

        std::cout << "preCICE expects " << participant.getMeshVertexSize(fluidMeshNodes) << " node vertices, adapter created " << vertexSize << "\n";

        std::cout << "preCICE expects " << participant.getMeshVertexSize(fluidMeshFaces) << " face vertices, adapter created " << facesSize << "\n";

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

            setDisplacements(&simulation, displacements);

            participant.startProfilingSection("Fluid Solver Step");
            simulation.advanceFluidSolver(dt);
            getForces(&simulation, forces);
            participant.stopLastProfilingSection();

            printVector1D("Forces (before write)", forces, 44);

            participant.writeData(fluidMeshFaces, force, faceIDs, forces);

            participant.advance(dt);

            if (participant.requiresReadingCheckpoint()) {
                simulation.reloadLastState();
            } else {
                simulation.updateSolid();
                currentTime += dt;
                int currentSec = static_cast<int>(currentTime);
                int lastSec = static_cast<int>(currentTime - dt);
                simulation.writeOutput(currentTime, currentSec, lastSec);
            }
        }

        participant.finalize();
        DEBUG(std::cout << "Finalizing\n");
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}
