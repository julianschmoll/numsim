#include "macros.h"
#include "precice/precice.hpp"
#include "settings.h"
#include "simulation/partitioning.h"
#include "simulation/simulation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <vector>

void printMesh(const std::string &label, const std::vector<precice::VertexID> &ids, const std::vector<double> &coords, int printWidth = 24,
               int printHeight = 10) {
    std::cout << "- [adapter] " << label << " vertices=" << ids.size() << std::endl;

    if (ids.empty())
        return;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < ids.size(); ++i) {
        double x = coords[3 * i + 0];
        double y = coords[3 * i + 1];
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
    }

    if (xmax == xmin)
        xmax = xmin + 1.0;
    if (ymax == ymin)
        ymax = ymin + 1.0;

    std::vector grid(printHeight, std::vector<std::string>(printWidth, "   + "));
    std::vector count(printHeight, std::vector(printWidth, 0));

    for (size_t i = 0; i < ids.size(); ++i) {
        double x = coords[3 * i + 0];
        double y = coords[3 * i + 1];

        int gx = (int)((x - xmin) / (xmax - xmin) * (printWidth - 1));
        int gy = (int)((y - ymin) / (ymax - ymin) * (printHeight - 1));
        gy = printHeight - 1 - gy;

        if (gx < 0 || gx >= printWidth || gy < 0 || gy >= printHeight)
            continue;

        count[gy][gx]++;

        if (count[gy][gx] == 1) {
            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << ids[i];
            grid[gy][gx] = ss.str();
        } else {
            grid[gy][gx] = "MULTI";
        }
    }

    std::cout << "  Vertex layout:\n";
    std::cout << "  +" << std::string(printWidth * 7, '-') << "+\n";

    for (int j = 0; j < printHeight; ++j) {
        std::cout << "  |";
        for (int i = 0; i < printWidth; ++i) {
            std::string cell = grid[j][i];
            if (cell.size() < 5)
                cell = std::string(5 - cell.size(), ' ') + cell;
            std::cout << "-" << cell << " ";
        }
        std::cout << "|\n";
    }
    std::cout << "  +" << std::string(printWidth * 7, '-') << "+\n";
}

void run(int ownRankNo, int nRanks, const std::string preciceConfigPath, Settings settings) {
    const precice::string_view fluidMesh = "Fluid-Mesh";
    const precice::string_view solidMesh = "Solid-Mesh";
    const precice::string_view displacementDelta = "DisplacementDelta";
    const precice::string_view force = "Force";

    precice::Participant participant("Fluid", preciceConfigPath, ownRankNo, nRanks);

    Simulation simulation{settings, "out"};
    auto partitioning = simulation.getPartitioning();
    auto diskOps = simulation.getDiscreteOperators();
    auto physicalSize = settings.physicalSize;
    int meshDim = participant.getMeshDimensions(fluidMesh);

    assert((meshDim == 2 || meshDim == 3) && "Adapter only supports 2D or 3D meshes");

    DEBUG(std::cout << "Mesh dimension: " << meshDim << "\n");

    auto verticesWidth = partitioning->nCellsLocal()[0];

    // number of vertices at wet surface
    auto vertexCount = 2 * verticesWidth;

    DEBUG(std::cout << "Vertex count: " << vertexCount << "\n");

    std::vector<double> vertexCoordinates(vertexCount * meshDim);
    std::vector<precice::VertexID> vertexIDs(vertexCount);

    // Both Meshes have the same layout
    double dx = static_cast<double>(physicalSize[0]) / partitioning->nCellsLocal()[0];

    int idx = 0;
    for (int row = 0; row < 2; ++row) {
        double xPos = (row == 0) ? 0.0 : static_cast<double>(physicalSize[1]);
        for (int i = 0; i < verticesWidth; ++i) {
            double yPos = (static_cast<double>(i) + 0.5) * dx;
            for (int d = 0; d < meshDim; ++d) {
                double v = (d == 0) ? yPos : (d == 1) ? xPos : 0.0;
                vertexCoordinates[idx] = v;
                idx++;
            }
        }
    }

    participant.setMeshVertices(fluidMesh, vertexCoordinates, vertexIDs);

    int displacementsDim = participant.getDataDimensions(fluidMesh, displacementDelta);
    std::vector displacements(vertexCount * displacementsDim, 0.0);

    int forcesDim = participant.getDataDimensions(fluidMesh, force);
    std::vector forces(vertexCount * forcesDim, 0.0);

    DEBUG(std::cout << forces.size() << " force entries\n");

    participant.initialize();

    const int meshSize = participant.getMeshVertexSize(solidMesh);
    const int solidMeshDim = participant.getMeshDimensions(solidMesh);

    std::vector<double> participantCoords(meshSize * solidMeshDim);
    std::vector<precice::VertexID> participantIDs(meshSize);

    participant.getMeshVertexIDsAndCoordinates(solidMesh, participantIDs, participantCoords);
    DEBUG(std::cout << participantCoords << "\n");

    double currentTime = 0.0;

    DEBUG(std::cout << "Expected " << participant.getMeshVertexSize(fluidMesh) << " vertices, created " << vertexCount << ".\n");
    assert(participant.getMeshVertexSize(fluidMesh) == vertexCount);

    // Set initial displacements
    std::cout << "- [adapter] Setting initial displacements...\n";

    const int topOffset = static_cast<int>(participantCoords.size() / 2) + 1;
    double topWallDispl = participantCoords[1];
    double bottomWallDispl = settings.physicalSize[1] - participantCoords[topOffset];

    // +1 again because we have a flat array and write to y component only
    const int bottomOffset = verticesWidth * meshDim + 1;
    for (int i = 0; i < verticesWidth; ++i) {
        int idxTop = 1 + i * meshDim;
        int idxBottom = bottomOffset + i * meshDim;

        if (idxTop >= 0 && idxTop < static_cast<int>(displacements.size()))
            displacements[idxTop] = topWallDispl;
        if (idxBottom >= 0 && idxBottom < static_cast<int>(displacements.size()))
            displacements[idxBottom] = -bottomWallDispl;
    }

    simulation.initializeDisplacements(displacements);

    participant.writeData(fluidMesh, force, vertexIDs, forces);

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            simulation.saveState();
        }

        double preciceDt = participant.getMaxTimeStepSize();
        auto tsInfo = simulation.computeTimeStepWidth();
        double dt = std::min(preciceDt, tsInfo.timeStepWidth);
        simulation.setTimeStepWidth(dt);

        participant.readData(fluidMesh, displacementDelta, vertexIDs, dt, displacements);

        double maxDispl = 0.0;
        for (size_t i = 0; i < displacements.size(); ++i) {
            maxDispl = std::max(maxDispl, std::abs(displacements[i]));
        }
        if (maxDispl > 100.0) {
            std::cerr << "- [adapter] ERROR: Received invalid displacement magnitude: " << maxDispl << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        simulation.setDisplacements(displacements);

        participant.startProfilingSection("Fluid Solver Step");
        simulation.advanceFluidSolver(dt);

        simulation.getForces(forces);
        participant.stopLastProfilingSection();

        participant.writeData(fluidMesh, force, vertexIDs, forces);

        participant.advance(dt);

        if (participant.requiresReadingCheckpoint()) {
            simulation.reloadLastState();
        } else {
            simulation.updateSolid();
            currentTime += dt;
            int currentSec = static_cast<int>(currentTime);
            int lastSec = static_cast<int>(currentTime - dt);
            assert(simulation.getCurrentTime() == currentTime);
            simulation.writeOutput(currentSec, lastSec);
        }
    }

    participant.finalize();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int ownRankNo = 0;
    int nRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    if (argc != 3) {
        std::cout << "- [adapter] The adapter was called with an incorrect number of arguments.\n";
        std::cout << "      Usage: ./numsim_adapter [precice config] [scenario file]\n\n";
        return EXIT_FAILURE;
    }
    const std::string preciceConfigPath(argv[1]);
    const std::string settingsPath(argv[2]);

    Settings settings;
    settings.loadFromFile(settingsPath);
    DEBUG(settings.printSettings());

    run(ownRankNo, nRanks, preciceConfigPath, settings);

    std::cout << "- [adapter] You are crazy! You did it! The simulation might have finished successfully! \n";
    MPI_Finalize();
    return EXIT_SUCCESS;
}
