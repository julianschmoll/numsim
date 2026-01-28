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

    // Create a local Partitioning to compute offsets / local cell counts
    Partitioning partitioning(settings.nCells);

    // Create the simulation (it will create its own Partitioning internally)
    Simulation simulation(settings, "out");

    std::string meshNameNodes = "Fluid-Mesh-Nodes";
    std::string meshNameFaces = "Fluid-Mesh-Faces";
    std::string readDataName = "DisplacementDelta";
    std::string writeDataName = "Force";

    precice::Participant participant("Fluid", preciceConfigPath, ownRankNo, nRanks);

    const int nCellsXGlobal = settings.nCells[0];
    const int nCellsYGlobal = settings.nCells[1];
    const double domainWidth = settings.physicalSize[0];
    const double domainHeight = settings.physicalSize[1];
    const double hx = domainWidth / static_cast<double>(nCellsXGlobal);

    // local number of boundary points in x-direction for this rank
    const int nLocal = partitioning.nCellsLocal()[0];
    const int offsetX = partitioning.nodeOffset()[0];

    // We have two boundary lines: bottom and top. For each local partition we
    // expose the local subset of vertices on bottom and top.
    const int numLocalVerticesPerLine = nLocal;
    const int totalLocalVertices = 2 * numLocalVerticesPerLine; // bottom + top

    // Build vertices: [bottom(0..nLocal-1) | top(0..nLocal-1)] with (x,y) coordinates.
    std::vector<double> vertices;
    vertices.reserve(static_cast<size_t>(totalLocalVertices) * 2);

    std::vector<int> vertexIdsNodes;
    vertexIdsNodes.reserve(static_cast<size_t>(totalLocalVertices));

    // bottom vertices (y = 0)
    for (int i = 0; i < numLocalVerticesPerLine; ++i) {
        int globalI = offsetX + i;
        double x = (globalI + 0.5) * hx; // cell-centered
        double y = 0.0;
        vertices.push_back(x);
        vertices.push_back(y);
        // global unique id for bottom line (use 0..nCellsXGlobal-1)
        vertexIdsNodes.push_back(globalI);
    }

    // top vertices (y = domainHeight)
    for (int i = 0; i < numLocalVerticesPerLine; ++i) {
        int globalI = offsetX + i;
        double x = (globalI + 0.5) * hx;
        double y = domainHeight;
        vertices.push_back(x);
        vertices.push_back(y);
        // global unique id for top line (we offset by nCellsXGlobal so top ids are distinct)
        vertexIdsNodes.push_back(globalI + nCellsXGlobal);
    }

    // For this adapter we use the same vertex set for reading displacements and writing forces.
    // setMeshVertices expects the local subset of the global mesh and global vertex IDs.
    participant.setMeshVertices(meshNameNodes, vertices, vertexIdsNodes);
    participant.setMeshVertices(meshNameFaces, vertices, vertexIdsNodes);

    participant.initialize();

    // Each vertex has 2 components (x,y). We have totalLocalVertices vertices.
    std::vector<double> receiveBuffer(static_cast<size_t>(totalLocalVertices) * 2, 0.0);
    std::vector<double> sendBuffer(static_cast<size_t>(totalLocalVertices) * 2, 0.0);

    std::vector<double> topDisplacements(numLocalVerticesPerLine);
    std::vector<double> bottomDisplacements(numLocalVerticesPerLine);
    std::vector<double> topForces(numLocalVerticesPerLine, 0.0);
    std::vector<double> bottomForces(numLocalVerticesPerLine, 0.0);

    double currentTime = 0.0;

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            simulation.saveState();
        }

        double preciceDt = participant.getMaxTimeStepSize();
        TimeSteppingInfo tsInfo = simulation.computeTimeStepWidth(currentTime);
        double dt = std::min(preciceDt, tsInfo.timeStepWidth);
        simulation.setTimeStepWidth(dt);

        // read vector data: ordering is [vx, vy] per vertex -> repeated for all local vertices
        participant.readData(meshNameNodes, readDataName, vertexIdsNodes, dt, receiveBuffer);

        // Unpack: receiveBuffer stores [bottom0.x, bottom0.y, bottom1.x, bottom1.y, ..., top0.x, top0.y, ...]
        // We only care about vertical (y) displacement here.
        for (int i = 0; i < numLocalVerticesPerLine; ++i) {
            bottomDisplacements[i] = receiveBuffer[2 * i + 1];                           // bottom i -> y
            topDisplacements[i] = receiveBuffer[2 * (numLocalVerticesPerLine + i) + 1];  // top i -> y
        }

        // provide the solver with the displacements (top, bottom)
        simulation.setDisplacements(topDisplacements, bottomDisplacements);

        participant.startProfilingSection("Fluid Solver Step");

        simulation.advanceFluidSolver(dt);

        // TODO: Simulation must expose a getter to read computed forces (topForces / bottomForces).
        // For now the vectors remain zero unless Simulation::getForces(...) is implemented.
        // Example intended API:
        // simulation.getForces(topForces, bottomForces);

        participant.stopLastProfilingSection();

        // Pack vector data for preCICE: each vertex expects [fx, fy]. We send only vertical force in fy.
        for (int i = 0; i < numLocalVerticesPerLine; ++i) {
            // bottom vertex i
            sendBuffer[2 * i + 0] = 0.0;
            sendBuffer[2 * i + 1] = bottomForces[i];
            // top vertex i (located at index numLocalVerticesPerLine + i)
            sendBuffer[2 * (numLocalVerticesPerLine + i) + 0] = 0.0;
            sendBuffer[2 * (numLocalVerticesPerLine + i) + 1] = topForces[i];
        }

        participant.writeData(meshNameFaces, writeDataName, vertexIdsNodes, sendBuffer);

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

    MPI_Finalize();
    return EXIT_SUCCESS;
}

// Note: leave getForces/setDisplacements helper functions removed from this file;
// prefer implementing accessors in Simulation to obtain forces and set displacements cleanly.
