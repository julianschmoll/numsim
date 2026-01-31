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
    std::cout << "[adapter-debug] " << label << " vertices=" << ids.size() << std::endl;

    if (ids.empty())
        return;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < ids.size(); ++i) {
        double x = coords[2 * i + 0];
        double y = coords[2 * i + 1];
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
        double x = coords[2 * i + 0];
        double y = coords[2 * i + 1];

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

void printVector(const std::string &label, const std::vector<double> &v, int maxElements = std::numeric_limits<int>::max()) {
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

    // Scope is here so MPI mem types are destroyed before MPI_Finalize
    {
        precice::string_view fluidMeshNodes = "Fluid-Mesh-Nodes";
        precice::string_view fluidMeshFaces = "Fluid-Mesh-Faces";
        precice::string_view displacementDelta = "DisplacementDelta";
        precice::string_view force = "Force";

        precice::Participant participant("Fluid", preciceConfigPath, ownRankNo, nRanks);

        // ToDo: Pass out folder maybe
        Simulation simulation{settings, "out"};
        auto partitioning = simulation.getPartitioning();
        auto diskOps = simulation.getDiscreteOperators();
        auto physicalSize = settings.physicalSize;
        int meshDim = participant.getMeshDimensions(fluidMeshNodes);

        //assert(meshDim == 2 && "Adapter only supports 2D meshes");

        DEBUG(std::cout << "Mesh dimension: " << meshDim << "\n");

        // +1 because of vertices, not faces
        auto verticesWidth = partitioning->nCellsLocal()[0];
        auto facesWidth = partitioning->nCellsLocal()[0];
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

        // Both Meshes have the same layout
        double dx = static_cast<double>(physicalSize[0]) / partitioning->nCellsLocal()[0];
        
        int idx = 0;
        for (int row = 0; row < 2; ++row) {
            double y_pos = (row == 0) ? 0.0 : static_cast<double>(physicalSize[1]);
            for (int i = 0; i < verticesWidth; ++i) {
                double x_pos = (static_cast<double>(i) + 0.5) * dx;
                for (int d = 0; d < meshDim; ++d) {
                    double v = (d == 0) ? x_pos : (d == 1) ? y_pos : 0.0;
                    nodeCoords[idx] = v;
                    faceCoords[idx] = v;
                    idx++;
                }
            }
        }

        participant.setMeshVertices(fluidMeshNodes, nodeCoords, nodeIDs);
        participant.setMeshVertices(fluidMeshFaces, faceCoords, faceIDs);

        printVector("Fluid Mesh", nodeCoords);
        printVector("Fluid Faces", faceCoords);

        int displacementsDim = participant.getDataDimensions(fluidMeshNodes, displacementDelta);
        std::vector displacements(vertexSize * displacementsDim, 0.15);

        int forcesDim = participant.getDataDimensions(fluidMeshFaces, force);
        std::vector forces(facesSize * forcesDim, 0.15);

        DEBUG(std::cout << forces.size() << " force entries\n");

        participant.initialize();

        int meshSize = participant.getMeshVertexSize("Solid-Mesh");
        int dim = participant.getMeshDimensions("Solid-Mesh");
        std::cout << "Participant Mesh Dimensions: size=" << meshSize << " dim=" << dim << "\n";

        std::vector<double> coords(meshSize * dim);
        std::vector<precice::VertexID> ids(meshSize);

        participant.getMeshVertexIDsAndCoordinates("Solid-Mesh", ids, coords);
        //printMesh("Solid-Mesh", ids, coords);

        std::cout << "Solid-Mesh " << coords << std::endl;

        double currentTime = 0.0;

        std::cout << "[adapter-debug] " << "Starting coupling loop\n";

        printVector("nodeCoords", nodeCoords, 50);
        printVector("faceCoords", faceCoords, 50);

        DEBUG(std::cout << "preCICE expects " << participant.getMeshVertexSize(fluidMeshNodes) << " node vertices, adapter created " << vertexSize
                        << "\n");
        assert(participant.getMeshVertexSize(fluidMeshNodes) == vertexSize);
        DEBUG(std::cout << "preCICE expects " << participant.getMeshVertexSize(fluidMeshFaces) << " face vertices, adapter created " << facesSize
                        << "\n");
        assert(participant.getMeshVertexSize(fluidMeshFaces) == facesSize);

        // Set initial displacements
        std::cout << "[adapter-debug] " << "Setting initial top wall displacements to " << settings.topWallDispl_ << " (top) and "
                  << settings.bottomWallDispl_ << " (bottom).\n";

        // +1 again because we have a flat array with x,y,x,y,...
        const int bottomOffset = verticesWidth * meshDim + 1;
        for (int i = 0; i < verticesWidth; ++i) {
            int idxTop = 1 + i * meshDim;
            int idxBottom = bottomOffset + i * meshDim;

            if (idxTop >= 0 && idxTop < static_cast<int>(displacements.size()))
                displacements[idxTop] = settings.topWallDispl_;
            if (idxBottom >= 0 && idxBottom < static_cast<int>(displacements.size()))
                displacements[idxBottom] = -settings.bottomWallDispl_;
        }
        printVector("Displacements (initial)", displacements, 44);
        simulation.setDisplacements(displacements);
        simulation.saveState();
        participant.writeData(fluidMeshFaces, force, faceIDs, forces);

        while (participant.isCouplingOngoing()) {
            if (participant.requiresWritingCheckpoint()) {
                simulation.saveState();
            }

            double preciceDt = participant.getMaxTimeStepSize();
            auto tsInfo = simulation.computeTimeStepWidth();
            double dt = std::min(preciceDt, tsInfo.timeStepWidth);
            simulation.setTimeStepWidth(dt);

            std::cout << "[adapter-debug] " << "Coupling Timestep " << dt << "\n";

            participant.readData(fluidMeshNodes, displacementDelta, nodeIDs, dt, displacements);
            printVector("Displacements (read)", displacements, 44);

            double maxDispl = 0.0;
            for (size_t i = 0; i < displacements.size(); ++i) {
                maxDispl = std::max(maxDispl, std::abs(displacements[i]));
            }
            if (maxDispl > 100.0) {
                std::cerr << "ERROR: Received invalid displacement magnitude: " << maxDispl << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            simulation.setDisplacements(displacements);

            participant.startProfilingSection("Fluid Solver Step");
            simulation.advanceFluidSolver(dt);
            
            simulation.getForces(forces);
            participant.stopLastProfilingSection();

            std::cout << "[adapter-debug] participant.writeData()" << std::endl;

            participant.writeData(fluidMeshFaces, force, faceIDs, forces);

            participant.advance(dt);

            if (participant.requiresReadingCheckpoint()) {
                simulation.reloadLastState();
            } else {
                //simulation.updateSolid();
                currentTime += dt;
                int currentSec = static_cast<int>(currentTime);
                int lastSec = static_cast<int>(currentTime - dt);
                simulation.writeOutput(currentSec, lastSec);
                forces.assign(forces.size(), 0.0);
            }
        }

        participant.finalize();
        DEBUG(std::cout << "You are crazy! You did it! The simulation might have finished successfully! \n");
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}
