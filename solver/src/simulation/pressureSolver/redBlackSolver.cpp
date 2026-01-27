#include "redBlackSolver.h"
#include "grid/dataField.h"
#include "macros.h"
#include "simulation/partitioning.h"
#include <limits>

RedBlackSolver::RedBlackSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings)
    : PressureSolver(std::move(grid), std::move(partitioning), settings) {}

void RedBlackSolver::solve(DataField &p) {
    int it = 0;

    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1 / dx2;
    const double invDy2 = 1 / dy2;
    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

    globalPressureError_ = std::numeric_limits<double>::max();

    const int beginI = p.beginI() + 1;
    const int endI = p.endI() - 1;

    const int beginJ = p.beginJ() + 1;
    const int endJ = p.endJ() - 1;

    const double epsilonSquared = settings_.epsilon * settings_.epsilon;
    const double oneMinusOmega = 1 - settings_.omega;
    const double omegaTimesScalingFactor = settings_.omega * scalingFactor;

    auto [i0, j0] = partitioning_->nodeOffset();
    const int globalOffset = (i0 + j0) % 2;

    while (globalPressureError_ > epsilonSquared && it < settings_.maximumNumberOfIterations) {
        localPressureError_ = 0.0;
        // rb = 0 is red, rb=1 is black pass
        for (int rb = 0; rb < 2; ++rb) {
#pragma omp parallel for reduction(+ : localPressureError_)
            for (int j = beginJ; j < endJ; ++j) {
                // offset is 0 or 1 depending on the color of the cell
                const int offset = (beginI + rb + globalOffset + j) & 1;
#pragma omp simd reduction(+ : localPressureError_)
                for (int i = beginI + offset; i < endI; i += 2) {
                    if (grid_->isSolid(i, j)) continue;
                    const double pDxx = (p(i - 1, j) + p(i + 1, j)) * invDx2;
                    const double pDyy = (p(i, j - 1) + p(i, j + 1)) * invDy2;
                    const double pNew = oneMinusOmega * p(i, j) + omegaTimesScalingFactor * (pDxx + pDyy - rhs(i, j));
                    localPressureError_ += (p(i, j) - pNew) * (p(i, j) - pNew);
                    p(i, j) = pNew;
                }
            }
            partitioning_->nonBlockingExchange(p);
            setBoundaryValues(p);
            setStructureBoundaries(p);
            partitioning_->waitForAllMPIRequests();
        }
        MPI_Allreduce(&localPressureError_, &globalPressureError_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        it++;
    }
}
