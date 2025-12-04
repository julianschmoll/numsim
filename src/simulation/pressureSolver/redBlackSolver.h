#pragma once

#include "grid/staggeredGrid.h"
#include "pressureSolver.h"
#include "simulation/partitioning.h"
#include <memory>

class RedBlackSolver final : public PressureSolver {
public:
    RedBlackSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, double epsilon, int maximumNumberOfIterations,
                   double omega);

    void solve() override;

private:
    // Relaxation factor
    double omega_;

    double localPressureError_ = 0;
    double globalPressureError_ = 0;
};
