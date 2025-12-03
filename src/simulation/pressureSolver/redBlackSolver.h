#pragma once

#include "grid/staggeredGrid.h"
#include "simulation/partitioning.h"
#include "pressureSolver.h"
#include <memory>

class RedBlackSolver : public PressureSolver {
public:
    RedBlackSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning,
        double epsilon, int maximumNumberOfIterations, double omega);

    void solve() final;

private:

    // Relaxation factor
    double omega_;

    double localPressureError_ = 0;
    double globalPressureError_ = 0;
};
