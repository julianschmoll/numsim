#pragma once

#include "grid/staggeredGrid.h"
#include "simulation/partitioning.h"
#include <memory>

class RedBlackSolver {
public:
    RedBlackSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning,
        double epsilon, int maximumNumberOfIterations, double omega);

    void solve();

private:
    void setBoundaryValues();

    std::shared_ptr<Partitioning> partitioning_;

    // object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;

    // Convergence threshold for the residual
    double epsilon_;

    // Maximum number of iterations allowed in the solver
    double maxNumberOfIterations_;

    // Relaxation factor
    double omega_;

    double localPressureError_ = 0;
    double globalPressureError_ = 0;
};
