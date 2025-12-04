#pragma once

#include "grid/staggeredGrid.h"
#include "simulation/partitioning.h"
#include <memory>

class PressureSolver {

public:
    PressureSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, double epsilon, int maximumNumberOfIterations);

    virtual ~PressureSolver() = default;

    PressureSolver(const PressureSolver &) = delete;

    PressureSolver(PressureSolver &&) = delete;

    virtual void solve() = 0;

protected:
    void setBoundaryValues();

    std::shared_ptr<Partitioning> partitioning_;

    // object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;

    // Convergence threshold for the residual
    double epsilon_;

    // Maximum number of iterations allowed in the solver
    double maxNumberOfIterations_;
};
