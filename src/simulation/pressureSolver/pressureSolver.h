#pragma once

#include <memory>
#include "grid/staggeredGrid.h"
#include <memory>

class PressureSolver {
public:
    virtual ~PressureSolver() = default;
    PressureSolver(std::shared_ptr<StaggeredGrid> grid, double epsilon, double maxNumberOfIterations, double omega);

    /**
     * solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
     */
    virtual void solve() = 0;

protected:
    // object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;

    // Convergence threshold for the residual
    double epsilon_;

    // Maximum number of iterations allowed in the solver
    double maxNumberOfIterations_;

    double omega_;
};
