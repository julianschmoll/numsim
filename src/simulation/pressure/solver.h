#pragma once

#include "grid/discretization.h"
#include <memory>

class solver {
    /**
     * Sets the boundary values to account for homogenous Neumann boundary conditions.
     *
     * This has to be called after every iteration.
     */
    void setBoundaryValues();

public:
    virtual ~solver() = default;
    solver(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations);

    /**
     * Solves the Poisson problem for the pressure.
     *
     * Uses the rhs and p field variables in the grid.
     */
    virtual void solve() = 0;

protected:
    // object holding the needed field variables for rhs and p
    std::shared_ptr<discretization> discretization_;
    double epsilon_;
    double maxNumberOfIterations_;
};
