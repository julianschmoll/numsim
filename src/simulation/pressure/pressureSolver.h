#pragma once

#include "grid/discretization.h"
#include <memory>

#define RESIDUAL_METHOD

class PressureSolver {
public:
    virtual ~PressureSolver() = default;
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, double maxNumberOfIterations, double omega);

    /**
     * Solves the Poisson problem for the pressure.
     *
     * Uses the rhs and p field variables in the grid.
     */
    void solve();

    /**
     * Sets the boundary values to account for homogenous Neumann boundary conditions.
     *
     * This has to be called after every iteration.
     */
    void setBoundaryValues();

    double calculateSquareResidual() const;

protected:
    // object holding the needed field variables for rhs and p
    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    double maxNumberOfIterations_;
    double omega_;
};
