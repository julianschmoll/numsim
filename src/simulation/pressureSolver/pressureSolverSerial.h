#pragma once

#include "grid/staggeredGrid.h"
#include <memory>

#define RESIDUAL_METHOD

/**
 * @class PressureSolverSerial
 * @brief Iterative solver for the Pressure.
 *
 * This class solves the pressure either using GaussSeidel or SOR,
 * depending on the relaxation factor omega. If omega equals 1,
 * the pressure is solved according to GaussSeidel, if omega does
 * not equal one, the pressure is solved according to SOR.
 */
class PressureSolverSerial {
public:
    /**
     * Constructs a PressureSolverSerial with the given grid and solver parameters.
     *
     * @param grid Shared pointer to the grid grid.
     * @param epsilon Convergence threshold for the residual.
     * @param maxNumberOfIterations Maximum number of iterations.
     * @param omega Relaxation factor (1.0 for Gauss-Seidel, otherwise SOR).
     */
    PressureSolverSerial(const std::shared_ptr<StaggeredGrid> &grid,
        double epsilon,
        double maxNumberOfIterations,
        double omega);

    /**
     * Destructs PressureSolverSerial object.
     */
    ~PressureSolverSerial() = default;

    /**
     * Solves the Poisson problem for the pressure.
     *
     * Uses the rhs and p field variables in the grid.
     */
    void solve();

private:
    /**
     * Sets the boundary values to account for homogenous Neumann boundary conditions.
     *
     * This has to be called after every iteration.
     */
    void setBoundaryValues();

    /**
     * Calculates the squared residual of the pressure field.
     *
     * @return Squared residual value.
     */
    [[nodiscard]] double calculateSquareResidual() const;

    int lastIterationCount_ = 0;

protected:
    // object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;

    // Convergence threshold for the residual
    double epsilon_;

    // Maximum number of iterations allowed in the solver
    double maxNumberOfIterations_;

    // Relaxation factor
    double omega_;
};
