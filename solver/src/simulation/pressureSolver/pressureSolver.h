#pragma once

#include "grid/staggeredGrid.h"
#include "simulation/partitioning.h"
#include <memory>

/**
 * @class PressureSolver
 * @brief Base class for implementation of pressure solvers.
 */
class PressureSolver {

public:
    /**
     * Constructs Pressure solver.
     *
     * @param grid Grid holding the needed data fields.
     * @param partitioning Information on how the grid is partitioned.
     * @param epsilon Convergence threshold for the residual.
     * @param maximumNumberOfIterations Maximum number of iterations.
     */
    PressureSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, double epsilon, int maximumNumberOfIterations);

    /**
     * Destructs Pressure solver.
     */
    virtual ~PressureSolver() = default;

    // We don't want pressure solver to be moved or copied
    PressureSolver(const PressureSolver &) = delete;
    PressureSolver(PressureSolver &&) = delete;

    /**
     * Virtual method to override in pressure solver implementations.
     */
    virtual void solve() = 0;

protected:
    /**
     * Sets Neumann boundary values on the pressure field.
     */
    void setBoundaryValues();

    /// Pointer to partitioning
    std::shared_ptr<Partitioning> partitioning_;

    /// object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;

    /// Convergence threshold for the residual
    double epsilon_;

    /// Maximum number of iterations allowed in the solver
    double maxNumberOfIterations_;
};
