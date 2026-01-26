#pragma once

#include "grid/staggeredGrid.h"
#include "simulation/partitioning.h"
#include "settings.h"
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
    PressureSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings);

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
    virtual void solve(DataField &p) = 0;

protected:
    /**
     * Sets Neumann boundary values on the pressure field.
     */
    void setBoundaryValues(DataField &p);

    /**
     * no slip at top and bottom structure boundary
     * call after other boundary methods
     */
    void setStructureBoundaries(DataField &p);

    /// Pointer to partitioning
    std::shared_ptr<Partitioning> partitioning_;

    // Configuration settings for the simulation
    Settings settings_;

    /// object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;
};
