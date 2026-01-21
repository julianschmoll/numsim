#pragma once

#include "grid/staggeredGrid.h"
#include "pressureSolver.h"
#include "simulation/partitioning.h"
#include "settings.h"
#include <memory>

/**
 * @class RedBlackSolver
 * @brief Iterative red black solver for the Pressure.
 *
 * This class solves the pressure either using Red Black GaussSeidel or Red Black  SOR,
 * depending on the relaxation factor omega. If omega equals 1,
 * the pressure is solved according to GaussSeidel, if omega does
 * not equal one, the pressure is solved according to SOR.
 */
class RedBlackSolver final : public PressureSolver {
public:
    /**
     * Constructs red black solver.
     *
     * @param grid Grid holding the needed data fields.
     * @param partitioning Information on how the grid is partitioned.
     * @param epsilon Convergence threshold for the residual.
     * @param maximumNumberOfIterations Maximum number of iterations.
     * @param omega Relaxation factor (1.0 for Gauss-Seidel, otherwise SOR).
     */
    RedBlackSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings);

    /**
     * Solves the poisson problem for the pressure.
     */
    void solve() override;

private:
    /// Relaxation factor
    double omega_;

    /// Error on the own partition.
    double localPressureError_ = 0;

    /// Global error sum.
    double globalPressureError_ = 0;
};
