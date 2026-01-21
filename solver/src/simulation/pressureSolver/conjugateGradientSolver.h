#pragma once

#include "grid/dataField.h"
#include "grid/staggeredGrid.h"
#include "pressureSolver.h"
#include "simulation/partitioning.h"
#include <memory>

// CG Algorithm according to Chapter 6, Page 11 ParNum Script
/**
 * @class ConjugateGradientSolver
 * @brief CG solver for pressure.
 */
class ConjugateGradientSolver final : public PressureSolver {

public:
    /**
     * Constructs CG Solver.
     *
     * @param grid Grid holding the needed data fields.
     * @param partitioning Information on how the grid is partitioned.
     * @param epsilon Convergence threshold for the residual.
     * @param maximumNumberOfIterations Maximum number of iterations.
     */
    ConjugateGradientSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings);

    /**
     * Solves pressure based on Conjugate Gradient Method.
     */
    void solve() override;

private:
    /**
     * Updates the residual field according to r := r - d.
     *
     * @param d The current search direction d.
     * @param alpha The calculated step size a.
     * @returns The squared L_2-norm of the new residual.
     */
    double decreaseResidual(const DataField &d, double alpha);

    /**
     * Applies the Laplace operator to the given DataField at a specific grid point.
     *
     * @param d The direction vector to apply the operator to.
     * @param i The x-index of the grid point.
     * @param j The y-index of the grid point.
     * @returns The value of the Laplace Operator at grid point.
     */
    double applyDiffusionOperator(const DataField &d, int i, int j) const;

    /**
     * Calculates the step size for the current iteration.
     *
     * @returns The calculated step size.
     */
    double calculateAlpha();

    /**
     * Updates the pressure field according to p := p + a d.
     *
     * @param alpha The calculated step size a.
     */
    void updatePressure(double alpha);

    /**
     * Updates the search direction field according to d:= r + bd, where r is the new residual.
     *
     * @param beta The calculated scaling factor.
     */
    void updateDirection(double beta);

    /// Store direction in this object: No need to reallocate each time we solve()
    DataField direction_;

    double dx2_;
    double dy2_;
    double invDx2_;
    double invDy2_;
};
