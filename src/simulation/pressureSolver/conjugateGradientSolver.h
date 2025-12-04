#pragma once

#include "grid/dataField.h"
#include "grid/staggeredGrid.h"
#include "pressureSolver.h"
#include "simulation/partitioning.h"
#include <memory>

// CG Algorithm according to Chapter 6, Page 11 ParNum Script
class ConjugateGradientSolver final : public PressureSolver {

public:
    ConjugateGradientSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, double epsilon,
                            int maximumNumberOfIterations);

    void solve() override;

private:
    /// rhs := rhs - alpha * Δd
    /// @returns ||rhs||²
    double decreaseResidual(const DataField &d, double alpha);

    /// Applies the Laplace operator (Δd = ∂²d/∂²x + ∂²d/∂²y, i.e. the matrix multiplication of the system: rhs = Δp)
    /// @returns Δd at (i,j)
    double applyDiffusionOperator(const DataField &d, int i, int j) const;

    /// @returns alpha := r^Td/(d^TΔd), r = rhs
    double calculateAlpha();

    //& p := p + alpha * d
    void updatePressure(double alpha);

    /// @returns d := r + beta * d, d = direction_
    void updateDirection(double beta);

    /// Store direction in this object: No need to reallocate each time we solve()
    DataField direction_;

    double dx2_;
    double dy2_;
    double invDx2_;
    double invDy2_;
};
