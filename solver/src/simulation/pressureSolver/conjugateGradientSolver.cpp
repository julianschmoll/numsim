#include "conjugateGradientSolver.h"
#include "grid/dataField.h"
#include <cmath>

ConjugateGradientSolver::ConjugateGradientSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings)
    : PressureSolver(std::move(grid), std::move(partitioning), settings), direction_(grid_->p().size(), grid_->meshWidth()),
      dx2_(grid_->dx() * grid_->dx()), dy2_(grid_->dy() * grid_->dy()), invDx2_(1.0 / dx2_), invDy2_(1.0 / dy2_) {}

double ConjugateGradientSolver::applyDiffusionOperator(const DataField &d, int i, int j) const {
    const double pij = d(i, j);
    const double pDxx = (d(i - 1, j) - 2.0 * pij + d(i + 1, j)) * invDx2_;
    const double pDyy = (d(i, j - 1) - 2.0 * pij + d(i, j + 1)) * invDy2_;
    return pDxx + pDyy;
}

void ConjugateGradientSolver::updatePressure(double alpha) {
    DataField &p = grid_->p();
    DataField &d = direction_;

    // paralleisierbar
    for (int j = d.beginJ() + 1; j < d.endJ() - 1; ++j) {
        for (int i = d.beginI() + 1; i < d.endI() - 1; ++i) {
            if (grid_->isSolid(i, j)) continue;
            p(i, j) += alpha * d(i, j);
        }
    }
    setBoundaryValues();
    partitioning_->exchange(p);
}

void ConjugateGradientSolver::updateDirection(double beta) {
    DataField &rhs = grid_->rhs();
    DataField &d = direction_;

    // paralleisierbar
    for (int j = d.beginJ() + 1; j < d.endJ() - 1; ++j) {
        for (int i = d.beginI() + 1; i < d.endI() - 1; ++i) {
            if (grid_->isSolid(i, j)) continue;
            d(i, j) = rhs(i, j) + beta * d(i, j);
        }
    }
    partitioning_->exchange(d);
}

double ConjugateGradientSolver::decreaseResidual(const DataField &d, double alpha) {
    double localSquareResidual = 0.0;

    DataField &rhs = grid_->rhs();

    // paralleisierbar
    for (int j = d.beginJ() + 1; j < d.endJ() - 1; ++j) {
        for (int i = d.beginI() + 1; i < d.endI() - 1; ++i) {
            if (grid_->isSolid(i, j)) continue;
            rhs(i, j) -= alpha * applyDiffusionOperator(d, i, j); // rhs becomes residual
            localSquareResidual += rhs(i, j) * rhs(i, j);
        }
    }
    partitioning_->exchange(rhs);
    double squareResidual = partitioning_->collectSum(localSquareResidual);

    return squareResidual;
}

double ConjugateGradientSolver::calculateAlpha() {

    DataField &rhs = grid_->rhs();
    DataField &d = direction_;

    double rdDotLocal = 0;
    double dAdDotLocal = 0;

    // paralleisierbar
    for (int j = d.beginJ() + 1; j < d.endJ() - 1; ++j) {
        for (int i = d.beginI() + 1; i < d.endI() - 1; ++i) {
            if (grid_->isSolid(i, j)) continue;
            rdDotLocal += rhs(i, j) * d(i, j);
            dAdDotLocal += d(i, j) * applyDiffusionOperator(d, i, j);
        }
    }
    double rdDot = partitioning_->collectSum(rdDotLocal);
    double dAdDot = partitioning_->collectSum(dAdDotLocal);
    return rdDot / dAdDot;
}

// alpha = gradient descent strength
// beta = direction update parameter
void ConjugateGradientSolver::solve() {

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();
    DataField &d = direction_;

    // initial residual: r = rhs - Δp
    const double r0 = decreaseResidual(p, 1);
    double rNew = r0;

    // initial descent direction equals residual
    for (int j = d.beginJ() + 1; j < d.endJ() - 1; ++j) {
        for (int i = d.beginI() + 1; i < d.endI() - 1; ++i) {
            if (grid_->isSolid(i, j)) continue;
            d(i, j) = rhs(i, j);
        }
    }
    partitioning_->exchange(d);

    while (rNew > settings_.epsilon * settings_.epsilon) {
        const double alpha = calculateAlpha();
        updatePressure(alpha);
        const double rOld = rNew;
        rNew = decreaseResidual(d, alpha);
        const double beta = rNew / rOld;
        updateDirection(beta);
    }
}
