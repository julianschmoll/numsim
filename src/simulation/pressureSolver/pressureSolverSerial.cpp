#include "simulation/pressureSolver/pressureSolverSerial.h"
#include "macros.h"
#include <iostream>

PressureSolverSerial::PressureSolverSerial(const std::shared_ptr<StaggeredGrid> &grid, const double epsilon, const double maxNumberOfIterations, const double omega)
    : grid_(grid), epsilon_(epsilon), maxNumberOfIterations_(maxNumberOfIterations), omega_(omega) {}

double PressureSolverSerial::calculateSquareResidual() const {
    double squareResidual = 0.0;

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1.0 / dx2;
    const double invDy2 = 1.0 / dy2;

    for (int j = p.beginJ() + 1; j < p.endJ() - 1; ++j) {
        for (int i = p.beginI() + 1; i < p.endI() - 1; ++i) {
            const double pij = p(i, j);
            const double pDxx = (p(i - 1, j) - 2.0 * pij + p(i + 1, j)) * invDx2;
            const double pDyy = (p(i, j - 1) - 2.0 * pij + p(i, j + 1)) * invDy2;
            const double diff = rhs(i, j) - pDxx - pDyy;
            squareResidual += diff * diff;
        }
    }

    const int nCellsX = grid_->nCells()[0] - 2;
    const int nCellsY = grid_->nCells()[1] - 2;
    const auto nCellsTotal = static_cast<double>(nCellsX * nCellsY);

    return squareResidual / nCellsTotal;
}

void PressureSolverSerial::solve() {
    int it = 0;

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1 / dx2;
    const double invDy2 = 1 / dy2;
    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

    const int residualStartThreshold = static_cast<int>(0.7 * lastIterationCount_);

    while (it < maxNumberOfIterations_) {
        for (int j = p.beginJ() + 1; j < p.endJ() - 1; j++) {
            for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
                const double pDxx = (p(i - 1, j) + p(i + 1, j)) * invDx2;
                const double pDyy = (p(i, j - 1) + p(i, j + 1)) * invDy2;
                p(i, j) = (1 - omega_) * p(i, j) + omega_ * scalingFactor * (pDxx + pDyy - rhs(i, j));
            }
        }

        setBoundaryValues();

        if (it > residualStartThreshold && calculateSquareResidual() < epsilon_ * epsilon_) {
            break;
        }

        it++;
    }

    DEBUG(std::cout << "PressureSolverSerial::solve():  it=" << it << ",  res²=" << calculateSquareResidual() << std::endl);

    lastIterationCount_ = it;
}

// TODO: Randwerte eventuell direkt in solver() Schleife setzen.
void PressureSolverSerial::setBoundaryValues() {
    DataField &p = grid_->p();
    for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
        // bottom
        p(i, p.beginJ()) = p(i, p.beginJ() + 1);
        // top
        p(i, p.endJ() - 1) = p(i, p.endJ() - 2);
    }

    for (int j = p.beginJ(); j < p.endJ(); j++) {
        // left
        p(p.beginI(), j) = p(p.beginI() + 1, j);
        // right
        p(p.endI() - 1, j) = p(p.endI() - 2, j);
    }
}