#include "simulation/pressureSolver.h"

#include <utility>

PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon,
                               const double maxNumberOfIterations, const double omega)
    : discretization_(std::move(discretization)), epsilon_(epsilon), maxNumberOfIterations_(maxNumberOfIterations),
      omega_(omega) {}

double PressureSolver::calculateSquareResidual() const {
    double squareResidual = 0;

    DataField& p = discretization_->p();
    DataField& rhs = discretization_->rhs();

    double dx2 = discretization_->dx() * discretization_->dx();
    double dy2 = discretization_->dy() * discretization_->dy();

    for (int j = p.beginJ() + 1; j < p.endJ() - 1; j++) {
        for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
            const double p_xx = (p(i - 1, j) - 2 * p(i, j) + p(i + 1, j)) / dx2;
            const double p_yy = (p(i, j - 1) - 2 * p(i, j) + p(i, j + 1)) / dy2;
            const double difference = rhs(i, j) - p_xx - p_yy;
            squareResidual += difference * difference;
        }
    }

    const int nCellsX = discretization_->nCells()[0] - 2;
    const int nCellsY = discretization_->nCells()[1] - 2;
    const auto nCellsTotal = static_cast<double>(nCellsX * nCellsY);

    return squareResidual / nCellsTotal;
}

void PressureSolver::solve() {
    int it = 0;
    double squareResidual = 0;

    DataField& p = discretization_->p();
    DataField& rhs = discretization_->rhs();

    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();

    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

#ifdef RESIDUAL_METHOD
    squareResidual = calculateSquareResidual();
#else
    squareResidual = std::numeric_limits<double>::max();
#endif

    while (it < maxNumberOfIterations_ && squareResidual > epsilon_ * epsilon_) {
        for (int j = p.beginJ() + 1; j < p.endJ() - 1; j++) {
            for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
                const double pDxx = (p(i - 1, j) + p(i + 1, j)) / dx2;
                const double pDyy = (p(i, j - 1) + p(i, j + 1)) / dy2;
                const double pNew = (1 - omega_) * p(i, j) + omega_ * scalingFactor * (pDxx + pDyy - rhs(i, j));
#ifndef RESIDUAL_METHOD
                squareResidual += (p(i, j) - pNew) * (p(i, j) - pNew);
#endif
                p(i, j) = pNew;
            }
        }

        setBoundaryValues();
#ifdef RESIDUAL_METHOD
        squareResidual = calculateSquareResidual();
#endif
        it++;
    }
}

void PressureSolver::setBoundaryValues() {
    DataField& p = discretization_->p();
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



