#include "simulation/pressure/PressureSolver.h"

PressureSolver::PressureSolver(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations, double omega)
    : discretization_(discretization), epsilon_(epsilon), maxNumberOfIterations_(maxNumberOfIterations), omega_(omega) 
{ }

double PressureSolver::calculateSquareResidual() const {
    double squareResidual = 0;

    dataField &p = discretization_->p();
    dataField &rhs = discretization_->rhs();

    double dx2 = discretization_->dx() * discretization_->dx();
    double dy2 = discretization_->dy() * discretization_->dy();

    for (int j = p.beginJ(); j < p.endJ(); j++) {
        for (int i = p.beginI(); i < p.endI(); i++) {
            const double p_xx = (p(i - 1, j) - 2 * p(i, j) + p(i + 1, j)) / dx2;
            const double p_yy = (p(i, j - 1) - 2 * p(i, j) + p(i, j + 1)) / dy2;
            const double difference = rhs(i, j) - p_xx - p_yy;
            squareResidual += difference * difference;
        }
    }
    return squareResidual;
}

void PressureSolver::solve() {
    int it = 0;
    double squareResidual = 0;

    dataField &p = discretization_->p();
    dataField &rhs = discretization_->rhs();

    double dx2 = discretization_->dx() * discretization_->dx();
    double dy2 = discretization_->dy() * discretization_->dy();

    double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

#ifdef RESIDUAL_METHOD
    squareResidual = calculateSquareResidual();
#else
    squareResidual = std::numeric_limits<double>::max();
#endif

    while (it < maxNumberOfIterations_ && squareResidual > epsilon_ * epsilon_) {
        squareResidual = 0;
        for (int j = p.beginJ(); j < p.endJ(); j++) {
            for (int i = p.beginI(); i < p.endI(); i++) {
                const double pDxx = (p(i - 1, j) + p(i + 1, j)) / dx2;
                const double pDyy = (p(i, j - 1) + p(i, j + 1)) / dy2;
                const double pNew = (1 - omega_) * p(i, j) + omega_ * scalingFactor * (pDxx + pDyy - rhs(i, j));
#ifndef RESIDUAL_METHOD
                squareResidual += (p(i, j) - pNew) * (p(i, j) - pNew);
#endif
                p(i, j) = pNew;
            }
        }
#ifdef RESIDUAL_METHOD
        squareResidual = calculateSquareResidual();
#endif
        // TODO: Here correct?
        setBoundaryValues();
        it++;
    }
}

void PressureSolver::setBoundaryValues() {
    dataField &p = discretization_->p();
    // bottom
    for (int i = p.beginI(); i < p.endI(); i++) {
        p(i, 0) = p(i, 1);
    }
    // left / right
    int n = p.endJ();
    for (int j = p.beginJ(); j < p.endJ(); j++) {
        p(0, j) = p(1, j);
        p(n, j) = p(n - 1, j);
    }
    // top
    n = p.endJ();
    for (int i = p.beginI(); i < p.endI(); i++) {
        p(i, n) = p(i, n - 1);
    }
}



