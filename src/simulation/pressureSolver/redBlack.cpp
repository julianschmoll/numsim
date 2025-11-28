#include "redBlack.h"

#include <limits>

RedBlack::RedBlack(const std::shared_ptr<StaggeredGrid> &grid,
                   const double epsilon,
                   const int maximumNumberOfIterations,
                   const double omega,
                   const std::shared_ptr<Partitioning> &partitioning)
    : partitioning_(partitioning), grid_(grid), epsilon_(epsilon), maxNumberOfIterations_(maximumNumberOfIterations), omega_(omega) {}


void RedBlack::solve() {
    int it = 0;

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1 / dx2;
    const double invDy2 = 1 / dy2;
    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

    const int residualStartThreshold = static_cast<int>(0.7 * lastIterationCount_);
    double squareResidual = std::numeric_limits<double>::max();

    while (it < maxNumberOfIterations_) {
        const int beginJ = p.beginJ() + 1;
        const int endJ = p.endJ() - 1;
        const int beginI = p.beginI() + 1;
        const int endI = p.endI() - 1;

        // rb = 0 is red, rb=1 is black pass
        for (int rb = 0; rb < 2; ++rb) {
            for (int j = beginJ; j < endJ; ++j) {
                // ((beginI + j + rb) & 1) calculates offset 0 or 1
                for (int i = beginI + ((beginI + j + rb) & 1); i < endI; i += 2) {
                    const double pDxx = (p(i - 1, j) + p(i + 1, j)) * invDx2;
                    const double pDyy = (p(i, j - 1) + p(i, j + 1)) * invDy2;
                    const double pNew = (1 - omega_) * p(i, j) + omega_ * scalingFactor * (pDxx + pDyy - rhs(i, j));
                    squareResidual += (p(i, j) - pNew) * (p(i, j) - pNew);
                    p(i, j) = pNew;
                }
            }
            exchange();
        }

        if (it > residualStartThreshold && squareResidual < epsilon_ * epsilon_) {
            break;
        }
        it++;
    }

}

void RedBlack::exchange() {
    // ToDo: This should exchange values with MPI
}

void RedBlack::setBoundaryValues(Direction direction) {
    DataField &p = grid_->p();

    const int beginI = p.beginI();
    const int endI = p.endI();
    const int beginJ = p.beginJ();
    const int endJ = p.endJ();

    switch (direction) {
    case Direction::Left: {
        const int iIndex = beginI;
        const int srcI = beginI + 1;
        for (int j = beginJ; j < endJ; ++j)
            p(iIndex, j) = p(srcI, j);
        break;
    }
    case Direction::Right: {
        const int iIndex = endI;
        const int srcI = endI - 2;
        for (int j = beginJ; j < endJ; ++j)
            p(iIndex, j) = p(srcI, j);
        break;
    }
    case Direction::Bottom: {
        const int jIndex = beginJ;
        const int srcJ = beginJ + 1;
        for (int i = beginI + 1; i < endI - 1; ++i)
            p(i, jIndex) = p(i, srcJ);
        break;
    }
    case Direction::Top: {
        const int jIndex = endJ;
        const int srcJ = endJ - 2;
        for (int i = beginI + 1; i < endI - 1; ++i)
            p(i, jIndex) = p(i, srcJ);
        break;
    }
    }
}
