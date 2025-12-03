#include "simulation/discreteOperators.h"
#include <cstdlib>

DiscreteOperators::DiscreteOperators(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, const Partitioning &partitioning,
                                     double alpha)
    : StaggeredGrid(nCells, meshWidth, partitioning), alpha_(alpha) {}

double DiscreteOperators::computeDu2Dx(const int i, const int j) const {
    const double uHalfRight = (u_(i + 1, j) + u_(i, j)) / 2;
    const double uHalfLeft = (u_(i - 1, j) + u_(i, j)) / 2;

    const double centralDifferenceDerivative = (uHalfRight * uHalfRight - uHalfLeft * uHalfLeft) / dx();

    const double uDiffRight = (u_(i, j) - u_(i + 1, j)) / 2;
    const double uDiffLeft = (u_(i - 1, j) - u_(i, j)) / 2;

    const double donorCellContribution = (std::abs(uHalfRight) * uDiffRight - std::abs(uHalfLeft) * uDiffLeft) / dx();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double DiscreteOperators::computeDv2Dy(const int i, const int j) const {
    const double vHalfUp = (v_(i, j + 1) + v_(i, j)) / 2;
    const double vHalfDown = (v_(i, j - 1) + v_(i, j)) / 2;

    const double centralDifferenceDerivative = (vHalfUp * vHalfUp - vHalfDown * vHalfDown) / dy();

    const double vDiffUp = (v_(i, j) - v_(i, j + 1)) / 2;
    const double vDiffDown = (v_(i, j - 1) - v_(i, j)) / 2;

    const double donorCellContribution = (std::abs(vHalfUp) * vDiffUp - std::abs(vHalfDown) * vDiffDown) / dy();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double DiscreteOperators::computeDuvDx(const int i, const int j) const {
    const double uHalfUp = (u_(i, j + 1) + u_(i, j)) / 2;
    const double uHalfUpLeft = (u_(i - 1, j + 1) + u_(i - 1, j)) / 2;

    const double vHalfRight = (v_(i + 1, j) + v_(i, j)) / 2;
    const double vHalfLeft = (v_(i - 1, j) + v_(i, j)) / 2;

    const double vDiffRight = (v_(i, j) - v_(i + 1, j)) / 2;
    const double vDiffLeft = (v_(i - 1, j) - v_(i, j)) / 2;

    const double centralDifferenceDerivative = (vHalfRight * uHalfUp - vHalfLeft * uHalfUpLeft) / dx();

    const double donorCellContribution = (std::abs(uHalfUp) * vDiffRight - std::abs(uHalfUpLeft) * vDiffLeft) / dx();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double DiscreteOperators::computeDuvDy(const int i, const int j) const {
    const double vHalfRight = (v_(i + 1, j) + v_(i, j)) / 2;
    const double vHalfRightDown = (v_(i, j - 1) + v_(i + 1, j - 1)) / 2;

    const double uHalfUp = (u_(i, j + 1) + u_(i, j)) / 2;
    const double uHalfDown = (u_(i, j) + u_(i, j - 1)) / 2;

    const double uDiffUp = (u_(i, j) - u_(i, j + 1)) / 2;
    const double uDiffDown = (u_(i, j - 1) - u_(i, j)) / 2;

    const double centralDifferenceDerivative = (vHalfRight * uHalfUp - vHalfRightDown * uHalfDown) / dy();

    const double donorCellContribution = (std::abs(vHalfRight) * uDiffUp - std::abs(vHalfRightDown) * uDiffDown) / dy();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double DiscreteOperators::computeD2uDx2(const int i, const int j) const {
    const double du2 = u_(i + 1, j) - 2 * u_(i, j) + u_(i - 1, j);
    return du2 / (dx() * dx());
}

double DiscreteOperators::computeD2uDy2(const int i, const int j) const {
    const double du2 = u_(i, j + 1) - 2 * u_(i, j) + u_(i, j - 1);
    return du2 / (dy() * dy());
}

double DiscreteOperators::computeD2vDx2(const int i, const int j) const {
    const double dv2 = v_(i + 1, j) - 2 * v_(i, j) + v_(i - 1, j);
    return dv2 / (dx() * dx());
}

double DiscreteOperators::computeD2vDy2(const int i, const int j) const {
    const double dv2 = v_(i, j + 1) - 2 * v_(i, j) + v_(i, j - 1);
    return dv2 / (dy() * dy());
}

double DiscreteOperators::computeDuDx(const int i, const int j) const {
    const double du = u_(i, j) - u_(i - 1, j);
    return du / dx();
}

double DiscreteOperators::computeDvDy(const int i, const int j) const {
    const double dv = v_(i, j) - v_(i, j - 1);
    return dv / dy();
}

double DiscreteOperators::computeDpDx(const int i, const int j) const {
    const double dp = p_(i + 1, j) - p_(i, j);
    return dp / dx();
}

double DiscreteOperators::computeDpDy(const int i, const int j) const {
    const double dp = p_(i, j + 1) - p_(i, j);
    return dp / dy();
}
