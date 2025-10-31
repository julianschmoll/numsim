#include "simulation/velocity/donorCell.h"

donorCell::donorCell(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, double alpha)
    : discretization(nCells, meshWidth), alpha_(alpha)
{ }

double donorCell::computeDu2Dx(int i, int j) const {
    const double uHalfRight = (u_(i + 1, j) + u_(i, j)) / 2;
    const double uHalfLeft  = (u_(i - 1, j) + u_(i, j)) / 2;

    const double centralDifferenceDerivative = (uHalfRight * uHalfRight - uHalfLeft * uHalfLeft) / dx();

    const double uDiffRight = (u_(i, j) - u_(i + 1, j)) / 2;
    const double uDiffLeft  = (u_(i - 1, j) - u_(i, j)) / 2;

    const double donorCellContribution = (std::abs(uHalfRight) * uDiffRight - std::abs(uHalfLeft) * uDiffLeft) / dx();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double donorCell::computeDv2Dy(int i, int j) const {
    const double vHalfUp   = (v_(i, j + 1) + v_(i, j)) / 2;
    const double vHalfDown = (v_(i, j - 1) + v_(i, j)) / 2;

    const double centralDifferenceDerivative = (vHalfUp * vHalfUp - vHalfDown * vHalfDown) / dx();

    const double vDiffUp   = (v_(i, j) - v_(i, j + 1)) / 2;
    const double vDiffDown = (v_(i, j - 1) - v_(i, j)) / 2;

    const double donorCellContribution = (std::abs(vHalfUp) * vDiffUp - std::abs(vHalfDown) * vDiffDown) / dx();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double donorCell::computeDuvDx(int i, int j) const {
    const double uHalfUp     = (u_(i, j + 1) + u_(i, j)) / 2;
    const double uHalfUpLeft = (u_(i - 1, j + 1) + u_(i - 1, j)) / 2;

    const double vHalfRight = (v_(i + 1, j) + v_(i, j)) / 2;
    const double vHalfLeft  = (v_(i - 1, j) + v_(i, j)) / 2;

    const double vDiffRight = (v_(i, j) - v_(i + 1, j)) / 2;
    const double vDiffLeft  = (v_(i - 1, j) - v_(i, j)) / 2;

    const double centralDifferenceDerivative = (vHalfRight * uHalfUp - vHalfLeft * uHalfUpLeft) / dx();

    const double donorCellContribution = (std::abs(uHalfUp) * vDiffRight - std::abs(uHalfUpLeft) * vDiffLeft)/ dx();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}

double donorCell::computeDuvDy(int i, int j) const {
    const double vHalfRight     = (v_(i + 1, j) + v_(i, j)) / 2;
    const double vHalfRightDown = (v_(i, j - 1) + v_(i + 1, j - 1)) / 2;

    const double uHalfUp   = (u_(i, j + 1) + u_(i, j)) / 2;
    const double uHalfDown = (u_(i, j) + u_(i, j - 1)) / 2;

    const double uDiffUp   = (u_(i, j) - u_(i, j + 1)) / 2;
    const double uDiffDown = (u_(i, j - 1) - u_(i, j)) / 2;

    const double centralDifferenceDerivative = (vHalfRight * uHalfUp - vHalfRightDown * uHalfDown) / dy();

    const double donorCellContribution = (fabs(vHalfRight) * uDiffUp - fabs(vHalfRightDown) * uDiffDown)/ dy();

    return centralDifferenceDerivative + alpha_ * donorCellContribution;
}