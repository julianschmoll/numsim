#include "simulation/velocity/centralDifferences.h"

CentralDifferences::CentralDifferences(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth)
    : Discretization(nCells, meshWidth)
{ }

double CentralDifferences::computeDu2Dx(int i, int j) const {
    const double uHalfRight = (u_(i + 1, j) + u_(i, j)) / 2;
    const double uHalfLeft  = (u_(i - 1, j) + u_(i, j)) / 2;

    return (uHalfRight * uHalfRight - uHalfLeft * uHalfLeft) / dx();
}

double CentralDifferences::computeDv2Dy(int i, int j) const {
    const double vHalfUp   = (v_(i, j + 1) + v_(i, j)) / 2;
    const double vHalfDown = (v_(i, j - 1) + v_(i, j)) / 2;

    return (vHalfUp * vHalfUp - vHalfDown * vHalfDown) / dy();
}

double CentralDifferences::computeDuvDx(int i, int j) const {
    const double vHalfRight = (v_(i + 1, j) + v_(i, j)) / 2;
    const double vHalfLeft  = (v_(i - 1, j) + v_(i, j)) / 2;

    const double uHalfUp     = (u_(i, j + 1) + u_(i, j)) / 2;
    const double uHalfUpLeft = (u_(i - 1, j + 1) + u_(i - 1, j)) / 2;
    
    return (vHalfRight * uHalfUp - vHalfLeft * uHalfUpLeft) / dx();
}

double CentralDifferences::computeDuvDy(int i, int j) const {
    const double vHalfRight     = (v_(i + 1, j) + v_(i, j)) / 2;
    const double vHalfRightDown = (v_(i, j - 1) + v_(i + 1, j - 1)) / 2;

    const double uHalfUp   = (u_(i, j + 1) + u_(i, j)) / 2;
    const double uHalfDown = (u_(i, j) + u_(i, j - 1)) / 2;

    return (vHalfRight * uHalfUp - vHalfRightDown * uHalfDown) / dy();
}
