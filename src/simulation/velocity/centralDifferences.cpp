#include "simulation/velocity/centralDifferences.h"

centralDifferences::centralDifferences(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth)
    : discretization(nCells, meshWidth)
{ }

double centralDifferences::computeDu2Dx(int i, int j) const {
    const double u1 = (u_(i + 1, j) + u_(i, j)) / 2;
    const double u2 = (u_(i - 1, j) + u_(i, j)) / 2;
    return (u1 * u1 - u2 * u2) / dx();
}

double centralDifferences::computeDv2Dy(int i, int j) const {
    const double v1 = (v_(i, j + 1) + v_(i, j)) / 2;
    const double v2 = (v_(i, j - 1) + v_(i, j)) / 2;
    return (v1 * v1 - v2 * v2) / dy();
}

double centralDifferences::computeDuvDx(int i, int j) const {
    const double v1 = (v_(i + 1, j) + v_(i, j)) / 2;
    const double v2 = (v_(i - 1, j) + v_(i, j)) / 2;

    const double u1 = (u_(i, j + 1) + u_(i, j)) / 2;
    const double u2 = (u_(i - 1, j + 1) + u_(i - 1, j)) / 2;
    
    return (v1 * u1 - v2 * u2) / dx();
}

double centralDifferences::computeDuvDy(int i, int j) const {
    const double v1 = (v_(i + 1, j) + v_(i, j)) / 2;
    const double v2 = (v_(i + 1, j - 1) + v_(i, j - 1)) / 2;

    const double u1 = (u_(i, j + 1) + u_(i, j)) / 2;
    const double u2 = (u_(i, j - 1) + u_(i, j)) / 2;
    
    return (v1 * u1 - v2 * u2) / dy();
}
