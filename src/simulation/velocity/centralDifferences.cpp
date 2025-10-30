#include "simulation/velocity/centralDifferences.h"

centralDifferences::centralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth) :
        discretization(nCells, meshWidth) {}

double centralDifferences::computeDu2Dx(int i, int j) {
    double uCenter = u(i, j);
    double uLeft = u(i - 1, j);
    double uRight = u(i + 1, j);
    double uHalfRight = (uRight + uCenter) / 2;
    double uHalfLeft = (uCenter + uLeft) / 2;

    return (uHalfRight * uHalfRight - uHalfLeft * uHalfLeft) / dx();
}

double centralDifferences::computeDv2Dy(int i, int j) {
    double vCenter = v(i, j);
    double vUp = v(i, j + 1);
    double vLow = v(i, j - 1);
    double vHalfUp = (vUp + vCenter) / 2;
    double vHalfLow = (vCenter + vLow) / 2;

    return (vHalfUp * vHalfUp - vHalfLow * vHalfLow) / dy();
}

double centralDifferences::computeDuvDx(int i, int j) {
    double vCenter = v(i, j);
    double vLeft = v(i - 1, j);
    double vRight = v(i + 1, j);

    double uLeftLow = u(i - 1, j);
    double uLeftUp = u(i - 1, j + 1);
    double uRightLow = u(i, j);
    double uRightUp = u(i, j + 1);

    double vHalfLeft = (vLeft + vCenter) / 2;
    double vHalfRight = (vCenter + vRight) / 2;
    double uHalfLeft = (uLeftLow + uLeftUp) / 2;
    double uHalfRight = (uRightLow + uRightUp) / 2;

    return (1 / dx()) * (vHalfRight * uHalfRight - vHalfLeft * uHalfLeft);
}

double centralDifferences::computeDuvDy(int i, int j) {
    double uCenter = u(i, j);
    double uLow = u(i, j - 1);
    double uUp = u(i, j + 1);

    double vUpLeft = v(i, j);
    double vUpRight = v(i + 1, j);
    double vLowLeft = v(i, j - 1);
    double vLowRight = v(i + 1, j - 1);

    double vHalfUp = (vUpLeft + vUpRight) / 2;
    double vHalfLow = (vLowLeft + vLowRight) / 2;
    double uHalfUp = (uUp + uCenter) / 2;
    double uHalfLow = (uLow + uCenter) / 2;

    return (1 / dy()) * (vHalfUp * uHalfUp - vHalfLow * uHalfLow);
}
