#include "donorCell.h"

donorCell::donorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha) :
        discretization(nCells, meshWidth),
        alpha_(alpha) {}

double donorCell::computeDu2Dx(int i, int j) {
    double uCenter = u(i, j);
    double uLeft = u(i - 1, j);
    double uRight = u(i + 1, j);
    double uHalfRight = (uCenter + uRight) / 2;
    double uHalfLeft = (uLeft + uCenter) / 2;

    double centDiffPart = (uHalfRight * uHalfRight - uHalfLeft * uHalfLeft) / dx();

    return centDiffPart +
           (alpha_ / dx()) *
           (abs(uHalfRight) * ((uCenter - uRight) / 2) - abs(uHalfLeft) * ((uLeft - uCenter) / 2));
}

double donorCell::computeDv2Dy(int i, int j) {
    double vCenter = v(i, j);
    double vUp = v(i, j + 1);
    double vLow = v(i, j - 1);
    double vHalfUp = (vCenter + vUp) / 2;
    double vHalfLow = (vLow + vCenter) / 2;

    double centDiffPart = (vHalfUp * vHalfUp - vHalfLow * vHalfLow) / dy();

    return centDiffPart +
           (alpha_ / dy()) *
           (abs(vHalfUp) * ((vCenter - vUp) / 2) - abs(vHalfLow) * ((vLow - vCenter) / 2));
}

double donorCell::computeDuvDx(int i, int j) {
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

    double centDiffPart = (1 / dx()) * (vHalfRight * uHalfRight - vHalfLeft * uHalfLeft);

    return centDiffPart +
           (alpha_ / dx()) *
           (abs(vHalfRight) * ((uRightLow - uRightUp) / 2) - abs(vHalfLeft) * ((uLeftLow - uLeftUp) / 2));
};

double donorCell::computeDuvDy(int i, int j) {
    double uCenter = u(i, j);
    double uLow = u(i, j - 1);
    double uUp = u(i, j + 1);

    double vUpLeft = v(i, j);
    double vUpRight = v(i + 1, j);
    double vLowLeft = v(i, j - 1);
    double vLowRight = v(i + 1, j - 1);

    double vHalfUp = (vUpLeft + vUpRight) / 2;
    double vHalfLow = (vLowLeft + vLowRight) / 2;
    double uHalfUp = (uCenter + uUp) / 2;
    double uHalfLow = (uLow + uCenter) / 2;

    double centDiffPart = (1 / dy()) * (vHalfUp * uHalfUp - vHalfLow * uHalfLow);

    return centDiffPart +
           (alpha_ / dy()) *
           (abs(vHalfUp) * ((uCenter - uUp) / 2) - abs(vHalfLow) * ((uLow - uCenter) / 2));
}
