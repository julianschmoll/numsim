#include "grid/dataField.h"
#include "grid/discretization.h"

#include <cassert>
#include <cmath>

#include <iostream>

dataField::dataField(const std::array<int, 2> size, const std::array<double, 2> meshWidth, const std::array<double, 2> offset, std::array<int, 2> iIndexRange, std::array<int, 2> jIndexRange)
    : array2d(size), meshWidth_(meshWidth), offset_(offset), iIndexRange_(iIndexRange), jIndexRange_(jIndexRange)
{}

int dataField::beginJ() const {
    return jIndexRange_[0];
}

int dataField::endJ() const {
    return jIndexRange_[1];
}

int dataField::beginI() const {
    return iIndexRange_[0];
}

int dataField::endI() const {
    return iIndexRange_[1];
}

double dataField::interpolateAt(double x, double y) const {

    assert(0 <= x && x < meshWidth_[0] * (size_[0] - 1));  // u, v, p brauchen gleiche größe
    assert(0 <= y && y < meshWidth_[1] * (size_[1] - 1));

    x /= meshWidth_[0];
    y /= meshWidth_[1];

    // Als Test, Druck:

    x += offset_[0]; // p: 0.5, u: 0
    y += offset_[1]; // p: 0.5, u: 0.5

    assert(0 <= x && x < size_[0]);
    assert(0 <= y && y < size_[1]);

    const int iLower = static_cast<int>(std::floor(x));
    const int iUpper = static_cast<int>(std::ceil(x));
    const double alpha = x - iLower;

    const int jLower = static_cast<int>(std::floor(y));
    const int jUpper = static_cast<int>(std::ceil(y));
    const double beta = y - jLower;

    assert(0 <= iLower && iLower < size_[0]);
    assert(0 <= iUpper && iUpper < size_[0]);
    assert(0 <= jLower && jLower < size_[1]);
    assert(0 <= jUpper && jUpper < size_[1]);

    const double d11 = data_.at(jLower * size_[0] + iLower);
    const double d21 = data_.at(jLower * size_[0] + iUpper);
    const double d12 = data_.at(jUpper * size_[0] + iLower);
    const double d22 = data_.at(jUpper * size_[0] + iUpper);

    const double dx1 = (1 - alpha) * d11 + alpha * d21;
    const double dx2 = (1 - alpha) * d12 + alpha * d22;

    return (1 - beta) * dx1 + beta * dx2;
}
