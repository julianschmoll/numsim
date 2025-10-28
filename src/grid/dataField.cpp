#include "grid/dataField.h"

#include <cassert>
#include <cmath>

dataField::dataField(const std::array<int, 2> size, const std::array<double, 2> meshWidth, const std::array<double, 2> offset)
    : array2d(size), meshWidth_(meshWidth), offset_(offset) {}

double dataField::interpolateAt(double x, double y) const {

    assert(0 <= x && x < meshWidth_[0]);
    assert(0 <= y && y < meshWidth_[1]);

    const double cellWidth = meshWidth_[0] / (size_[0] - 2);  // u, v, p brauchen gleiche größe
    const double cellHeight = meshWidth_[1] / (size_[1] - 2);

    x /= cellWidth;
    y /= cellHeight;

    // Als Test, Druck:

    x += offset_[0]; // p: 0.5, u: 0
    y += offset_[1]; // p: 0.5, u: 0.5

    assert(0 <= x && x < size_[0]);
    assert(0 <= y && y < size_[1]);

    const size_t i1 = static_cast<size_t>(x);
    const size_t i2 = static_cast<size_t>(x + 0.5);
    const double alpha = x - i1;

    const size_t j1 = static_cast<size_t>(y);
    const size_t j2 = static_cast<size_t>(y + 0.5);
    const double beta = y - j1;

    // TODO: Randwerte?
    const double d11 = data_.at(j1 * size_[0] + i1);
    const double d21 = data_.at(j1 * size_[0] + i2);
    const double d12 = data_.at(j2 * size_[0] + i1);
    const double d22 = data_.at(j2 * size_[0] + i2);

    const double dj1 = (1 - alpha) * d11 + alpha * d21;
    const double dj2 = (1 - alpha) * d12 + alpha * d22;

    return (1 - beta) * dj1 + beta * dj2;
}
