#include "dataCell.h"
#include <assert.h>

dataCell::dataCell(const std::array<unsigned int, 2> size, const std::array<double, 2> meshWidth) : array2d(size),
    meshWidth_(meshWidth) {}

double dataCell::interpolateAt(double x, double y) const {
    const double interpolated = 1.0;

    x /= meshWidth_[0] * (size_[0] - 1); // x := x / h, h = meshWidth_[0] / (size_[0] - 1)
    y /= meshWidth_[1] * (size_[1] - 1);

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
