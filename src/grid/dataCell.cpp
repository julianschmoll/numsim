#include "dataCell.h"


dataCell::dataCell(const std::array<unsigned int, 2> size, const std::array<double, 2> meshWidth) : array2d(size),
    meshWidth_(meshWidth) {}

double dataCell::interpolateAt(double x, double y) const {
    const double interpolated = 1.0;
    return interpolated;
}
