#include "grid/dataField.h"

#include <cassert>
#include <cmath>
#include <mpi.h>


DataField::DataField() : 
    meshWidth_({0, 0}),
    offset_({0, 0}),
    fieldID_(-1)
{}


DataField::DataField(const std::array<int, 2> size, const std::array<double, 2> meshWidth, const std::array<double, 2> offset, int fieldID)
    : Array2d(size), meshWidth_(meshWidth), offset_(offset), fieldID_(fieldID)
{}

int DataField::rows() const {
    return size_[1];
}

int DataField::cols() const {
    return size_[0];
}

int DataField::beginJ() const {
    return -1;
}

// TODO: end hat mich jetzt schon mehrfach verarscht...
int DataField::endJ() const {
    return size_[1] - 1; // size = cells + 2
}

int DataField::beginI() const {
    return -1;
}

int DataField::endI() const {
    return size_[0] - 1;
}

void DataField::setToZero() {
    for (auto &entry : data_) entry = 0;
}

int DataField::getID() const {
    return fieldID_;
}

double DataField::interpolateAt(double x, double y) const {
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
