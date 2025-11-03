#include "grid/array2d.h"

#include <algorithm>
#include <cassert>
#include <cmath>

array2d::array2d(const std::array<int, 2> size) : size_(size) { data_.resize(size_[0] * size_[1], 0.0); }

std::array<int, 2> array2d::size() const { return size_; }

double& array2d::operator()(const int i, const int j) {
    const int index = (j +1) * size_[0] + i +1;
    assert(index < static_cast<int>(data_.size()));
    return data_[index];
}

double array2d::operator()(const int i, const int j) const {
    const int index = (j +1) * size_[0] + i +1;
    assert(index < static_cast<int>(data_.size()));
    return data_[index];
}

double array2d::absMax() const {
    return std::fabs(
        *std::max_element(
            data_.begin(), data_.end(),
            [](const double a, const double b) { return std::fabs(a) < std::fabs(b); }));
}
