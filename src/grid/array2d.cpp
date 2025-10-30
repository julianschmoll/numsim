#include "grid/array2d.h"

#include <algorithm>
#include <cassert>

array2d::array2d(const std::array<int, 2> size) : size_(size) { data_.resize(size_[0] * size_[1], 0.0); }

std::array<int, 2> array2d::size() const { return size_; }

double& array2d::operator()(const int i, const int j) {
    const int index = j * size_[0] + i;
    assert(index < static_cast<int>(data_.size()));
    return data_[index];
}

double array2d::operator()(const int i, const int j) const {
    const int index = j * size_[0] + i;
    assert(index < static_cast<int>(data_.size()));
    return data_[index];
}

// ToDo: Not sure if we actually need absMax or max. absmax could be retrieved with lambda function
/*
 *
 * #include <cmath>
 *
 * double array2d::absmax() const {
 *   return *std::max_element(data_.begin(), data_.end(),
 *       [](double a, double b) {
 *           return std::fabs(a) < std::fabs(b);
 *       });
 * }
 */
double array2d::max() const {
    return *std::max_element(data_.begin(), data_.end());
}
