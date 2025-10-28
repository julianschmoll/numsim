#include "grid/array2d.h"

#include <cassert>

array2d::array2d(const std::array<unsigned int, 2> size) : size_(size) { data_.resize(size_[0] * size_[1], 0.0); }

std::array<unsigned int, 2> array2d::size() const { return size_; }

double& array2d::operator()(const unsigned int x, const unsigned int y) {
    const unsigned int index = y * size_[0] + x;
    assert(index < static_cast<unsigned int>(data_.size()));
    return data_[index];
}

double array2d::operator()(const unsigned int x, const unsigned int y) const {
    const unsigned int index = y * size_[0] + x;
    assert(index < static_cast<unsigned int>(data_.size()));
    return data_[index];
}
