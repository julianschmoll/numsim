#include "array2d.h"

#include <cassert>

array2d::array2d(std::array<int,2> size): size_(size)
{
    //allocate data and initialize to 0
    assert(size_[0] > 0 && size_[1] > 0);
    data_.resize(size_[0] * size_[1], 0.0);
}

std::array<int,2> array2d::size() const
{
    return size_;
}

double &array2d::operator()(int i, int j)
{
    const int index = j * size_[0] + i;

    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j*size_[0] + i < static_cast<int>(data_.size());

    return data_[index];
}
