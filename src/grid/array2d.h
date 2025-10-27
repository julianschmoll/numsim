#pragma once

#include <vector>
#include <array>


class array2d
{
public:
    //! constructor arguments: width, heigth
    array2d(std::array<int, 2> size);

    //! returns the size <width, height>
    std::array<int, 2> size() const;

    //! access changeable
    double &operator()(int i, int j);

    //! access const
    double &operator()(int i, int j) const;

protected:
    std::vector<double> data_; //row major array values
    const std::array<int,2> size_; //width, height
};
