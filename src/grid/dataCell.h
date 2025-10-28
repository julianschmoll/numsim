#pragma once

#include <memory>
#include <cmath>
#include <cassert>

#include "grid/array2d.h"


class dataCell :
    public array2d {
public:
    dataCell(std::array<unsigned int, 2> size, std::array<double, 2> meshWidth);

    double interpolateAt(double x, double y) const;

private:
    const std::array<double, 2> meshWidth_;
};
