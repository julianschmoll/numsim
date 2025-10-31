#pragma once

#include <memory>
#include <cmath>
#include <cassert>

#include "grid/array2d.h"


class dataField : public array2d {
public:
    dataField(std::array<int, 2> size, std::array<double, 2> meshWidth, std::array<double, 2> offset);

    double interpolateAt(double x, double y) const;

private:
    const std::array<double, 2> meshWidth_;
    const std::array<double, 2> offset_;
};
