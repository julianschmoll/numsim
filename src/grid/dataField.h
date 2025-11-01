#pragma once

#include <memory>
#include <cmath>
#include <cassert>

#include "grid/array2d.h"


class dataField : public array2d {
public:
    dataField(std::array<int, 2> size, std::array<double, 2> meshWidth, std::array<double, 2> offset, std::array<int, 2> iIndexRange, std::array<int, 2> jIndexRange);

    double interpolateAt(double x, double y) const;

    int beginJ() const;
    int endJ() const;
    int beginI() const;
    int endI() const;

private:
    const std::array<double, 2> meshWidth_;

    /**
     * Describes the offset of the position of the data point from the lower left cell edge.
     * TODO: unintuitive
     * 
     * (0,1)---(1,1)
     *   |       |
     *   |       |
     * (0,0)---(1,0)
     * 
     * - Pressure:                0.5, 0.5
     * - Velocity in x direction: 0,   0.5
     * - Velocity in y direction: 0.5, 0
     */
    const std::array<double, 2> offset_;

    const std::array<int, 2> iIndexRange_;
    const std::array<int, 2> jIndexRange_;
};
