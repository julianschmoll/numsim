#pragma once

#include "grid/array2d.h"

/**
 * @class DataField
 * @brief 2D Data representation containing method for interpolation.
 */
class DataField : public array2d {
public:
    DataField(std::array<int, 2> size, std::array<double, 2> meshWidth, std::array<double, 2> offset,
              std::array<int, 2> iIndexRange, std::array<int, 2> jIndexRange);

    /**
     * Interpolates the value at a given coordinate (x, y).
     *
     * Uses bilinear interpolation based on mesh width and offset.
     *
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @return Interpolated value at (x, y).
     */
    [[nodiscard]] double interpolateAt(double x, double y) const;

    /**
     * Returns the starting index in the j-direction.
     *
     * @return First valid j-index.
     */
    [[nodiscard]] int beginJ() const;

    /**
     * Returns the ending index in the j-direction.
     *
     * @return Last valid j-index.
     */
    [[nodiscard]] int endJ() const;

    /**
     * Returns the starting index in the i-direction.
     *
     * @return First valid i-index.
     */
    [[nodiscard]] int beginI() const;

    /**
     * Returns the ending index in the i-direction.
     *
     * @return Last valid i-index.
     */
    [[nodiscard]] int endI() const;

private:
    // Grid Spacing
    const std::array<double, 2> meshWidth_;

    /**
     * Describes the offset of the position of the data point from the lower left cell edge.
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

    // Valid index range in i direction
    const std::array<int, 2> iIndexRange_;

    // Valid index range in j direction
    const std::array<int, 2> jIndexRange_;
};
