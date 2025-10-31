#pragma once

#include <vector>
#include <array>


class array2d {
public:
    /**
     * Constructor for 2D Array storing data in a vector.
     *
     * This class allows indexed access in two dimensions,
     * while simply storing the data in a std::array.
     *
     * @param size Dimensions of the 2d array
     */
    explicit array2d(std::array<int, 2> size);

    /**
     * Gets the size of the 2D Array in x any y direction.
     *
     * @return Dimensions of the 2d Array.
     */
    [[nodiscard]] std::array<int, 2> size() const;

    /**
     * Returns the value of 2D Array at specified coordinates.
     *
     * Since this is not const, the value can be overridden.
     *
     * @param i X-Coordinate of value to return.
     * @param j Y-Coordinate of value to return.
     * @return Value
     */
    double& operator()(int i, int j);

    /**
     * Returns the value of 2D Array at specified coordinates.
     *
     * Value is returned as const and cannot be overridden.
     *
     * @param i X-Coordinate of the value.
     * @param j Y-Coordinate of the value.
     * @return The value.
     */
    double operator()(int i, int j) const;

protected:
    std::vector<double> data_;
    std::array<int, 2> size_;
};
