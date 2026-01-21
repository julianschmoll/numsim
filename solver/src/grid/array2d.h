#pragma once

#include <array>
#include <vector>

/**
 * @class Array2d
 * @brief Stores Data in 2D array, while internally storing values in a flat std::vector.
 */
class Array2d {
public:
    /**
     * Constructor for 2D Array storing data in a vector.
     *
     * This class allows indexed access in two dimensions,
     * while simply storing the data in a std::array.
     *
     * @param size Dimensions of the 2d array
     */
    explicit Array2d(std::array<int, 2> size);

    /// initialize size with 0 and data_ via the default constructor
    Array2d();

    /**
     * Destructor for Array2d.
     */
    virtual ~Array2d() = default;

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
    double &operator()(int i, int j);

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

    /**
     * Returns the max value in the 2D Array using c++ standard library.
     *
     * Convenience method, so that iteration over both dimensions is not
     * necessary to retrieve max value.
     *
     * @return Max value.
     */
    [[nodiscard]] double absMax() const;

    double *data();

protected:
    /// Flat storage for 2D array values
    std::vector<double> data_;

    /// Dimensions: width, height
    std::array<int, 2> size_;
};
