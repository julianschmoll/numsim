#pragma once

#include <array>
#include <iomanip>
#include <vector>
#include <cassert>
#include <iostream>

// TODO: This is because a std::vector<bool> specialization requires some serious template magic to overload the operators.
enum class CellType {
    Fluid, Solid
};

/**
 * @class Array2d
 * @brief Stores Data in 2D array, while internally storing values in a flat std::vector.
 */
template<typename T>
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
    T &operator()(int i, int j);

    /**
     * Returns the value of 2D Array at specified coordinates.
     *
     * Value is returned as const and cannot be overridden.
     *
     * @param i X-Coordinate of the value.
     * @param j Y-Coordinate of the value.
     * @return The value.
     */
    T operator()(int i, int j) const;

    /**
     * Returns the max value in the 2D Array using c++ standard library.
     *
     * Convenience method, so that iteration over both dimensions is not
     * necessary to retrieve max value.
     *
     * @return Max value.
     */
    [[nodiscard]] double absMax() const;

    T *data();

    template<typename S>
    friend std::ostream &operator<<(std::ostream& out, const Array2d<S> &arr);

protected:
    /// Flat storage for 2D array values
    std::vector<T> data_;

    /// Dimensions: width, height
    std::array<int, 2> size_;
};


template<typename T>
T Array2d<T>::operator()(int i, int j) const {
    const int index = (j + 1) * size_[0] + i + 1;
    assert(index < static_cast<int>(data_.size()));
    return data_[index];
}

template<typename T>
T &Array2d<T>::operator()(int i, int j) {
    const int index = (j + 1) * size_[0] + i + 1;
    assert(index < static_cast<int>(data_.size()));
    return data_[index];
}

template<typename T>
Array2d<T>::Array2d(const std::array<int, 2> size) : size_(size) {
    data_.resize(size_[0] * size_[1], T{});
}

template<typename T>
Array2d<T>::Array2d() : size_({0, 0}) {}

template<typename T>
std::array<int, 2> Array2d<T>::size() const {
    return size_;
}

template<typename T>
T *Array2d<T>::data() {
    return data_.data();
}

template<typename T>
std::ostream &operator<<(std::ostream& out, const Array2d<T> &arr) {
    for (int j = arr.size_[1] - 1; j >= 0; --j) {
        for (int i = 0; i < arr.size_[0]; ++i) {
            out << std::setw(6) << std::setprecision(4) << arr.data_[j * arr.size_[0] + i] << " ";
        }
        if (j > 0) out << "\n";
    }
    return out;
}

template<>
inline std::ostream &operator<<(std::ostream& out, const Array2d<CellType> &arr) {
    for (int j = arr.size_[1] - 1; j >= 0; --j) {
        for (int i = 0; i < arr.size_[0]; ++i) {
            out << (arr.data_[j * arr.size_[0] + i] == CellType::Solid ? "# " : ". ");
        }
        if (j > 0) out << "\n";
    }
    return out;
}


