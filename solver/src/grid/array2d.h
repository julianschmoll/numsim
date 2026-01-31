#pragma once

#include <array>
#include <iomanip>
#include <vector>
#include <cassert>
#include <iostream>
#include <macros.h>

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

    /// Array size in x direction. Equivalent to size()[0].
    int sizeI() const;
    /// Array size in y direction. Equivalent to size()[1].
    int sizeJ() const;

    /// Maximum valid index in x direction. Indexing begins at -1.
    int maxI() const;
    /// Minimum valid index in y direction. Indexing begins at -1.
    int maxJ() const;

    /// Minimum valid index in x direction. Indexing begins at -1.
    constexpr int minI() const;
    /// Minimum valid index in y direction. Indexing begins at -1.
    constexpr int minJ() const;

    /**
     * Returns the value of 2D Array at specified coordinates.
     *
     * Since this is not const, the value can be overridden.
     *
     * @param i X-Coordinate of value to return.
     * @param j Y-Coordinate of value to return.
     * @return Value
     */
    std::vector<T>::reference operator()(int i, int j) {
        DEBUG(checkIndices(i, j));
        return data_[ (j + 1) * size_[0] + i + 1];
    }

    /**
     * Returns the value of 2D Array at specified coordinates.
     *
     * Value is returned as const and cannot be overridden.
     *
     * @param i X-Coordinate of the value.
     * @param j Y-Coordinate of the value.
     * @return The value.
     */
    T operator()(int i, int j) const {
        DEBUG(checkIndices(i, j));
        return data_[ (j + 1) * size_[0] + i + 1];
    }

    /**
     * Returns the max value in the 2D Array using c++ standard library.
     *
     * Convenience method, so that iteration over both dimensions is not
     * necessary to retrieve max value.
     *
     * @return Max value.
     */
    [[nodiscard]] double absMax() const;

    /**
     * Set all array entries to `constant`.
     */
    void fill(const T& constant);

    T *data();

    template<typename S>
    friend std::ostream &operator<<(std::ostream& out, const Array2d<S> &arr);

protected:
    /// Flat storage for 2D array values
    std::vector<T> data_;

    /// Dimensions: width, height
    std::array<int, 2> size_;

private:
    void checkIndices(int i, int j) const;
};


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
void Array2d<T>::fill(const T &constant) {
    for (int j = 0; j < sizeJ(); ++j) {
        for (int i = 0; i < sizeI(); ++i) {
            data_[j * sizeI() + i] = constant;
        }
    }
}


template<typename T>
void Array2d<T>::checkIndices(int i, int j) const {
    assert(-1 <= i);
    assert(-1 <= j);
    assert(i < size_[0] - 1);
    assert(j < size_[1] - 1);
}

template<typename T>
inline std::ostream &operator<<(std::ostream& out, const std::vector<T> &arr) {
    out << "size=" << arr.size() << ": ";
    out << "[ ";
    for (double d : arr) out << d << " ";
    out << "]";
    return out;
}

template<typename T>
std::ostream &operator<<(std::ostream& out, const Array2d<T> &arr) {
    if (arr.sizeI() == 1 || arr.sizeJ() == 1) {
        out << arr.data_;
        return out;
    }
    for (int j = arr.sizeJ() - 1; j >= 0; --j) {
        for (int i = 0; i < arr.sizeI(); ++i) {
            out << std::setw(6) << std::setprecision(4) << arr.data_[j * arr.sizeI() + i] << " ";
        }
        if (j > 0) out << "\n";
    }
    return out;
}

template<>
inline std::ostream &operator<<(std::ostream& out, const Array2d<bool> &arr) {
    for (int j = arr.sizeJ() - 1; j >= 0; --j) {
        for (int i = 0; i < arr.sizeI(); ++i) {
            out << (arr.data_[j * arr.sizeI() + i] ? "# " : ". ");
        }
        if (j > 0) out << "\n";
    }
    return out;
}

template<typename T>
int Array2d<T>::sizeI() const { return size_[0]; };

template<typename T>
int Array2d<T>::sizeJ() const { return size_[1]; };

template<typename T>
int Array2d<T>::maxI() const { return size_[0] - 2; };

template<typename T>
int Array2d<T>::maxJ() const { return size_[1] - 2; };

template<typename T>
constexpr int Array2d<T>::minI() const { return -1; };

template<typename T>
constexpr int Array2d<T>::minJ() const { return -1; };