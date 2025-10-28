#pragma once

#include <array>

#include "staggeredGrid.h"

class discretization : public staggeredGrid {
public:
    virtual ~discretization() = default;

    //! construct the object with given numbers of cells in x and y direction
    discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth);

    //! compute 1st derivative p / x
    virtual double computeDpDx(int i, int j) const = 0;

    //! compute 2nd derivative u / x^2
    virtual double computeD2uDx2(int i, int j) const = 0;

    //! compute 2nd derivative u / y^2
    virtual double compute2DuDy2(int i, int j) const = 0;

    //! compute the 1st derivative u^2 / x
    virtual double computeDu2Dx(int i, int j) const = 0;

    //! compute the 1st derivative (uv) / y
    virtual double computeDuvDy(int i, int j) const = 0;

    //! compute the 1st derivative u / x
    virtual double computeDuDx(int i, int j) const = 0;

    //! compute the 1st derivative v / y
    virtual double computeDvDy(int i, int j) const = 0;


    //! compute the 1st derivative v^2 / x
    virtual double computeDu2Dy(int i, int j) const = 0;

    //! compute the 1st derivative (uv) / x
    virtual double computeDuvDx(int i, int j) const = 0;

    //! compute 1st derivative p / y
    virtual double computeDpDy(int i, int j) const = 0;

    /** missing methods:
     * meshWidth
     * nCells
     * p
     * u
     * v
     * interpolateAt
    */
};
