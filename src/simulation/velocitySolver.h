#pragma once

#include "grid/discretization.h"

class velocitySolver final : public Discretization {

public:
    velocitySolver(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, double alpha);

    double computeDu2Dx(const int i, const int j) const;

    double computeDv2Dy(const int i, const int j) const;

    double computeDuvDx(const int i, const int j) const;

    double computeDuvDy(const int i, const int j) const;

    /// compute d²u/dx²
    double computeD2uDx2(const int i, const int j) const;

    /// compute d²u/y²
    double computeD2uDy2(const int i, const int j) const;

    /// compute d²v/dx²
    double computeD2vDx2(const int i, const int j) const;

    /// compute d²v/y²
    double computeD2vDy2(const int i, const int j) const;

    /// compute dp/x
    double computeDpDx(const int i, const int j) const;

    /// compute du/x
    double computeDuDx(const int i, const int j) const;

    /// compute dv/y
    double computeDvDy(const int i, const int j) const;

    /// compute dp/y
    double computeDpDy(const int i, const int j) const;

private:
    double alpha_;
};
