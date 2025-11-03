#pragma once

#include "grid/discretization.h"

class velocitySolver final : public Discretization {

public:
    velocitySolver(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, double alpha);

    double computeDu2Dx(int i, int j) const;

    double computeDv2Dy(int i, int j) const;

    double computeDuvDx(int i, int j) const;

    double computeDuvDy(int i, int j) const;

    /// compute d²u/dx²
    double computeD2uDx2(int i, int j) const;

    /// compute d²u/y²
    double computeD2uDy2(int i, int j) const;

    /// compute d²v/dx²
    double computeD2vDx2(int i, int j) const;

    /// compute d²v/y²
    double computeD2vDy2(int i, int j) const;

    /// compute dp/x
    double computeDpDx(int i, int j) const;

    /// compute du/x
    double computeDuDx(int i, int j) const;

    /// compute dv/y
    double computeDvDy(int i, int j) const;

    /// compute dp/y
    double computeDpDy(int i, int j) const;

private:
    double alpha_;
};
