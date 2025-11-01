#pragma once

#include "grid/discretization.h"

class CentralDifferences final : public Discretization {

public:
    CentralDifferences(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth);

    double computeDv2Dy(int i, int j) const override;

    double computeDu2Dx(int i, int j) const override;

    double computeDuvDy(int i, int j) const override;

    double computeDuvDx(int i, int j) const override;
};
