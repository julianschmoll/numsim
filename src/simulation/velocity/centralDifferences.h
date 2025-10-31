#pragma once

#include "grid/discretization.h"

class centralDifferences final : public discretization {

public:
    centralDifferences(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth);

    double computeDu2Dx(int i, int j) const override;

    double computeDv2Dy(int i, int j) const override;

    double computeDuvDx(int i, int j) const override;

    double computeDuvDy(int i, int j) const override;

};
