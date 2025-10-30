#pragma once

#include "grid/discretization.h"

class centralDifferences final : public discretization {
    using discretization::discretization;

public:

    centralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth);

    double computeDu2Dx(int i, int j);

    double computeDv2Dy(int i, int j);

    double computeDuvDx(int i, int j);

    double computeDuvDy(int i, int j);
};
