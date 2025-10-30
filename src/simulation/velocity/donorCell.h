#pragma once

#include "grid/discretization.h"

class donorCell final : public discretization {
    using discretization::discretization;

public:

    donorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha);

    double computeDu2Dx(int i, int j);

    double computeDv2Dy(int i, int j);

    double computeDuvDx(int i, int j);

    double computeDuvDy(int i, int j);

private:
    double alpha_;
};
