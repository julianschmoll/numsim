#pragma once

#include "grid/discretization.h"

class DonorCell final : public Discretization {

public:
    DonorCell(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, double alpha);

    double computeDu2Dx(int i, int j) const override;

    double computeDv2Dy(int i, int j) const override;

    double computeDuvDx(int i, int j) const override;

    double computeDuvDy(int i, int j) const override;

private:
    double alpha_;
};
