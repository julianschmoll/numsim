#pragma once

#include "grid/discretization.h"

class centralDifferences final : public discretization {
    using discretization::discretization;

/*
    double computeDu2Dx(int i, int j) const override;

    double computeDu2Dy(int i, int j) const override;

    double computeDuvDx(int i, int j) const override;

    double computeDuvDy(int i, int j) const override;
*/
};
