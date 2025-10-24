#pragma once

#include <settings.h>
#include <memory>
#include <vector>

struct Cell {
    double u;
    double v;
    double p;
};

struct Grid {

    std::vector<Cell> grid;

    const Settings &settings;

    Grid(const Settings &settings);

    Cell& operator()(int i, int j);

    double uInterpolateAt(int i, int j);
    double vInterpolateAt(int i, int j);
    double pInterpolateAt(int i, int j);

    double uDerivative(int i, int j);
    double vDerivative(int i, int j);
};
