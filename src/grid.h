#pragma once

#include "settings.h"
#include <memory>
#include <vector>

#include <cassert> // TODO: Make sure, asserts are not included in release build.

struct Cell {
    double u;
    double v;
    double p;
};

struct IntermediateCell {
    double F;
    double G;
};

struct VelocityEnv {
    double u[3];

    inline double operator()(const int i) const {
        assert(-2 < i && i < 2);
        return u[1 + i];
    }

    double derivative();
};

struct PressureEnv {
    double p[5]; //[ p_{i,j-1}, p_{i-1,j}, p_{i,j}, p_{i,j}, p_{i,j+1} ]

    inline double operator()(const int i, const int j) const {
        assert(-2 < i && i < 2 && -2 < j && j < 2 && std::abs(i) + std::abs(j) < 2);  // TODO: test all cases
        return p[(j + 1) * 2 + i];
    }

    double xDerivative();
    double yDerivative();
};

struct CellEnvironment {
    VelocityEnv u;
    VelocityEnv v;
    PressureEnv p;
};


struct Grid {

    std::vector<Cell> grid;
    std::vector<IntermediateCell> intermediateGrid; // TODO: like this or combine with cell?

    const Settings &settings;

    Grid(const Settings &settings);

    Cell& operator()(int i, int j);

    double uInterpolateAt(int i, int j);
    double vInterpolateAt(int i, int j);
    double pInterpolateAt(int i, int j);

    double uDerivative(int i, int j);
    double vDerivative(int i, int j);

    CellEnvironment get(int i, int j);
};
