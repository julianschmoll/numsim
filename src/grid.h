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

    std::shared_ptr<Settings> settings;

    //Cell& operator (i, j)


    double uDerivative(int i, int j);
    double vDerivative(int i, int j);


};