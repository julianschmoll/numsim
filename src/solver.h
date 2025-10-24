#pragma once
//#include <grid.h>
#include "settings.h"
class Solver {

public:
    double time = 0.0;
    Settings settings_;

    Solver(Settings settings) : settings_(settings) {};
    ~Solver() = default;

    int AdvanceTimeStep();

    //Grid grid;
};