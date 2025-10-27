#pragma once

#include "settings.h"

class Solver {

public:
    double time = 0.0;
    Settings settings_;

    Solver(const Settings &settings);
    ~Solver() = default;

    double CalculateTimeStepWidth();

    void setBoundaryValues();

    int AdvanceTimeStep();

};
