#pragma once

#include "simulation/velocity/centralDifferences.h"
#include "settings.h"
#include "grid/discretization.h"
# include "simulation/pressure/solver.h"

class Simulation {
public:

    Simulation(Settings settings);
    ~Simulation() = default;

    void run();

private:
    // double CalculateTimeStepWidth();

    // void setBoundaryValues();

    std::shared_ptr<discretization> discretization_;
    std::shared_ptr<solver> solver_;
    Settings settings_;
    double timeStepWidth_ = 0.0;
};
