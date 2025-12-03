#pragma once
#include "simulation.h"
#include "partitioning.h"
#include "pressureSolver/redBlack.h"

struct TimeSteppingInfo {
    double convectiveConstraint;
    double diffusiveConstraint;
    double maxVelocity;
    double timeStepWidth;
};

// ToDo: Do we really want/need inheritance here?
class ParallelSimulation final : public Simulation {
public:
    /**
     * Runs the simulation.
     */
    void run();

    using Simulation::Simulation;

    ParallelSimulation(const Settings &settings);

private:

    /**
    * Sets boundary values of u and v.
     */
    void setBoundaryUV();

    /**
     * Sets boundary values of F and G.
     */
    void setBoundaryFG();

    /**
    * Computes the time step width dt from maximum velocities.
    */
    TimeSteppingInfo computeTimeStepWidth(double currentTime);

    void printConsoleInfo(double currentTime, const TimeSteppingInfo &timeSteppingInfo) const;

    std::shared_ptr<Partitioning> partitioning_;

    // Solver for the pressure
    std::unique_ptr<RedBlack> pressureSolver_;

};
