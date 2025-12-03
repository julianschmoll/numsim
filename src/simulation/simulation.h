#pragma once
#include "settings.h"
#include "partitioning.h"
#include "pressureSolver/redBlack.h"
#include "simulation/discreteOperators.h"
#include "outputWriter/outputWriter.h"

struct TimeSteppingInfo {
    double convectiveConstraint;
    double diffusiveConstraint;
    double maxVelocity;
    double timeStepWidth;
};

// ToDo: Do we really want/need inheritance here?
class Simulation  {
public:
    /**
     * Runs the simulation.
     */
    void run();

    Simulation(const Settings &settings);

private:
    // Grid width in x and y directions
    std::array<double, 2> meshWidth_{};

    // Discrete operators for grid data
    std::shared_ptr<DiscreteOperators> discOps_;

    // Writer for exporting simulation results for Paraview
    std::unique_ptr<OutputWriter> outputWriterParaview_;

    // Writer for exporting simulation results in plain text
    std::unique_ptr<OutputWriter> outputWriterText_;

    // Configuration settings for the simulation
    Settings settings_;

    // Time step size used in the simulation loop
    double timeStepWidth_ = 0.1;
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
    void setPreliminaryVelocities();
    void setRightHandSide();
    void setVelocities();

    void printConsoleInfo(double currentTime, const TimeSteppingInfo &timeSteppingInfo) const;

    std::shared_ptr<Partitioning> partitioning_;

    // Solver for the pressure
    std::unique_ptr<RedBlack> pressureSolver_;

};
