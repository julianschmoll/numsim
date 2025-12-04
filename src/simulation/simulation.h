#pragma once
#include "outputWriter/outputWriter.h"
#include "partitioning.h"
#include "settings.h"
#include "simulation/discreteOperators.h"
#include "simulation/pressureSolver/pressureSolver.h"

/**
 * @struct TimeSteppingInfo
 * @brief Holds all necessary data for calculating a stable time step.
 */
struct TimeSteppingInfo {
    double convectiveConstraint;
    double diffusiveConstraint;
    double maxVelocity;
    double timeStepWidth;
};

/**
 * @class Simulation
 * @brief Entry point for fluid simulation.
 */
class Simulation {
public:
    /**
     * Runs the simulation.
     */
    void run();

    /**
     * Constructs simulation object.
     *
     * @param settings Settings to run simulation with.
     */
    explicit Simulation(const Settings &settings);

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

    /**
     * Sets preliminary velocities.
     */
    void setPreliminaryVelocities();

    /**
     * Sets rhs of the poisson equation to solve pressure with.
     */
    void setRightHandSide();

    /**
     * Updates final velocities based on solved pressure.
     */
    void setVelocities();

    /**
     * Prints current progress of the simulation.
     *
     * @param currentTime Current time of the simulation.
     * @param timeSteppingInfo Struct storing time stepping information.
     */
    void printConsoleInfo(double currentTime, const TimeSteppingInfo &timeSteppingInfo) const;

    // Partitioning for the grid on multiple ranks.
    std::shared_ptr<Partitioning> partitioning_;

    // Solver for the pressure
    std::unique_ptr<PressureSolver> pressureSolver_;
};
