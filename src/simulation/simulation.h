#pragma once

#include "grid/staggeredGrid.h"
#include "outputWriter/outputWriterParaview.h"
#include "outputWriter/outputWriterText.h"
#include "settings.h"
#include "simulation/discreteOperators.h"
#include "simulation/pressureSolver.h"
#include <memory>

/**
 * @class Simulation
 * @brief Class responsible for running the simulation.
 *
 * This class contains the simulation loop and calls all the necessary
 * functions to compute velocities and pressure.
 *
 */
class Simulation {
public:
    /**
     * Constructs a Simulation object with given Settings.
     *
     * @param settings Configuration parameters for simulation.
     */
    explicit Simulation(Settings settings);

    /**
     * Destructs Simulation object.
     */
    ~Simulation() = default;

    /**
     * Runs the simulation.
     */
    void run();

private:
    // Grid width in x and y directions
    std::array<double, 2> meshWidth_{};

    // Discrete operators for grid data
    std::shared_ptr<DiscreteOperators> discOps_;

    // Writer for exporting simulation results for Paraview
    std::unique_ptr<OutputWriter> outputWriterParaview_;

    // Writer for exporting simulation results in plain text
    std::unique_ptr<OutputWriter> outputWriterText_;

    // Solver for the pressure
    std::unique_ptr<PressureSolver> pressureSolver_;

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
     * Computes the preliminary velocities, F and G.
     */
    void setPreliminaryVelocities();

    /**
     * Computes the right hand side of the Poisson equation.
     */
    void setRightHandSide();

    /**
     * Computes the time step width dt from maximum velocities.
     */
    void computeTimeStepWidth();

    /**
     * computes the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p.
     */
    void setVelocities();
};
