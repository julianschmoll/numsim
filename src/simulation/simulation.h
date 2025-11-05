#pragma once

#include "settings.h"
#include "grid/discretization.h"
#include "outputWriter/outputWriterText.h"
#include "outputWriter/outputWriterParaview.h"
#include "simulation/pressureSolver.h"
#include "simulation/discreteOperators.h"
#include <memory>

class Simulation {
public:
    explicit Simulation(const Settings& settings);
    ~Simulation() = default;

    /**
     * Runs the simulation.
     *
     * @returns Status code whether simulation was successful
     */
    int run();

private:
    /**
     * Sets boundary values of u and v.
     */
    void setBoundaryValues();

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

    std::array<double, 2> meshWidth_{};
    std::shared_ptr<discreteOperators> discOps_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    Settings settings_;
    double timeStepWidth_ = 0.1;
};
