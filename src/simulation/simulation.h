#pragma once

#include "settings.h"
#include "grid/discretization.h"
#include "outputWriter/outputWriterText.h"
#include "outputWriter/outputWriterParaview.h"
#include "simulation/pressure/solver.h"
#include "simulation/velocity/centralDifferences.h"
#include "simulation/velocity/donorCell.h"
#include "simulation/pressure/gaussSeidel.h"
#include "simulation/pressure/sor.h"

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
    void computePreliminaryVelocities();

    /**
     * Solves the Poisson equation for the pressure.
     */
    void computePressure();

    /**
     * Computes the right hand side of the Poisson equation.
     */
    void computeRightHandSide();

    /**
     * Computes the time step width dt from maximum velocities.
     *
     * @return Time step width.
     */
    double computeTimeStepWidth() const;

    /**
     * computes the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p.
     */
    void computeVelocities();

    std::shared_ptr<discretization> discretization_;
    double meshWidth_ = 0.0;
    std::unique_ptr<outputWriterParaview> outputWriterParaview_;
    std::unique_ptr<outputWriterText> outputWriterText_;
    std::unique_ptr<solver> solver_;
    Settings settings_;
};
