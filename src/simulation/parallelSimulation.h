#pragma once
#include "simulation.h"
#include "partitioning.h"
#include "outputWriter/outputWriterParaviewParallel.h"
#include "outputWriter/outputWriterTextParallel.h"

class ParallelSimulation final : public Simulation {
public:
    /**
     * Runs the simulation.
     */
    void run();

    using Simulation::Simulation;

    void initialize(const Settings &settings);

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
    void computeTimeStepWidth();

    std::shared_ptr<Partitioning> partitioning_;


};
