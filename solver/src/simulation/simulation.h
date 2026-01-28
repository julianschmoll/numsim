#pragma once
#include "outputWriter/outputWriter.h"
#include "partitioning.h"
#include "settings.h"
#include "simulation/discreteOperators.h"
#include "simulation/pressureSolver/pressureSolver.h"
#include "grid/dataField.h"

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
     * Saves current state of u, v and p.
     */
    void saveState();

    /**
     * Reloads saved states of u,v and p.
     */
    void reloadLastState();

    /**
     * Constructs simulation object.
     *
     * @param settings Settings to run simulation with.
     */
    explicit Simulation(const Settings &settings, const std::string &folderName);

    void writeOutput(double currentTime, int currentSec, int lastSec) const;

    /**
     * Updates final velocities based on solved pressure.
     */
    void setVelocities();

    /**
     * Corrects velocities after boundary movement based on unphysical q pressure field.
     */
    void correctVelocities();

    void calculateForces();

    /**
     * no slip at top and bottom structure boundary
     * call after other boundary methods
     */
    void setStructureBoundaries();

    /**
     * Sets the simulation timestep width (e.g. to the time step calculated by the precice adapter)
     * @param dt Timestep width
     */
    void setTimeStepWidth(double dt);

    // preCICE Interface:

    /**
     * set fluid boundaries, solve for preliminary velocities, pressure and then velocities
     * calculate forces
     */
    void advanceFluidSolver(double dt);

    void updateSolid();
    
    /**
     * Computes the time step width dt from maximum velocities.
     */
    TimeSteppingInfo computeTimeStepWidth(double currentTime);

    /**
     * Displacements must be set for every vertical partition.
     */
    void setDisplacements(const std::vector<double> &topDisplacements, const std::vector<double> &bottomDisplacements);

    void test();

    std::shared_ptr<Partitioning> getPartitioning() const noexcept;

    std::shared_ptr<DiscreteOperators> getDiscreteOperators() const noexcept {
        return discOps_;
    }

    
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

    // Old state of u to reload in precice
    DataField uCheckpoint_;

    // Old state of v to reload with precice
    DataField vCheckpoint_;

    // Old state of p to reload with precice
    DataField pCheckpoint_;

    /**
     * Sets boundary values of u and v.
     */
    void setBoundaryUV(double currentTime);

    /**
     * Sets boundary values of F and G.
     */
    void setBoundaryFG();


    /**
     * Sets preliminary velocities.
     */
    void setPreliminaryVelocities();

    /**
     * Sets rhs of the poisson equation to solve pressure with.
     */
    void setRightHandSide();

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
