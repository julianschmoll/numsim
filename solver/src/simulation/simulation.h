#pragma once
#include "grid/array2d.h"
#include "grid/dataField.h"
#include "outputWriter/outputWriter.h"
#include "partitioning.h"
#include "settings.h"
#include "simulation/discreteOperators.h"
#include "simulation/pressureSolver/pressureSolver.h"

#include <vector>

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
     * Initializes the displacement. Should be called before Simulation starts.
     *
     * @param displacements Flat vector of displacements for top and bottom structure boundary.
     */
    void initializeDisplacements(std::vector<double> &displacements);

    /**
     * Initializes the displacement. Should be called before Simulation starts.
     *
     * @param topDisplacements Vector of displacements for top structure boundary.
     * @param bottomDisplacements Vector of displacements for bottom structure boundary.
     */
    void initializeDisplacements(const std::vector<double> &topDisplacements, const std::vector<double> &bottomDisplacements);

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

    /**
     * Writes output files.
     *
     * @param currentSec Current simulation second.
     * @param lastSec Last simulation second when output was written.
     * @param always True if output should be written for every timestep.
     */
    void writeOutput(const int currentSec, const int lastSec, bool always = true) const;

    /**
     * Updates final velocities based on solved pressure.
     */
    void setVelocities();

    /**
     * Corrects velocities after boundary movement based on unphysical q pressure field.
     */
    void correctVelocities();

    /**
     * Calculates forces at top and bottom of the fluid domain.
     */
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

    /**
     * Updates the solid structure in the simulation.
     */
    void updateSolid();

    /**
     * Computes the time step width dt from maximum velocities.
     */
    TimeSteppingInfo computeTimeStepWidth();

    /**
     * Displacements must be set for every vertical partition.
     */
    void setDisplacements(const std::vector<double> &topDisplacements, const std::vector<double> &bottomDisplacements);

    /**
     * Test function for debugging.
     */
    void test();

    /**
     * Getter function for partitioning pointer to be used somewhere else.
     * @return Pointer to partitioning
     */
    std::shared_ptr<Partitioning> getPartitioning() const noexcept;

    /**
     * Getter function for discrete operators pointer to be used somewhere else.
     * @return Pointer to discrete operators
     */
    std::shared_ptr<DiscreteOperators> getDiscreteOperators() const noexcept {
        return discOps_;
    }

    /**
     * Gets the forces at the top and bottom boundaries.
     *
     * @param forces Flat vector to store the forces.
     */
    void getForces(std::vector<double> &forces);

    /**
     * Gets the displacements at the top and bottom boundaries.
     *
     * @param displacements Flat vector to store the displacements.
     */
    void setDisplacements(std::vector<double> &displacements);

    /**
     * Gets the current simulation time.
     *
     * @return Current simulation time.
     */
    double getCurrentTime();

    static void writeLineToFile(const std::string &filePath, const std::string &line);

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

    // Current time of the simulation.
    double currentTime_ = 0;

    // Old state of u to reload in precice
    DataField uCheckpoint_;

    // Old state of v to reload with precice
    DataField vCheckpoint_;

    // Old state of p to reload with precice
    DataField pCheckpoint_;

    // Old state of q to reload with precice
    DataField qCheckpoint_;

    // Old state of f to reload with precice
    DataField fCheckpoint_;

    // Old state of g to reload with precice
    DataField gCheckpoint_;

    // Time step width at checkpoint
    double timeStepWidthCheckpoint_ = 0;

    // Structure boundary positions at checkpoint
    std::vector<double> topBoundaryPositionCheckpoint_;
    std::vector<double> bottomBoundaryPositionCheckpoint_;

    // Displacements at checkpoint
    std::vector<double> displacementsTopCheckpoint_;
    std::vector<double> displacementsBottomCheckpoint_;

    // Structure presence at checkpoint
    Array2d<bool> structureCheckpoint_;

    // Time at checkpoint
    double checkpointTime_ = 0;

    /**
     * Sets boundary values of u, v, f and g.
     */
    void setOuterVelocityBoundaries();

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
    void printConsoleInfo(const TimeSteppingInfo &timeSteppingInfo) const;

    // Partitioning for the grid on multiple ranks.
    std::shared_ptr<Partitioning> partitioning_;

    // Solver for the pressure
    std::unique_ptr<PressureSolver> pressureSolver_;
};
