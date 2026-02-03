#pragma once

#include "../simulation/partitioning.h"
#include "grid/array2d.h"
#include "grid/dataField.h"
#include "settings.h"
#include <array>
#include <vector>

/**
 * @class StaggeredGrid
 * @brief Staggered Grid for numerical solving of poission problem.
 */
class StaggeredGrid {

    /// Cell types for structure representation
    enum CellType { Fluid = false, Solid = true };

    /// number of cells in x, y direction not including boundary
    const std::array<int, 2> nCells_;

    /// Mesh width of the grid
    const std::array<double, 2> meshWidth_;

    /// Reference to partitioning information
    const Partitioning &partitioning_;

protected:
    /// Field for Velocity in y-direction.
    DataField u_;

    /// Field for Velocity in x-direction.
    DataField v_;

    /// Field of the pressure.
    DataField p_;

    /// Field for preliminary velocity in x-direction.
    DataField f_;

    /// Field for preliminary velocity in y-direction.
    DataField g_;

    /// Field for storing rhs of the poission equation.
    DataField rhs_;

    /// Field for unphysical corrective pressure for solid cell movement correction
    DataField q_;

    /// Forces on the structure boundaries
    DataField fTop_;
    DataField fBottom_;

public:
    /// Positions of the structure boundaries
    std::vector<double> bottomBoundaryPosition_;
    std::vector<double> topBoundaryPosition_;

    /// Displacements of the structure boundaries
    std::vector<double> displacementsTop_;
    std::vector<double> displacementsBottom_;

    /// Field to represent structure presence in cells
    Array2d<bool> structure_;

    /**
     * Destructs Staggered Grid instance.
     */
    virtual ~StaggeredGrid() = default;

    /**
     * Constructs Staggered Grid instance.
     *
     * @param nCells Dimensions of the domain,
     * @param meshWidth Mesh width in x and y directions.
     * @param partitioning Object containing information on how the domain is partitioned.
     */
    StaggeredGrid(const Settings &settings, const Partitioning &partitioning);

    /**
     * Gets Mesh width.
     *
     * @return Mesh Width in x and y direction.
     */
    const std::array<double, 2> &meshWidth() const;

    /**
     * Gets the number of cells on Staggered Grid.
     *
     * @return Number of cells in x and y direction.
     */
    const std::array<int, 2> &nCells() const;

    /**
     * Gets the value of u at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double u(int i, int j);

    /**
     * Gets the value of v at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double v(int i, int j);

    /**
     * Gets the value of p at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double p(int i, int j);

    /**
     * Gets the value of q at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double q(int i, int j);

    /**
     * Gets the value of f at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double f(int i, int j);

    /**
     * Gets the value of g at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double g(int i, int j);

    /**
     * Gets the value of rhs at a specified point on grid.
     *
     * @param i i index of the point.
     * @param j j index of the point.
     * @return The value at specified position.
     */
    double rhs(int i, int j);

    /// DataField for velocity in x direction.
    DataField &u();

    /// DataField for velocity in y direction.
    DataField &v();

    /// DataField for pressure.
    DataField &p();

    DataField &q();

    /// DataField for preliminary velocity in x direction.
    DataField &f();

    /// DataField for preliminary velocity in y direction.
    DataField &g();

    /// DataField for rhs of the poisson problem.
    DataField &rhs();

    /// DataField for force at top boundary
    DataField &bottomF();

    /// DataField for force at bottom boundary
    DataField &topF();

    /// lower edge of cell in vertical physical coordinates
    double globalDomainPosJ(int j);

    /**
     * Getter for force at index i.
     * @param i Index to get force from.
     * @return the force at index i.
     */
    double &bottomF(int i);

    /**
     * Getter for bottom force at index i.
     * @param i Index to get force from.
     * @return the force at index i.
     */
    double &topF(int i);

    /**
     * Getter for displacement at index i.
     * @param i Index to get displacement from.
     * @return the displacement at index i.
     */
    double &bottomDisplacement(int i);

    /**
     * Getter for displacement at index i.
     * @param i Index to get displacement from.
     * @return the displacement at index i.
     */
    double &topDisplacement(int i);

    /**
     * Gets mesh width in x direction.
     *
     * @return Mesh width in x direction.
     */
    double dx() const;

    /**
     * Gets mesh width in y direction.
     *
     * @return Mesh width in y direction.
     */
    double dy() const;

    /**
     * Checks if cell at (i,j) is a fluid cell.
     *
     * @param i i index of the cell.
     * @param j j index of the cell.
     * @return true if cell is fluid, false otherwise.
     */
    bool isFluid(int i, int j) const;

    /**
     * Checks if cell at (i,j) is a solid cell.
     *
     * @param i i index of the cell.
     * @param j j index of the cell.
     * @return true if cell is solid, false otherwise.
     */
    bool isSolid(int i, int j) const;

    /**
     * Updates the structure cells based on the current boundary positions.
     *
     * @param dt Time step width.
     * @return true if structure cells were updated, false otherwise.
     */
    bool updateStructureCells(double dt);

    /**
     * Initializes the structure field based on initial boundary positions.
     */
    void initializeStructureField();

    /**
     * Getter for bottom boundary position at index i.
     * @param i Index to get position from.
     * @return the bottom boundary position at index i.
     */
    double &bottomBoundaryPosition(int i);

    /**
     * Getter for top boundary position at index i.
     * @param i Index to get position from.
     * @return the top boundary position at index i.
     */
    double &topBoundaryPosition(int i);

    /**
     * Test function for debugging.
     * @param settings Configuration settings for the simulation.
     */
    void test(const Settings &settings);
};
