#pragma once

#include "../simulation/partitioning.h"
#include "grid/array2d.h"
#include "grid/dataField.h"
#include <array>
#include <vector>

/**
 * @class StaggeredGrid
 * @brief Staggered Grid for numerical solving of poission problem.
 */
class StaggeredGrid {

    enum CellType {
        Fluid = false, Solid = true
    };

    /// Mesh width of the grid
    const std::array<double, 2> meshWidth_;

    /// number of cells in x, y direction not including boundary
    const std::array<int, 2> nCells_;

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

    Array2d<bool> structure_;

    // TODO: Bessere Möglichkeit? Wir brauchen die Zeitableitung...
    std::vector<double> oldDisplacementsTop_;
    std::vector<double> oldDisplacementsBottom_;

    std::vector<double> newDisplacementsTop_;
    std::vector<double> newDisplacementsBottom_;

    std::vector<double> fTop_;
    std::vector<double> fBottom_;

public:
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
    StaggeredGrid(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, const Partitioning &partitioning);

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

    // for accessing with the same indexing (-1, 0,..., size() - 2)
    double &fTop(int i);
    // for accessing with the same indexing (-1, 0,..., size() - 2)
    double &fBottom(int i);

    std::vector<double> &fTop();
    std::vector<double> &fBottom();

    /// DataField for velocity in x direction.
    DataField &u();

    /// DataField for velocity in y direction.
    DataField &v();

    /// DataField for pressure.
    DataField &p();

    /// DataField for preliminary velocity in x direction.
    DataField &f();

    /// DataField for preliminary velocity in y direction.
    DataField &g();

    /// DataField for rhs of the poisson problem.
    DataField &rhs();

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

    bool isFluid(int i, int j) const;
    bool isSolid(int i, int j) const;
};
