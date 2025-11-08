#pragma once

#include "grid/discretization.h"

class DiscreteOperators final : public Discretization {

public:
    /**
     * Constructs a discreteOperators object that provides finite-difference approximations of
     * spatial derivatives for the grid data.
     *
     * @param nCells: Number of cells in x- and y-direction
     * @param meshWidth: Physical size of a grid cell in x and y-direction
     * @param alpha: Donor cell contribution factor, range from 0 to 1. If set 0, pure central
     * differences discretization is used.
     */
    DiscreteOperators(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, double alpha);

    /**
     * Computes du²/dx
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of u² in x-direction
     */
    double computeDu2Dx(int i, int j) const;

    /**
     * Computes dv²/dx
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of v² in y-direction
     */
    double computeDv2Dy(int i, int j) const;

    /**
     * Computes d(uv)/dx
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of uv in x-direction
     */
    double computeDuvDx(int i, int j) const;

    /**
     * Computes d(uv)/dx
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of uv in y-direction
     */
    double computeDuvDy(int i, int j) const;

    /**
     * Computes d²u²/dx²
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Second derivative of u in x-direction
     */
    double computeD2uDx2(int i, int j) const;

    /**
     * Computes d²u/dy²
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Second derivative of u in y-direction
     */
    double computeD2uDy2(int i, int j) const;

    /**
     * Computes d²v/dx²
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Second derivative of v in x-direction
     */
    double computeD2vDx2(int i, int j) const;

    /**
     * Computes d²v/dx²
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Second derivative of v in x-direction
     */
    double computeD2vDy2(int i, int j) const;

    /**
     * Computes dp/dx
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of p in x-direction
     */
    double computeDpDx(int i, int j) const;

    /**
     * Computes du/dx
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of u in x-direction
     */
    double computeDuDx(int i, int j) const;

    /**
     * Computes dv/dy
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return Derivative of v in y-direction
     */
    double computeDvDy(int i, int j) const;

    /**
     * Computes dp/dy
     *
     * @param i Index in x-direction
     * @param j Index in y-direction
     * @return
     */
    double computeDpDy(int i, int j) const;

private:
    double alpha_;
};
