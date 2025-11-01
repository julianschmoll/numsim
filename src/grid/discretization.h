#pragma once

#include <array>
#include <stdlib.h>
#include "grid/dataField.h"


class Discretization {
    const std::array<double, 2> meshWidth_;
    /// number of cells in x, y direction not including boundary
    const std::array<int, 2> nCells_;

protected:
    DataField u_;
    DataField v_;
    DataField p_;
    DataField f_;
    DataField g_;
    DataField rhs_;

public:
    virtual ~Discretization() = default;

    //! construct the object with given numbers of cells in x and y direction
    Discretization(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth);

    const std::array<double, 2> &meshWidth() const;
    const std::array<int, 2> &nCells() const;


    /// compute du²/dx
    virtual double computeDu2Dx(int i, int j) const = 0;

    /// compute dv²/dy
    virtual double computeDv2Dy(int i, int j) const = 0;

    /// compute duv/dx
    virtual double computeDuvDx(int i, int j) const = 0;

    /// compute duv/dy
    virtual double computeDuvDy(int i, int j) const = 0;


    /// compute d²u/dx²
    double computeD2uDx2(int i, int j) const;

    /// compute d²u/y²
    double computeD2uDy2(int i, int j) const;

    /// compute d²v/dx²
    double computeD2vDx2(int i, int j) const;

    /// compute d²v/y²
    double computeD2vDy2(int i, int j) const;

    /// compute dp/x
    double computeDpDx(int i, int j) const;

    /// compute du/x
    double computeDuDx(int i, int j) const;

    /// compute dv/y
    double computeDvDy(int i, int j) const;

    /// compute dp/y
    double computeDpDy(int i, int j) const;

    double u(int i, int j);
    double v(int i, int j);
    double p(int i, int j);
    double f(int i, int j);
    double g(int i, int j);
    double rhs(int i, int j);

    DataField& u();
    DataField& v();
    DataField& p();
    DataField& f();
    DataField& g();
    DataField& rhs();

    double dx() const;
    double dy() const;

};
