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
