#pragma once

#include <array>
#include <stdlib.h>
#include "grid/dataField.h"

class discretization {
    const std::array<double, 2> meshWidth_;
    /// number of cells in x, y direction not including boundary
    const std::array<int, 2> nCells_;

protected:
    dataField u_;
    dataField v_;
    dataField p_;
    dataField f_;
    dataField g_;
    dataField rhs_;

public:
    virtual ~discretization() = default;

    //! construct the object with given numbers of cells in x and y direction
    discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth);

    const std::array<double, 2> &meshWidth() const;
    const std::array<int, 2> &nCells() const;

    //! compute the 1st derivative u^2 / x
    virtual double computeDu2Dx(int i, int j) const = 0;

    //! compute the 1st derivative v^2 / x
    virtual double computeDu2Dy(int i, int j) const = 0;

    //! compute the 1st derivative (uv) / x
    virtual double computeDuvDx(int i, int j) const = 0;

    //! compute the 1st derivative (uv) / y
    virtual double computeDuvDy(int i, int j) const = 0;

    //! compute 2nd derivative u / x^2
    double computeD2uDx2(int i, int j);

    //! compute 1st derivative p / x
    double computeDpDx(int i, int j);

    //! compute 2nd derivative u / y^2
    double compute2DuDy2(int i, int j);

    //! compute the 1st derivative u / x
    double computeDuDx(int i, int j);

    //! compute the 1st derivative v / y
    double computeDvDy(int i, int j);

    //! compute 1st derivative p / y
    double computeDpDy(int i, int j);

    int pIBegin() const;
    int pJBegin() const;
    int pIEnd() const;
    int pJEnd() const;

    int vIBegin() const;
    int vJBegin() const;
    int vIEnd() const;
    int vJEnd() const;

    int uIBegin() const;
    int uJBegin() const;
    int uIEnd() const;
    int uJEnd() const;

    double u(int i, int j);
    double v(int i, int j);
    double p(int i, int j);
    double f(int i, int j);
    double g(int i, int j);
    double rhs(int i, int j);

    dataField& u();
    dataField& v();
    dataField& p();
    dataField& f();
    dataField& g();
    dataField& rhs();

    double dx() const;
    double dy() const;

};
