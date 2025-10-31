#include "grid/discretization.h"

discretization::discretization(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth)
    : meshWidth_(meshWidth),
    nCells_(nCells),
    u_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.0, 0.5}),
    v_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.0}),
    p_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}),
    f_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.0, 0.5}),
    g_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.0}),
    rhs_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}) {
}

const std::array<double, 2> &discretization::meshWidth() const {
    return meshWidth_;
}

const std::array<int, 2> &discretization::nCells() const {
    return nCells_;
}

int discretization::pIBegin() const {
    return 1;
}

int discretization::pJBegin() const {
    return 1;
}

int discretization::pIEnd() const {
    return nCells_[0] + 1;
}
int discretization::pJEnd() const {
    return nCells_[1] + 1;
}

int discretization::vIBegin() const {
    return 1;
}
int discretization::vJBegin() const {
    return 1;
}
int discretization::vIEnd() const {
    return nCells_[0] + 1;
}
int discretization::vJEnd() const {
    return nCells_[1];
};

int discretization::uIBegin() const {
    return 1;
}
int discretization::uJBegin() const {
    return 1;
}
int discretization::uIEnd() const {
    return nCells_[0];
}
int discretization::uJEnd() const {
    return nCells_[1] + 1;
}

double discretization::u(int i, int j) {
    return u_(i,j);
}

double discretization::v(int i, int j) {
    return v_(i,j);
}

double discretization::p(int i, int j) {
    return p_(i,j);
}

double discretization::f(int i, int j) {
    return f_(i,j);
}

double discretization::g(int i, int j) {
    return g_(i,j);
}

double discretization::rhs(int i, int j) {
    return rhs_(i,j);
}

dataField &discretization::u() {
    return u_;
}

dataField& discretization::v() {
    return v_;
}

dataField& discretization::p() {
    return p_;
}

dataField& discretization::f() {
    return f_;
}

dataField& discretization::g() {
    return g_;
}

dataField& discretization::rhs() {
    return rhs_;
}


double discretization::dx() const {
    return meshWidth_[0];
}

double discretization::dy() const {
    return meshWidth_[1];
}



double discretization::computeD2uDx2(int i, int j) const {
    const double du2 = u_(i + 1, j) - 2 * u_(i, j) + u_(i - 1, j);
    return du2 / (dx() * dx());
}

double discretization::compute2DuDy2(int i, int j) const {
    const double du2 = u_(i, j + 1) - 2 * u_(i, j) + u_(i, j - 1);
    return du2 / (dy() * dy());
}

double discretization::computeDuDx(int i, int j) const {
    const double du = u_(i, j) - u_(i - 1, j);
    return du / dx();
}

double discretization::computeDvDy(int i, int j) const {
    const double dv = v_(i, j) - v_(i, j - 1);
    return dv / dy();
}

double discretization::computeDpDx(int i, int j) const {
    const double dp = p_(i + 1, j) - p_(i, j);
    return dp / dx();
}

double discretization::computeDpDy(int i, int j) const {
    const double dp = p_(i, j + 1) - p_(i, j);
    return dp / dy();
}
