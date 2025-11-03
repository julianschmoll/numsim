#include "grid/discretization.h"

Discretization::Discretization(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth)
    : meshWidth_(meshWidth),
      nCells_(nCells),
      u_({nCells[0] + 1, nCells[1] + 2}, meshWidth, {0.0, 0.5}, {-1, nCells_[0]}, {-1, nCells_[1] +1}),
      v_({nCells[0] + 2, nCells[1] + 1}, meshWidth, {0.5, 0.0}, {-1, nCells_[0] +1}, {-1, nCells_[1]}),
      p_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}, {-1, nCells_[0] +1}, {-1, nCells_[1]+1}),
      f_({nCells[0] + 1, nCells[1] + 2}, meshWidth, {0.0, 0.5}, {-1, nCells_[0]}, {-1, nCells_[1] +1}),
      g_({nCells[0] + 2, nCells[1] + 1}, meshWidth, {0.5, 0.0}, {-1, nCells_[0] +1}, {-1, nCells_[1]}),
      rhs_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}, {-1, nCells_[0] +1}, {-1, nCells_[1]+1}) {
}

const std::array<double, 2> &Discretization::meshWidth() const {
    return meshWidth_;
}

const std::array<int, 2> &Discretization::nCells() const {
    return nCells_;
}

double Discretization::u(int i, int j) {
    return u_(i,j);
}

double Discretization::v(int i, int j) {
    return v_(i,j);
}

double Discretization::p(int i, int j) {
    return p_(i,j);
}

double Discretization::f(int i, int j) {
    return f_(i,j);
}

double Discretization::g(int i, int j) {
    return g_(i,j);
}

double Discretization::rhs(int i, int j) {
    return rhs_(i,j);
}

DataField &Discretization::u() {
    return u_;
}

DataField& Discretization::v() {
    return v_;
}

DataField& Discretization::p() {
    return p_;
}

DataField& Discretization::f() {
    return f_;
}

DataField& Discretization::g() {
    return g_;
}

DataField& Discretization::rhs() {
    return rhs_;
}


double Discretization::dx() const {
    return meshWidth_[0];
}

double Discretization::dy() const {
    return meshWidth_[1];
}

double Discretization::computeD2uDx2(int i, int j) const {
    const double du2 = u_(i + 1, j) - 2 * u_(i, j) + u_(i - 1, j);
    return du2 / (dx() * dx());
}

double Discretization::computeD2uDy2(int i, int j) const {
    const double du2 = u_(i, j + 1) - 2 * u_(i, j) + u_(i, j - 1);
    return du2 / (dy() * dy());
}

double Discretization::computeD2vDx2(int i, int j) const {
    const double dv2 = v_(i + 1, j) - 2 * v_(i, j) + v_(i - 1, j);
    return dv2 / (dx() * dx());
}

double Discretization::computeD2vDy2(int i, int j) const {
    const double dv2 = v_(i, j + 1) - 2 * v_(i, j) + v_(i, j - 1);
    return dv2 / (dy() * dy());
}

double Discretization::computeDuDx(int i, int j) const {
    const double du = u_(i, j) - u_(i - 1, j);
    return du / dx();
}

double Discretization::computeDvDy(int i, int j) const {
    const double dv = v_(i, j) - v_(i, j - 1);
    return dv / dy();
}

double Discretization::computeDpDx(int i, int j) const {
    const double dp = p_(i + 1, j) - p_(i, j);
    return dp / dx();
}

double Discretization::computeDpDy(int i, int j) const {
    const double dp = p_(i, j + 1) - p_(i, j);
    return dp / dy();
}
