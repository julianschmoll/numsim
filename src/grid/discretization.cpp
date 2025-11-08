#include "grid/discretization.h"

Discretization::Discretization(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth)
    : meshWidth_(meshWidth), nCells_(nCells), u_({nCells[0] + 1, nCells[1] + 2}, meshWidth, {0.0, 0.5}, {-1, nCells_[0]}, {-1, nCells_[1] + 1}),
      v_({nCells[0] + 2, nCells[1] + 1}, meshWidth, {0.5, 0.0}, {-1, nCells_[0] + 1}, {-1, nCells_[1]}),
      p_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}),
      f_({nCells[0] + 1, nCells[1] + 2}, meshWidth, {0.0, 0.5}, {-1, nCells_[0]}, {-1, nCells_[1] + 1}),
      g_({nCells[0] + 2, nCells[1] + 1}, meshWidth, {0.5, 0.0}, {-1, nCells_[0] + 1}, {-1, nCells_[1]}),
      rhs_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}) {}

const std::array<double, 2> &Discretization::meshWidth() const {
    return meshWidth_;
}

const std::array<int, 2> &Discretization::nCells() const {
    return nCells_;
}

double Discretization::u(int i, int j) {
    return u_(i, j);
}

double Discretization::v(int i, int j) {
    return v_(i, j);
}

double Discretization::p(int i, int j) {
    return p_(i, j);
}

double Discretization::f(int i, int j) {
    return f_(i, j);
}

double Discretization::g(int i, int j) {
    return g_(i, j);
}

double Discretization::rhs(int i, int j) {
    return rhs_(i, j);
}

DataField &Discretization::u() {
    return u_;
}

DataField &Discretization::v() {
    return v_;
}

DataField &Discretization::p() {
    return p_;
}

DataField &Discretization::f() {
    return f_;
}

DataField &Discretization::g() {
    return g_;
}

DataField &Discretization::rhs() {
    return rhs_;
}

double Discretization::dx() const {
    return meshWidth_[0];
}

double Discretization::dy() const {
    return meshWidth_[1];
}
