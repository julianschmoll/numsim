#include "grid/staggeredGrid.h"

StaggeredGrid::StaggeredGrid(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth) // TODO: Je nachdem, ob an globalem Rand oder nicht
    : meshWidth_(meshWidth), nCells_(nCells), u_({nCells[0] + 1, nCells[1] + 2}, meshWidth, {0.0, 0.5}, {-1, nCells_[0]}, {-1, nCells_[1] + 1}),
      v_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.0}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}),
      p_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}),
      f_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.0, 0.5}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}),
      g_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.0}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}),
      rhs_({nCells[0] + 2, nCells[1] + 2}, meshWidth, {0.5, 0.5}, {-1, nCells_[0] + 1}, {-1, nCells_[1] + 1}) {}

const std::array<double, 2> &StaggeredGrid::meshWidth() const {
    return meshWidth_;
}

const std::array<int, 2> &StaggeredGrid::nCells() const {
    return nCells_;
}

double StaggeredGrid::u(const int i, const int j) {
    return u_(i, j);
}

double StaggeredGrid::v(const int i, const int j) {
    return v_(i, j);
}

double StaggeredGrid::p(const int i, const int j) {
    return p_(i, j);
}

double StaggeredGrid::f(const int i, const int j) {
    return f_(i, j);
}

double StaggeredGrid::g(const int i, const int j) {
    return g_(i, j);
}

double StaggeredGrid::rhs(const int i, const int j) {
    return rhs_(i, j);
}

DataField &StaggeredGrid::u() {
    return u_;
}

DataField &StaggeredGrid::v() {
    return v_;
}

DataField &StaggeredGrid::p() {
    return p_;
}

DataField &StaggeredGrid::f() {
    return f_;
}

DataField &StaggeredGrid::g() {
    return g_;
}

DataField &StaggeredGrid::rhs() {
    return rhs_;
}

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}
