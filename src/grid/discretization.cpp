#include "grid/discretization.h"

discretization::discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth)
    : meshWidth_(meshWidth),
    nCells_(nCells),
    u_(nCells, meshWidth, {1.0, 0.5}),
    v_(nCells, meshWidth, {0.5, 1.0}),
    p_(nCells, meshWidth, {0.5, 0.5}),
    f_(nCells, meshWidth, {1.0, 0.5}),
    g_(nCells, meshWidth, {0.5, 1.0}),
    rhs_(nCells, meshWidth, {0.5, 0.5}) {

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
