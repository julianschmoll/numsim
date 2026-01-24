#include "grid/staggeredGrid.h"
#include "grid/array2d.h"
#include "grid/dataField.h"
#include "macros.h"
#include "simulation/partitioning.h"
#include <array>
#include <iostream>

StaggeredGrid::StaggeredGrid(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, const Partitioning &partitioning)
    : meshWidth_(meshWidth), nCells_(nCells) {
    int vWidth = nCells[0] + 2;
    int vHeight = nCells[1] + 1;

    int uWidth = nCells[0] + 1;
    int uHeight = nCells[1] + 2;

    const int pWidth = nCells[0] + 2;
    const int pHeight = nCells[1] + 2;

    // Add additional row/col as a new boundary to the top/right, otherwise the owner of the cells on the partition boundary cannot compute new
    // values.
    if (!partitioning.ownContainsBoundary<Direction::Top>()) {
        vHeight += 1;
    }
    if (!partitioning.ownContainsBoundary<Direction::Right>()) {
        uWidth += 1;
    }

    structure_ = Array2d<CellType>({pWidth, pHeight});
    if (partitioning.ownContainsBoundary<Direction::Top>()) {
        newDisplacementsTop_.resize(pWidth, 0);
        oldDisplacementsTop_.resize(pWidth, 0);
    }
    if (partitioning.ownContainsBoundary<Direction::Bottom>()) {
        newDisplacementsBottom_.resize(pWidth, 0);
        oldDisplacementsBottom_.resize(pWidth, 0);
    }

    fTop_.resize(pWidth, 0);
    fBottom_.resize(pWidth, 0);

    // We initialize this field with a border of solid.
    // This should also not be changed by the displacements since it would break the velocity boundarie conditions.
    for (int j = -1; j <= nCells[1]; j++) { // TODO: move the begin and end methods to Array2d? They don't use DataField specific information.
        structure_(-1, j) = CellType::Solid;
        structure_(nCells[0], j) = CellType::Solid;
    }
    for (int i = -1; i <= nCells[0]; i++) {
        structure_(i, -1) = CellType::Solid;
        structure_(i, nCells[1]) = CellType::Solid;
    }

    // TODO: Testcode:
    for (int j = -1; j <= nCells[1]; j++) { 
        for (int i = -1; i <= nCells[0]; i++) {
            if (j <= 0.2 * i) {
                structure_(i, j) = CellType::Solid;
            }
        }
    }
    std::cout << structure_ << std::endl;

    p_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, P_ID);
    rhs_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, RHS_ID);

    u_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5}, U_ID);
    f_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5}, F_ID);

    v_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0}, V_ID);
    g_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0}, G_ID);
}

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

double &StaggeredGrid::fTop(int i) {
    return fTop_[i];
}

double &StaggeredGrid::fBottom(int i) {
    return fBottom_[i];
}

std::vector<double> &StaggeredGrid::fTop() {
    return fTop_;
}

std::vector<double> &StaggeredGrid::fBottom() {
    return fBottom_;
}

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}

bool StaggeredGrid::isFluid(int i, int j) const {
    return structure_(i, j) == CellType::Fluid;
}

bool StaggeredGrid::isSolid(int i, int j) const {
    return structure_(i, j) == CellType::Solid;
}
