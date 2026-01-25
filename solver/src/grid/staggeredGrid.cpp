#include "grid/staggeredGrid.h"
#include "grid/array2d.h"
#include "grid/dataField.h"
#include "macros.h"
#include "simulation/partitioning.h"
#include <array>
#include <cassert>
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

    verticalNodeOffset_ =  partitioning.nodeOffset()[1];

    // TODO: think about parallelism
    structure_ = Array2d<bool>({pWidth, pHeight});
    if (partitioning.ownContainsBoundary<Direction::Top>()) {
        displacementsTop_.resize(pWidth, 0);
        topBoundaryPosition_.resize(pWidth, 0); // TODO: initialize
    }
    if (partitioning.ownContainsBoundary<Direction::Bottom>()) {
        displacementsBottom_.resize(pWidth, 0);
        bottomBoundaryPosition_.resize(pWidth, 0); // TODO: initialize
    }

    fTop_.resize(pWidth, 0);
    fBottom_.resize(pWidth, 0);

    // TODO: only initialize domain boundary with solid?!
    for (int j = -1; j <= nCells[1]; j++) {
        structure_(-1, j) = Solid;
        structure_(nCells[0], j) = Solid;
    }
    for (int i = -1; i <= nCells[0]; i++) {
        structure_(i, -1) = Solid;
        structure_(i, nCells[1]) = Solid;
    }

    // TODO: Testcode:
    for (int j = -1; j <= nCells[1]; j++) { 
        for (int i = -1; i <= nCells[0]; i++) {
            if (j <= 0.2 * i) {
                structure_(i, j) = Solid;
            }
        }
    }
    std::cout << structure_ << std::endl;

    p_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, P_ID);
    q_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, Q_ID);
    rhs_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, RHS_ID);

    u_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5}, U_ID);
    f_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5}, F_ID);

    v_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0}, V_ID);
    g_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0}, G_ID);
}

double StaggeredGrid::globalDomainPosJ(int j) {
    return (verticalNodeOffset_ + j) * dy(); // TODO: drüber nachdenken.
}

void StaggeredGrid::applyDisplacementsToBoundary()  {
    // TODO: wir sollten diese iterationsmethode irgendwie ändern. Code duplikation, unschön und fehleranfällig...
    // top
    for (int i = v_.beginI(); i < v_.endI(); ++i) {
        for (int j = v_.endJ() + 1; j >= v_.beginJ(); --j) {
            if (globalDomainPosJ(j) + 1 >= topBoundaryPosition_[i + 1]) { // + hier richtig?
                structure_(i, j) = Solid;
            } else if (isSolid(i, j)) { // nur dann müssen wir tatsächlich arbeit investieren
                structure_(i, j) = Fluid;
                // update cell contents?
            }
        }
    }
    // bottom
    for (int i = v_.beginI(); i < v_.endI(); ++i) {
        for (int j = v_.beginJ(); j < v_.endJ() - 1; ++j) {
            // aufpassen, dass wir nicht den echten simulationsrand verändern!
            if (globalDomainPosJ(j) <= bottomBoundaryPosition_[i + 1]) {
                structure_(i, j) = Solid;
            } else if (isSolid(i, j)) {
                structure_(i, j) = Fluid;
                // update cell contents?
            }
        }
    }
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

DataField &StaggeredGrid::q() {
    return q_;
}

DataField &StaggeredGrid::rhs() {
    return rhs_;
}

double &StaggeredGrid::fTop(int i) {
    assert(0 <= i + 1 && i + 1 < static_cast<int>(fTop_.size()));
    return fTop_[i + 1];
}

double &StaggeredGrid::fBottom(int i) {
    assert(0 <= i + 1 && i + 1 < static_cast<int>(fBottom_.size()));
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
    return !structure_(i, j);
}

bool StaggeredGrid::isSolid(int i, int j) const {
    return structure_(i, j);
}
