#include "grid/staggeredGrid.h"
#include "grid/array2d.h"
#include "grid/dataField.h"
#include "macros.h"
#include "settings.h"
#include "simulation/partitioning.h"
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

StaggeredGrid::StaggeredGrid(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, const Partitioning &partitioning)
    : meshWidth_(meshWidth), nCells_(nCells), partitioning_(partitioning) {
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
        displacementsTop_.resize(nCells[0], 0);
        topBoundaryPosition_.resize(nCells[0], 0); // TODO: initialize, settings.physicalSize[1]
    }
    if (partitioning.ownContainsBoundary<Direction::Bottom>()) {
        displacementsBottom_.resize(nCells[0], 0);
        bottomBoundaryPosition_.resize(nCells[0], 0);
    }

    fTop_.resize(pWidth, 0);
    fBottom_.resize(pWidth, 0);

    initializeStructureField();

    p_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, P_ID);
    q_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, Q_ID);
    rhs_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5}, RHS_ID);

    u_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5}, U_ID);
    f_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5}, F_ID);

    v_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0}, V_ID);
    g_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0}, G_ID);
}

double StaggeredGrid::globalDomainPosJ(int j) {
    return (partitioning_.nodeOffset()[1] + j) * dy(); // TODO: drüber nachdenken.
}

void StaggeredGrid::updateStructureCells()  {
    // top
    const int beginI = -1;
    const int beginJ = -1;
    const int endI = structure_.size()[0] - 2;
    const int endJ = structure_.size()[1] - 2;

    for (int j = beginJ; j <= endJ; ++j) { // inneres Feld

        for (int i = beginI + 1; i <= endI - 1; ++i) {
            assert(topBoundaryPosition_[i] > bottomBoundaryPosition_[i]);
            
            // double previousLowerCellEdge = globalDomainPosJ(j - 1);
            // double nextUpperCellEdge = globalDomainPosJ(j + 2);
            double lowerCellEdge = globalDomainPosJ(j);
            double upperCellEdge = globalDomainPosJ(j + 1);

            if (lowerCellEdge < bottomBoundaryPosition_[i] || upperCellEdge > topBoundaryPosition_[i]) {
                structure_(i, j) = Solid;
            } else {
                const bool isLowerDomainBoundary = (j == beginJ && partitioning_.ownContainsBoundary<Direction::Bottom>());
                const bool isUpperDomainBoundary = (j == endJ && partitioning_.ownContainsBoundary<Direction::Top>());
                // shouldn't be necessary, but there might be slight numerical deviances when we set the bottom & top boundaries
                if (isLowerDomainBoundary || isUpperDomainBoundary) {
                    continue;
                }
                if (isSolid(i, j) && (j == beginJ || isSolid(i, j - 1) /* previousLowerCellEdge <= bottomBoundaryPosition_[i] */)) { // Bottom: Solid -> Fluid
                    v_(i, j) = displacementsBottom_[i];
                    u_(i, j) = 0;
                    p_(i, j) = p_(i, j + 1); // TODO: notwendig?
                } else if (isSolid(i, j) && (j == endJ || isSolid(i, j + 1) /* nextUpperCellEdge >= topBoundaryPosition_[i] */)) { // Top: Solid -> Fluid
                    if (j < endJ) { // j == endJ is handled at the bottom of the partition above if it exists
                        v_(i, j) = displacementsBottom_[i];
                    }
                    u_(i, j) = 0;
                    p_(i, j) = p_(i, j + 1); // TODO: notwendig?
                }
                structure_(i, j) = Fluid;
            }
        }
    }
}

void StaggeredGrid::initializeStructureField() {
    // TODO: only initialize domain boundary with solid?!
    const int endI = structure_.size()[0] - 2;
    const int endJ = structure_.size()[1] - 2;

    for (int j = -1; j <= endJ; ++j) { // inneres Feld
        for (int i = -1; i <= endI - 1; ++i) {
            structure_(i, j) = Fluid;
        }
    }

    if (partitioning_.ownContainsBoundary<Direction::Top>()) {
        for (int i = -1; i <= endJ; i++) {
            structure_(i, endJ) = Solid;
        }
    }
    if (partitioning_.ownContainsBoundary<Direction::Bottom>()) {
        for (int i = -1; i <= endJ; i++) {
            structure_(i, -1) = Solid;
        }
    }
    if (partitioning_.ownContainsBoundary<Direction::Left>()) {
        for (int j = -1; j <= endI; j++) {
            structure_(-1, j) = Solid;
        } 
    }
    if (partitioning_.ownContainsBoundary<Direction::Right>()) {
        for (int j = -1; j <= endI; j++) {
            structure_(endI, j) = Solid;
        }
    }
}

void StaggeredGrid::test(const Settings &settings) {
    // TODO: Testcode: lid_driven_cavity.txt

    constexpr double eps = 1e-10; // TODO: dieses epsilon in denn Strukturcode verschieben?

    initializeStructureField();
    for (size_t i = 0; i < bottomBoundaryPosition_.size(); i++) {
        bottomBoundaryPosition_[i] = 0.1 * i - eps;
        topBoundaryPosition_[i] = settings.physicalSize[1];
    }
    updateStructureCells();
    std::cout << structure_ << "\n\n";

    initializeStructureField();
    for (size_t i = 0; i < bottomBoundaryPosition_.size(); i++) {
        bottomBoundaryPosition_[i] = 0;
        topBoundaryPosition_[i] = settings.physicalSize[1] - 0.1 * i + eps;
    }
    updateStructureCells();
    std::cout << structure_ << "\n\n";

    initializeStructureField();
    for (size_t i = 0; i < bottomBoundaryPosition_.size(); i++) {
        bottomBoundaryPosition_[i] = 0.2 + 0.2 * std::sin(i * 0.3);
        topBoundaryPosition_[i] = settings.physicalSize[1] - 0.2 + 0.2 * std::sin(i * 0.3);
    }
    updateStructureCells();
    std::cout << structure_ << "\n\n";
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
