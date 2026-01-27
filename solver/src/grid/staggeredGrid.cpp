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
#include <vector>

StaggeredGrid::StaggeredGrid(const Settings &settings, const Partitioning &partitioning)
    : nCells_(settings.nCells),
      meshWidth_({settings.physicalSize[0] / nCells_[0], settings.physicalSize[1] / nCells_[1]}), 
      partitioning_(partitioning) {

    int vWidth = nCells_[0] + 2;
    int vHeight = nCells_[1] + 1;

    int uWidth = nCells_[0] + 1;
    int uHeight = nCells_[1] + 2;

    const int pWidth = nCells_[0] + 2;
    const int pHeight = nCells_[1] + 2;

    // Add additional row/col as a new boundary to the top/right, otherwise the owner of the cells on the partition boundary cannot compute new
    // values.
    if (!partitioning.ownContainsBoundary<Direction::Top>()) {
        vHeight += 1;
    }
    if (!partitioning.ownContainsBoundary<Direction::Right>()) {
        uWidth += 1;
    }

    structure_ = Array2d<bool>({pWidth, pHeight});

    displacementsTop_ = std::vector<double>(nCells_[0], 0);
    topBoundaryPosition_ = std::vector<double>(nCells_[0], settings.physicalSize[1]); // TODO: initialize, settings.physicalSize[1]

    displacementsBottom_ = std::vector<double>(nCells_[0], 0);
    bottomBoundaryPosition_ = std::vector<double>(nCells_[0], 0);

    fTop_ = std::vector<double>(pWidth, 0);
    fBottom_ = std::vector<double>(pWidth, 0);

    initializeStructureField();

    p_ = DataField({pWidth, pHeight}, meshWidth_, {0.5, 0.5}, P_ID);
    q_ = DataField({pWidth, pHeight}, meshWidth_, {0.5, 0.5}, Q_ID);
    rhs_ = DataField({pWidth, pHeight}, meshWidth_, {0.5, 0.5}, RHS_ID);

    u_ = DataField({uWidth, uHeight}, meshWidth_, {0.0, 0.5}, U_ID);
    f_ = DataField({uWidth, uHeight}, meshWidth_, {0.0, 0.5}, F_ID);

    v_ = DataField({vWidth, vHeight}, meshWidth_, {0.5, 0.0}, V_ID);
    g_ = DataField({vWidth, vHeight}, meshWidth_, {0.5, 0.0}, G_ID);
}

double StaggeredGrid::globalDomainPosJ(int j) {
    return (partitioning_.nodeOffset()[1] + j) * dy(); // TODO: drüber nachdenken.
}

void StaggeredGrid::updateStructureCells()  {
    const int beginI = structure_.minI();
    const int beginJ = structure_.minJ();
    const int endI = structure_.maxI();
    const int endJ = structure_.maxJ();

    // We use a small offset to make sure a cell is actually fully solid.
    // This way we also avoid labeling the domain boundary as fluid.
    constexpr double eps = 1e-10;

    // We iterate over all vertical coordinates since there is no communication of the structure array.
    for (int j = beginJ; j <= endJ; ++j) {

        for (int i = beginI + 1; i <= endI - 1; ++i) {
            assert(topBoundaryPosition_[i] > bottomBoundaryPosition_[i]);
            
            double lowerCellEdge = globalDomainPosJ(j);
            double upperCellEdge = globalDomainPosJ(j + 1);

            if (lowerCellEdge < bottomBoundaryPosition_[i] - eps || upperCellEdge > topBoundaryPosition_[i] + eps) {
                structure_(i, j) = Solid;
            } else {
                if (isSolid(i, j) && (j == beginJ || isSolid(i, j - 1))) { // Bottom: Solid -> Fluid
                    v_(i, j) = displacementsBottom_[i];
                    u_(i, j) = 0;
                    p_(i, j) = p_(i, j + 1);
                } else if (isSolid(i, j) && (j == endJ || isSolid(i, j + 1))) { // Top: Solid -> Fluid
                    if (j < endJ) { 
                        // j == endJ is handled at the bottom of the partition above if it exists.
                        v_(i, j) = displacementsTop_[i];
                    }
                    u_(i, j) = 0;
                    p_(i, j) = p_(i, j - 1);
                }
                structure_(i, j) = Fluid;
            }
        }
    }
}

void StaggeredGrid::initializeStructureField() {
    // TODO: only initialize domain boundary with solid?!
    const int endI = structure_.maxI();
    const int endJ = structure_.maxJ();

    for (int j = -1; j <= endJ; ++j) { // inneres Feld
        for (int i = -1; i <= endI - 1; ++i) {
            structure_(i, j) = Fluid;
        }
    }

    if (partitioning_.ownContainsBoundary<Direction::Top>()) {
        for (int i = -1; i <= endI; i++) {
            structure_(i, endJ) = Solid;
        }
    }
    if (partitioning_.ownContainsBoundary<Direction::Bottom>()) {
        for (int i = -1; i <= endI; i++) {
            structure_(i, -1) = Solid;
        }
    }
    /*
    if (partitioning_.ownContainsBoundary<Direction::Left>()) {
        for (int j = -1; j <= endJ; j++) {
            structure_(-1, j) = Solid;
        } 
    }
    if (partitioning_.ownContainsBoundary<Direction::Right>()) {
        for (int j = -1; j <= endJ; j++) {
            structure_(endI, j) = Solid;
        }
    }
    */
}

void StaggeredGrid::test(const Settings &settings) {
    // TODO: Testcode: lid_driven_cavity.txt

    constexpr double eps = 1e-10; // TODO: dieses epsilon in denn Strukturcode verschieben?

    const int iOffset = partitioning_.nodeOffset()[0];

    initializeStructureField();
    for (size_t i = 0; i < bottomBoundaryPosition_.size(); i++) {
        bottomBoundaryPosition_[i] = 0.1 * (i + iOffset);
        topBoundaryPosition_[i] = settings.physicalSize[1];
    }
    updateStructureCells();
    std::cout << structure_ << "\n\n";

    initializeStructureField();
    for (size_t i = 0; i < bottomBoundaryPosition_.size(); i++) {
        bottomBoundaryPosition_[i] = 0;
        topBoundaryPosition_[i] = settings.physicalSize[1] - 0.1 * (i + iOffset);
    }
    updateStructureCells();
    std::cout << structure_ << "\n\n";

    if (partitioning_.nRanks() == 1) {
        initializeStructureField();
        for (size_t i = 0; i < bottomBoundaryPosition_.size(); i++) {
            bottomBoundaryPosition_[i] = 0.2 + 0.2 * std::sin((i + iOffset) * 0.3);
            topBoundaryPosition_[i] = settings.physicalSize[1] - 0.2 + 0.2 * std::sin((i + iOffset) * 0.3);
        }
        updateStructureCells();
        std::cout << structure_ << "\n\n";
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
