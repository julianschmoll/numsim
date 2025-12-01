#include "grid/staggeredGrid.h"
#include "grid/dataField.h"
#include "simulation/partitioning.h"
#include <array>

//TODO: Rename meshWidth -> meshResolution
StaggeredGrid::StaggeredGrid(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, const Partitioning &partitioning)
    : meshWidth_(meshWidth), nCells_(nCells)
{   
    int vWidth = nCells[0] + 2;
    int vHeight = nCells[1] + 1;

    int uWidth = nCells[0] + 1;
    int uHeight = nCells[1] + 2;

    const int pWidth = nCells[0] + 2;
    const int pHeight = nCells[1] + 2;

    // Add additional row/col as a new boundary to the top/right, otherwise the owner of the cells on the partition boundary cannot compute new values.
    bool topBorderPartition = partitioning.ownPartitionContainsBoundary(Direction::Top);
    bool rightBorderPartition = partitioning.ownPartitionContainsBoundary(Direction::Right);
    if (!topBorderPartition) {
        vHeight += 1;
    }
    if (!rightBorderPartition) {
        uWidth += 1;
    }

    p_   = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5});
    rhs_ = DataField({pWidth, pHeight}, meshWidth, {0.5, 0.5});

    u_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5});
    f_ = DataField({uWidth, uHeight}, meshWidth, {0.0, 0.5});

    v_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0});
    g_ = DataField({vWidth, vHeight}, meshWidth, {0.5, 0.0});

    int ownRankNo = partitioning.ownRankNo();
    auto rankCoords = partitioning.getCurrentRankCoords();

    std::cout << "rank " << ownRankNo << ": (" << rankCoords[0] << ", " << rankCoords[1] << ")\n";
    std::cout << " -- p_: size = [" << pWidth << ", " << pHeight << "]\n";
    std::cout << " -- v_: size = [" << vWidth << ", " << vHeight << "]\n";
    std::cout << " -- u_: size = [" << uWidth << ", " << uHeight << "]\n";
    std::cout << std::endl;
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

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}
