#include "simulation/velocity/donorCell.h"

donorCell::donorCell(const std::array<int, 2> &nCells, const std::array<double, 2> &meshWidth, double alpha)
    : discretization(nCells, meshWidth), alpha_(alpha)
{ }

double donorCell::computeDu2Dx(int i, int j) const {
    return 0;
}

double donorCell::computeDv2Dy(int i, int j) const {
    return 0;
}

double donorCell::computeDuvDx(int i, int j) const {
    return 0;
}

double donorCell::computeDuvDy(int i, int j) const {
    return 0;
}