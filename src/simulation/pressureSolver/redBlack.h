#pragma once
#include "pressureSolver.h"
#include "simulation/partitioning.h"


class RedBlack final : PressureSolver {
public:
    RedBlack(std::shared_ptr<StaggeredGrid> grid,
        double epsilon,
        int maximumNumberOfIterations,
        double omega,
        std::shared_ptr<Partitioning> partitioning);

    void solve() override;

protected:
    std::shared_ptr<Partitioning> partitioning_;

};
