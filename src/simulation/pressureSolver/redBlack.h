#pragma once

#include "simulation/partitioning.h"
#include "grid/staggeredGrid.h"
#include <memory>


class RedBlack final {
public:
    RedBlack(const std::shared_ptr<StaggeredGrid> &grid,
        double epsilon,
        int maximumNumberOfIterations,
        double omega,
        const std::shared_ptr<Partitioning> &partitioning);

    ~RedBlack();

    void solve();

private:
    void updatePressureBoundaries();

    void setBoundaryValues();

protected:
    std::shared_ptr<Partitioning> partitioning_;

    // object holding the needed field variables for rhs and p
    std::shared_ptr<StaggeredGrid> grid_;

    // Convergence threshold for the residual
    double epsilon_;

    // Maximum number of iterations allowed in the solver
    double maxNumberOfIterations_;

    // Relaxation factor
    double omega_;

    // MPI Requests buffer
    std::vector<MPI_Request> requests_;

    double localPressureError_ = 0;

    double globalPressureError_ = 0;

    MPI_Datatype mpiColumn_;
};
