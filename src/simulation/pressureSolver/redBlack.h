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

    ~RedBlack() = default;

    void solve();

private:
    void exchange();

    void setBoundaryValues(Direction direction);
    template <class Func>
    void copyBoundary(Func set, int begin, int end);
    void setVerticalBoundaryValues(int iIndex, int offset, DataField &df);
    void setHorizontalBoundaryValues(int jIndex, int offset, DataField &df);

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

    int lastIterationCount_ = 0;
};
