#pragma once

#include "grid/discretization.h"

class solver {
    virtual void solve() = 0;

    void setBoundaryConditions();

public:
    virtual ~solver() = default;
    solver(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations);

protected:
    std::shared_ptr<discretization> discretization_;
    double epsilon_;
    double maxNumberOfIterations_;
};
