#pragma once

#include "simulation/pressure/solver.h"
#include "grid/discretization.h"
#include <memory>

class sor final : public solver {
public:
    sor(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations, double omega);

    void solve() override;

private:
    std::shared_ptr<discretization> discretization_;
    double omega_;
};
