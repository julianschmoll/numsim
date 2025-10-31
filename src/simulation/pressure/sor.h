#pragma once

#include "simulation/pressure/solver.h"
#include "grid/discretization.h"
#include <memory>

class sor final : public solver {
public:
    sor(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations, double omega);

    /**
     * Solves the system of the Poisson equation for pressure.
     *
     * Overrides virtual solver method with successive over-relaxation solver.
     */
    void solve() override;

private:
    double omega_;
};
