#pragma once

#include "simulation/pressure/solver.h"

class gaussSeidel final : public solver {
public:
    using solver::solver;

    /**
     * Solves the system of the Poisson equation for pressure.
     */
    void solve() override;
};
