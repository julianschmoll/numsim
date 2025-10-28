#pragma once

#include "simulation/pressure/solver.h"

class gaussSeidel final : public solver {
public:
    using solver::solver;

    void solve() override;
};
