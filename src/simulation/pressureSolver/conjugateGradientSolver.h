#pragma once

class ConjugateGradientSolver {

public:
    ConjugateGradient(const std::shared_ptr<StaggeredGrid> &grid, double epsilon, int maximumNumberOfIterations, double omega,
             const std::shared_ptr<Partitioning> &partitioning);

    ~ConjugateGradient() = default;

    void solve();

};
