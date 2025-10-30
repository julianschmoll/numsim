#include "simulation/pressure/solver.h"

solver::solver(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations)
    : discretization_(discretization), epsilon_(epsilon), maxNumberOfIterations_(maxNumberOfIterations) {}

void solver::setBoundaryValues() {
    // TODO: implement boundary condition logic
}
