#include "pressureSolver.h"

PressureSolver::PressureSolver(std::shared_ptr<StaggeredGrid> grid, double epsilon, double maxNumberOfIterations, double omega)
    : grid_(grid), epsilon_(epsilon), maxNumberOfIterations_(maxNumberOfIterations), omega_(omega) {}


