#include "simulation/pressure/sor.h"

sor::sor(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations, double omega)
    : solver(discretization, epsilon, maxNumberOfIterations), omega_(omega) {}

void sor::solve() {
    // TODO: implement SOR pressure solver
}
