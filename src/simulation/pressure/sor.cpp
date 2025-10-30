#include "simulation/pressure/sor.h"

#include <iostream>

sor::sor(std::shared_ptr<discretization> discretization, double epsilon, double maxNumberOfIterations, double omega)
    : solver(discretization, epsilon, maxNumberOfIterations), omega_(omega) {}

void sor::solve() {
#ifndef NDEBUG
    std::cout << "Solving pressure using SOR" << std::endl;
#endif
    // TODO: implement SOR pressure solver
}
