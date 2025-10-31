#include "simulation/pressure/gaussSeidel.h"

#include <iostream>

void gaussSeidel::solve() {
#ifndef NDEBUG
    std::cout << "Solving pressure using Gauss Seidel" << std::endl;
#endif
    // TODO: implement Gauss-Seidel pressure solver
    int it = 0;
    double residual = 100.0;

    while (it<maxNumberOfIterations_ && residual > epsilon_) {
        setBoundaryValues();
        it++;
    }
}
