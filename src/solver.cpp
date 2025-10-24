#include "solver.h"


int Solver::AdvanceTimeStep() {
    // this is just some boilerplate code, here we would actually call the solve method
    time += settings_.maximumDt;
    std::cout << "current time: " << time << " seconds" << std::endl;
    return 0;
}
