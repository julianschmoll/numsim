#include "solver.h"


int Solver::AdvanceTimeStep() {
    // this is just some boilerplate code, here we would actually call the solve method
    time += settings_.maximumDt;
    std::cout << "current time: " << time << " seconds" << std::endl;

    // compute_fg
    // for(int j = 0; j < grid.ySize(); j++) {
        //for(int i = 0; i < grid.xSize(); i++) {
            // CellEnvironment cell = grid.get(i, j);
            // F_ij := ...
            // G_ij := ...
        //}
    //}

    // sor
    // do {
        // RHS_ij := grid.getRHS(i, j) // compute on the fly?
        // p_ij := sor step 
        // res := res + (RHS_ij - discretized poisson equation LHS_ij)^2
    // } while (res > eps * eps  && iter < maxIterations)



    return 0;
}
