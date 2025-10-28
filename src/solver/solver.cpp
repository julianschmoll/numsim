#include "solver.h"

#include <utility>

Solver::Solver(Settings  settings) : time(0.0), settings_(std::move(settings)) {
	std::cout << "Initializing Solver..." << std::endl;

	if (settings_.useDonorCell) {
		std::cout << "Using Donor Cell." << std::endl;
		// discretization_ = std::make_shared<DonorCell>(settings_.nCells, settings_.physicalSize, settings_.alpha);
	} else {
		std::cout << "Using Central Differences." << std::endl;
		// discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, settings_.physicalSize);
	}

	if (settings_.pressureSolver == "SOR") {
		std::cout << "Using SOR solver." << std::endl;
		// pressureSolver_ = std::make_shared<SOR>(discretization_, settings_.omega);
	} else {
		// ?
	}

	// initialize fields with zero


	// randwerte außerhalb des grids einmal berechnen
}


int Solver::AdvanceTimeStep() {
	// this is just some boilerplate code, here we would actually call the solve method

	// compute stable timestep
	// apply boundary conditions
	// compute velocities
	// compute pressure
	// correct velocities
	// apply boundary conditions

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
