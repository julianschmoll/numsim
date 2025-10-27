// #include "outputWriter.h"
#include "settings.h"
#include "solver/solver.h"

#include <iostream>

int main(int argc, char *argv[]) {
    // we need an input file being specified
	if (argc == 1) {
		std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
 		return EXIT_FAILURE;
 	}
	std::string filename = argv[1];

	Settings settings;
 	settings.loadFromFile(filename);

#ifndef NDEBUG
	settings.printSettings();
#endif

	Solver solver(settings);

	// outputWriter needs solver discretization
	// outputWriterParaview outputWriter(solver.discretization());

	double outputTime = solver.time;
	double step = 0.1;

	// Simulation loop
	// might have to be do while so we get last output
	while (solver.time <= settings.endTime){
		solver.AdvanceTimeStep();

		if (solver.time >= outputTime) {
			// outputWriter.writeFile(solver.time);
			// outputTime += step;
		}
	}


	//cell.v.derivative();
	//cell.u.derivative();

	//cell.p.xDerivative();
	//cell.p.yDerivative();

	return EXIT_SUCCESS;
}
