// #include "outputWriter.h"
#include "settings.h"
#include "solver.h"

#include "grid.h"

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
	//OutputWriterParaview outputWriter;

	// Simulation loop
	while (solver.time <= settings.endTime){
		solver.AdvanceTimeStep();
		//outputWriter.writeFile();
	}


	CellEnvironment cell;

	cell.u.u[0] = 1;
	cell.u.u[1] = 2;
	cell.u.u[2] = 3;

	cell.p.p[0] = 1;
	cell.p.p[1] = 2;
	cell.p.p[2] = 3;
	cell.p.p[3] = 4;
	cell.p.p[4] = 5;

	int i = 0, j = 0;

	//cell.v.derivative();
	//cell.u.derivative();

	//cell.p.xDerivative();
	//cell.p.yDerivative();

	std::cout << cell.u(-1) << " " << cell.u(0) << " " << cell.u(1) << std::endl;

	std::cout << "  "                 << cell.p(0, -1)       << std::endl;
	std::cout << cell.p(-1, 0) << " " << cell.p(0, 0) << " " << cell.p(1, 0) << std::endl;
	std::cout << "  "                 << cell.p(0, 1)        << std::endl;


	return EXIT_SUCCESS;
}
