// #include "outputWriter.h"
#include "settings.h"
#include "simulation/simulation.h"

#include <iostream>

int main(const int argc, char* argv[]) {
	// we need an input file being specified
	if (argc == 1) {
		std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
		return EXIT_FAILURE;
	}
	const std::string filename = argv[1];

	Settings settings;
	settings.loadFromFile(filename);

#ifndef NDEBUG
	settings.printSettings();
#endif

	Simulation simulation(settings);

	// outputWriter needs solver discretization
	// outputWriterParaview outputWriter(solver.discretization());

	simulation.run();


	//cell.v.derivative();
	//cell.u.derivative();

	//cell.p.xDerivative();
	//cell.p.yDerivative();

	return EXIT_SUCCESS;
}
