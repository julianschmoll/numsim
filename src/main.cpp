// #include "outputWriter.h"
#include "settings.h"
#include "solver.h"

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

	return EXIT_SUCCESS;
}
