#include "simulation.h"

#include "velocity/donorCell.h"


Simulation::Simulation(const Settings& settings) : settings_(settings) {
	std::cout << "Initializing Simulation..." << std::endl;

	if (settings_.useDonorCell) {
		std::cout << "Using Donor Cell." << std::endl;
		discretization_ = std::make_unique<donorCell>(settings_.nCells, settings_.physicalSize, settings_.alpha);
	}
	else {
		std::cout << "Using Central Differences." << std::endl;
		discretization_ = std::make_unique<centralDifferences>(settings_.nCells, settings_.physicalSize);
	}

	if (settings_.pressureSolver == "SOR") {
		std::cout << "Using SOR solver." << std::endl;
		solver_ = std::make_unique<sor>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations,
		                                settings_.omega);
	}
	else if (settings_.pressureSolver == "GaussSeidel") {
		solver_ = std::make_unique<
			gaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
	}
	else { throw std::runtime_error("Unknown pressure solver."); }
}

int Simulation::run() {
	double currentTime = 0.0;

	// simulation loop
	// ToDo: what if we reach more than the end time?
	while (currentTime < settings_.endTime) {
		setBoundaryValues();

		computeTimeStepWidth();

		currentTime += timeStepWidth_;

		computePreliminaryVelocities();

		computeRightHandSide();

		computePressure();

		computeVelocities();

		outputWriterParaview_->writeFile(currentTime);
		outputWriterText_->writeFile(currentTime);
	}
	return 0;
}

void Simulation::computePreliminaryVelocities() {
#ifndef NDEBUG
	std::cout << "Computing preliminary velocities" << std::endl;
#endif
}

void Simulation::computePressure() {
#ifndef NDEBUG
	std::cout << "Computing pressure" << std::endl;
#endif
}

void Simulation::computeRightHandSide() {
#ifndef NDEBUG
	std::cout << "Computing right hand side" << std::endl;
#endif
}

void Simulation::computeTimeStepWidth() {
#ifndef NDEBUG
	std::cout << "Computing time step width" << std::endl;
#endif
}

void Simulation::computeVelocities() {
#ifndef NDEBUG
	std::cout << "Computing velocities" << std::endl;
#endif
}

void Simulation::setBoundaryValues() {
#ifndef NDEBUG
	std::cout << "Computing boundary values" << std::endl;
#endif
}
