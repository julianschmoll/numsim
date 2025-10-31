#include "simulation.h"
#include "velocity/donorCell.h"
#include "velocity/centralDifferences.h"

#include "outputWriter/outputWriterParaview.h"
#include "outputWriter/outputWriterText.h"

#include <memory>


#include "simulation/pressure/gaussSeidel.h"
#include "simulation/pressure/sor.h"

Simulation::Simulation(Settings settings) : settings_(settings) {
	std::cout << "Initializing Simulation..." << std::endl;

	if (settings_.useDonorCell) {
		std::cout << "Using Donor Cell." << std::endl;
		//discretization_ = std::make_shared<//>(settings_.nCells, settings_.physicalSize, settings_.alpha);
	}
	else {
		std::cout << "Using Central Differences." << std::endl;
		//discretization_ = std::make_shared<centralDifferences>(settings_.nCells, settings_.physicalSize);
	}

	if (settings_.pressureSolver == "SOR") {
		std::cout << "Using SOR solver." << std::endl;
		//solver_ = std::make_shared<sor>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations,
		//                                settings_.omega);
	}
	else if (settings_.pressureSolver == "GaussSeidel") {
		// solver_ = std::make_shared<
			// gaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
	} else {
		throw std::runtime_error("Unknown pressure solver.");
	}
	std::array<int, 2> a = {10, 10};
	discretization_ = std::make_shared<centralDifferences>(a, settings_.physicalSize);
}

void Simulation::run() {

	for (int i = discretization_->uIBegin(); i < std::min(discretization_->uIEnd(), discretization_->uJEnd()); i++) {
		discretization_->u()(i, i) = 1;
	}

	for (int i = discretization_->vIBegin(); i < std::min(discretization_->vIEnd(), discretization_->vJEnd()); i++) {
		discretization_->v()(i, i) = 1;
	}

	for (int i = discretization_->pIBegin(); i < std::min(discretization_->pIEnd(), discretization_->pJEnd()); i++) {
		discretization_->p()(i, i) = 1;
	}

	outputWriterParaview outputWriterParaview(discretization_);
	outputWriterText outputWriterText(discretization_);

	outputWriterParaview.writeFile(1);
	outputWriterText.writeFile(1);
}
