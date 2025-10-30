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

	outputWriterParaview_ = std::make_unique<outputWriterParaview>(outputWriterParaview(discretization_));
	outputWriterText_ = std::make_unique<outputWriterText>(outputWriterText(discretization_));
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

		return 0;
	}
	return 0;
}

void Simulation::setBoundaryValues() {
#ifndef NDEBUG
	std::cout << "Setting boundary values" << std::endl;
#endif
	const auto uIBegin = discretization_->uIBegin();
	const auto uIEnd = discretization_->uIEnd();
	const auto uJBegin = discretization_->uJBegin();
	const auto uJEnd = discretization_->uJEnd();

	const auto vIBegin = discretization_->vIBegin();
	const auto vIEnd = discretization_->vIEnd();
	const auto vJBegin = discretization_->vJBegin();
	const auto vJEnd = discretization_->vJEnd();

	const auto uBottom = settings_.dirichletBcBottom[0];
	const auto uTop = settings_.dirichletBcTop[0];
	const auto uLeft = settings_.dirichletBcLeft[0];
	const auto uRight = settings_.dirichletBcRight[0];

	const auto vBottom = settings_.dirichletBcBottom[1];
	const auto vTop = settings_.dirichletBcTop[1];
	const auto vLeft = settings_.dirichletBcLeft[1];
	const auto vRight = settings_.dirichletBcRight[1];

	// storing references here is probably faster than
	// always accessing through vtable
	auto& u = discretization_->u();
	auto& v = discretization_->v();
	auto& f = discretization_->f();
	auto& g = discretization_->g();

	for (int i = uIBegin; i < uIEnd; ++i) {
		u(i, uJBegin) = 2. * uBottom - u(i, uJBegin + 1);
		u(i, uJEnd - 1) = 2. * uTop - u(i, uJEnd - 2);
		f(i, uJBegin) = u(i, uJBegin);
		f(i, uJEnd - 1) = u(i, uJEnd - 1);
	}

	for (int j = uJBegin; j < uJEnd; ++j) {
		u(uIBegin, j) = uLeft;
		u(uIEnd - 1, j) = uRight;
		f(uIBegin, j) = u(uIBegin, j);
		f(uIEnd - 1, j) = u(uIEnd - 1, j);
	}

	for (int i = vIBegin; i < vIEnd; ++i) {
		v(i, vJBegin) = vBottom;
		v(i, vJEnd - 1) = vTop;
		g(i, vJBegin) = v(i, vJBegin);
		g(i, vJEnd - 1) = v(i, vJEnd - 1);
	}

	for (int j = vJBegin; j < vJEnd; ++j) {
		v(vIBegin, j) = 2. * vLeft - v(vIBegin + 1, j);
		v(vIEnd - 1, j) = 2. * vRight - v(vIEnd - 2, j);
		g(vIBegin, j) = v(vIBegin, j);
		g(vIEnd - 1, j) = v(vIEnd - 1, j);
	}
}

void Simulation::computeTimeStepWidth() {
#ifndef NDEBUG
	std::cout << "Computing time step width" << std::endl;
#endif
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

void Simulation::computeVelocities() {
#ifndef NDEBUG
	std::cout << "Computing velocities" << std::endl;
#endif
}