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

	// simulation loop, do while so we definitely reach end time
	do {
		setBoundaryValues();

		currentTime += computeTimeStepWidth();

		computePreliminaryVelocities();

		computeRightHandSide();

		computePressure();

		computeVelocities();

		outputWriterParaview_->writeFile(currentTime);
		outputWriterText_->writeFile(currentTime);

	} while (currentTime < settings_.endTime);
	return 0;
}

// ToDo: Do we need to set F and G?
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

	for (int i = uIBegin; i < uIEnd; ++i) {
		u(i, uJBegin) =  uBottom;
		u(i, uJEnd - 1) = uTop;
	}

	for (int j = uJBegin; j < uJEnd; ++j) {
		u(uIBegin, j) = uLeft;
		u(uIEnd - 1, j) = uRight;
	}

	for (int i = vIBegin; i < vIEnd; ++i) {
		v(i, vJBegin) = vBottom;
		v(i, vJEnd - 1) = vTop;
	}

	for (int j = vJBegin; j < vJEnd; ++j) {
		v(vIBegin, j) = vLeft;
		v(vIEnd - 1, j) = vRight;
	}
}

double Simulation::computeTimeStepWidth() const {
#ifndef NDEBUG
	std::cout << "Computing time step width" << std::endl;
#endif
	const double uMax  = discretization_->u().max();
	const double vMax = discretization_->v().max();

	const double hx = settings_.physicalSize[0] / settings_.nCells[0];
	const double hy = settings_.physicalSize[1] / settings_.nCells[1];

	const double diffusive = (((hx * hx) + (hy * hy)) / 2.0) * (settings_.re / 2.0);

	// we need to handle velocities being 0
	const double convectiveU = (uMax > 0) ? hx / uMax : std::numeric_limits<double>::max();
	const double convectiveV = (vMax > 0) ? hy / vMax : std::numeric_limits<double>::max();

	const double dt = std::min({diffusive, convectiveU, convectiveV}) * settings_.tau;

#ifndef NDEBUG
	std::cout << "Computed time step width: " << dt << std::endl;
#endif

	return std::min(dt, settings_.maximumDt);
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