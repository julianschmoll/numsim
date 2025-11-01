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
		pressureSolver_ = std::make_unique<sor>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations,
		                                        settings_.omega);
	}
	else if (settings_.pressureSolver == "GaussSeidel") {
		pressureSolver_ = std::make_unique<
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

		computeTimeStepWidth();
		currentTime += timeStepWidth_;

		setPreliminaryVelocities();

		setRightHandSide();

		// computes the pressure
		pressureSolver_->solve();

		setVelocities();

		outputWriterParaview_->writeFile(currentTime);
		outputWriterText_->writeFile(currentTime);
	}
	while (currentTime < settings_.endTime);
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
		u(i, uJBegin) = uBottom;
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

void Simulation::computeTimeStepWidth() {
#ifndef NDEBUG
	std::cout << "Computing time step width" << std::endl;
#endif
	const double uMax = discretization_->u().absMax();
	const double vMax = discretization_->v().absMax();

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

	timeStepWidth_ = std::min(dt, settings_.maximumDt);
}

void Simulation::setPreliminaryVelocities() {

	const double invRe = 1.0 / settings_.re;

	auto& u = discretization_->u();
	auto& v = discretization_->v();
	auto& f = discretization_->f();
	auto& g = discretization_->g();

	const auto& uJBegin = discretization_->uJBegin();
	const auto& uJEnd = discretization_->uJEnd();
	const auto& uIBegin = discretization_->uIBegin();
	const auto& uIEnd = discretization_->uIEnd();

	const auto& vJBegin = discretization_->vJBegin();
	const auto& vJEnd = discretization_->vJEnd();
	const auto& vIBegin = discretization_->vIBegin();
	const auto& vIEnd = discretization_->vIEnd();

	for (int j = uJBegin; j < uJEnd; j++) {
		for (int i = uIBegin; i < uIEnd; i++) {
			double d2udx2 = discretization_->computeD2uDx2(i, j);
			double du2dx = discretization_->computeDu2Dx(i, j);
			double d2udy2 = discretization_->computeD2uDy2(i, j);
			double duvdy = discretization_->computeDuvDy(i, j);

			double convection = du2dx + duvdy;
			double diffusion = d2udx2 + d2udy2;

			f(i, j) = u(i, j) + timeStepWidth_ * (invRe * diffusion - convection + 0); // TODO external forces
		}
	}

	for (int j = vJBegin; j < vJEnd; j++) {
		for (int i = vIBegin; i < vIEnd; i++) {
			double d2vdx2 = discretization_->computeD2vDx2(i, j);
			double dv2dy = discretization_->computeDv2Dy(i, j);
			double d2vdy2 = discretization_->computeD2vDy2(i, j);
			double duvdx = discretization_->computeDuvDx(i, j);

			double convection = dv2dy + duvdx;
			double diffusion = d2vdx2 + d2vdy2;

			g(i, j) = v(i, j) + timeStepWidth_ * (invRe * diffusion - convection + 0); // TODO external forces
		}
	}
}

void Simulation::setRightHandSide() {

	const auto& pJBegin = discretization_->pJBegin();
	const auto& pJEnd = discretization_->pJEnd();
	const auto& pIBegin = discretization_->pIBegin();
	const auto& pIEnd = discretization_->pIEnd();

	const auto& dx = discretization_->dx();
	const auto& dy = discretization_->dy();

	const auto& f = discretization_->f();
	const auto& g = discretization_->g();
	auto& rhs = discretization_->rhs();

	const double invTimeStep = 1.0 / timeStepWidth_;
	const double invDx = 1.0 / dx;
	const double invDy = 1.0 / dy;

	for (int j = pJBegin; j < pJEnd; j++) {
		for (int i = pIBegin; i < pIEnd; i++) {
			const double diffF = (f(i, j) - f(i - 1, j)) * invDx;
			const double diffG = (g(i, j) - g(i, j - 1)) * invDy;
			rhs(i, j) = invTimeStep * (diffF + diffG);
		}
	}
}

void Simulation::setVelocities() {

	const auto& uJBegin = discretization_->uJBegin();
	const auto& uJEnd = discretization_->uJEnd();
	const auto& uIBegin = discretization_->uIBegin();
	const auto& uIEnd = discretization_->uIEnd();

	const auto& vJBegin = discretization_->vJBegin();
	const auto& vJEnd = discretization_->vJEnd();
	const auto& vIBegin = discretization_->vIBegin();
	const auto& vIEnd = discretization_->vIEnd();

	auto& u = discretization_->u();
	auto& v = discretization_->v();
	auto& f = discretization_->f();
	auto& g = discretization_->g();

	for (int j = uJBegin; j < uJEnd; j++) {
		for (int i = uIBegin; i < uIEnd; i++) {
			u(i, j) = f(i, j) - timeStepWidth_ * discretization_->computeDpDx(i, j);
		}
	}

	for (int j = vJBegin; j < vJEnd; j++) {
		for (int i = vIBegin; i < vIEnd; i++) {
			v(i, j) = g(i, j) - timeStepWidth_ * discretization_->computeDpDy(i, j);
		}
	}
}
