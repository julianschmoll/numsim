#include "simulation.h"

#include "velocity/donorCell.h"


Simulation::Simulation(const Settings& settings) : settings_(settings) {
	std::cout << "Initializing Simulation..." << std::endl;

	for (int i = 0; i < 2; ++i) { meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i]; }

	if (settings_.useDonorCell) {
		std::cout << "Using Donor Cell." << std::endl;
		discretization_ = std::make_unique<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
	}
	else {
		std::cout << "Using Central Differences." << std::endl;
		discretization_ = std::make_unique<CentralDifferences>(settings_.nCells, meshWidth_);
	}

	if (settings_.pressureSolver == "SOR") {
		std::cout << "Using SOR solver." << std::endl;
		pressureSolver_ = std::make_unique<PressureSolver>(discretization_, settings_.epsilon,
		                                                   settings_.maximumNumberOfIterations, settings_.omega);
	}
	else if (settings_.pressureSolver == "GaussSeidel") {
		pressureSolver_ = std::make_unique<PressureSolver>(discretization_, settings_.epsilon,
		                                                   settings_.maximumNumberOfIterations, 1);
	}
	else { throw std::runtime_error("Unknown pressure solver."); }

	outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
	outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
}

int Simulation::run() {
	double currentTime = 0.0;

	while (currentTime < settings_.endTime) {
		setBoundaryValues();

		computeTimeStepWidth();

		if (currentTime + timeStepWidth_ > settings_.endTime) { timeStepWidth_ = settings_.endTime - currentTime; }
		currentTime += timeStepWidth_;

		setPreliminaryVelocities();

		setRightHandSide();

		// computes the pressure
		pressureSolver_->solve();

		setVelocities();

		outputWriterParaview_->writeFile(currentTime);
		outputWriterText_->writeFile(currentTime);
	}
	return 0;
}

void Simulation::setBoundaryValues() {
#ifndef NDEBUG
	std::cout << "Setting boundary values" << std::endl;
#endif
	const auto uBottom = settings_.dirichletBcBottom[0];
	const auto uTop = settings_.dirichletBcTop[0];
	const auto uLeft = settings_.dirichletBcLeft[0];
	const auto uRight = settings_.dirichletBcRight[0];

	const auto vBottom = settings_.dirichletBcBottom[1];
	const auto vTop = settings_.dirichletBcTop[1];
	const auto vLeft = settings_.dirichletBcLeft[1];
	const auto vRight = settings_.dirichletBcRight[1];

	auto& u = discretization_->u();
	auto& v = discretization_->v();
	auto& f = discretization_->f();
	auto& g = discretization_->g();

	for (int i = u.beginI(); i < u.endI(); ++i) {
		u(i, u.beginJ()) = 2.0 * uBottom - u(i, u.beginJ() + 1);
		u(i, u.endJ() - 1) = 2.0 * uTop - u(i, u.endJ() - 2);
	}
	for (int i = v.beginI(); i < v.endI(); ++i) {
		v(i, v.beginJ()) = vBottom;
		v(i, v.endJ() - 1) = vTop;
	}

	for (int j = u.beginJ(); j < u.endJ(); ++j) {
		u(u.beginI(), j) = uLeft;
		u(u.endI() - 1, j) = uRight;
	}
	for (int j = v.beginJ(); j < v.endJ(); ++j) {
		v(v.beginI(), j) = 2.0 * vLeft - v(v.beginI() + 1, j);
		v(v.endI() - 1, j) = 2.0 * vRight - v(v.endI() - 2, j);
	}

	for (int i = f.beginI(); i < f.endI(); ++i) {
		f(i, f.beginJ()) = u(i, u.beginJ());
		f(i, f.endJ() - 1) = u(i, u.endJ() - 1);
	}
	for (int i = g.beginI(); i < g.endI(); ++i) {
		g(i, g.beginJ()) = v(i, v.beginJ());
		g(i, g.endJ() - 1) = v(i, v.endJ() - 1);
	}

	for (int j = f.beginJ(); j < f.endJ(); ++j) {
		f(f.beginI(), j) = u(u.beginI(), j);
		f(f.endI() - 1, j) = u(u.endI() - 1, j);
	}
	for (int j = g.beginJ(); j < g.endJ(); ++j) {
		g(g.beginI(), j) = v(v.beginI(), j);
		g(g.endI() - 1, j) = v(v.endI() - 1, j);
	}
}

void Simulation::computeTimeStepWidth() {
#ifndef NDEBUG
	std::cout << "Computing time step width" << std::endl;
#endif
	const double uMax = discretization_->u().absMax();
	const double vMax = discretization_->v().absMax();

	const double hx = discretization_->dx();
	const double hy = discretization_->dy();

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

	for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
		for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
			double uDxx = discretization_->computeD2uDx2(i, j);
			double u2Dx = discretization_->computeDu2Dx(i, j);
			double uDyy = discretization_->computeD2uDy2(i, j);
			double uvDy = discretization_->computeDuvDy(i, j);

			double convection = u2Dx + uvDy;
			double diffusion = uDxx + uDyy;

			f(i, j) = u(i, j) + timeStepWidth_ * (invRe * diffusion - convection + 0); // TODO external forces
		}
	}

	for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
		for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
			double vDxx = discretization_->computeD2vDx2(i, j);
			double v2Dy = discretization_->computeDv2Dy(i, j);
			double vDyy = discretization_->computeD2vDy2(i, j);
			double uvDx = discretization_->computeDuvDx(i, j);

			double convection = v2Dy + uvDx;
			double diffusion = vDxx + vDyy;

			g(i, j) = v(i, j) + timeStepWidth_ * (invRe * diffusion - convection + settings_.g[1]);
		}
	}
}

void Simulation::setRightHandSide() {
	const auto& dx = discretization_->dx();
	const auto& dy = discretization_->dy();

	const auto& f = discretization_->f();
	const auto& g = discretization_->g();
	auto& rhs = discretization_->rhs();

	const double invTimeStep = 1.0 / timeStepWidth_;
	const double invDx = 1.0 / dx;
	const double invDy = 1.0 / dy;

	for (int j = rhs.beginJ() + 1; j < rhs.endJ() - 1; j++) {
		for (int i = rhs.beginI() + 1; i < rhs.endI() - 1; i++) {
			const double diffF = (f(i, j) - f(i - 1, j)) * invDx;
			const double diffG = (g(i, j) - g(i, j - 1)) * invDy;
			rhs(i, j) = invTimeStep * (diffF + diffG);
		}
	}
}

void Simulation::setVelocities() {
	auto& u = discretization_->u();
	auto& f = discretization_->f();

	for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
		for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
			u(i, j) = f(i, j) - timeStepWidth_ * discretization_->computeDpDx(i, j);
		}
	}

	auto& v = discretization_->v();
	auto& g = discretization_->g();

	for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
		for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
			v(i, j) = g(i, j) - timeStepWidth_ * discretization_->computeDpDy(i, j);
		}
	}
}
