#include "simulation.h"

#include <utility>

Simulation::Simulation(Settings settings) : settings_(std::move(settings)) {
	std::cout << "Initializing Simulation..." << std::endl;

	for (int i = 0; i < 2; ++i) { meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i]; }

	if (settings_.useDonorCell) {
		std::cout << "Using Donor Cell." << std::endl;
		discOps_ = std::make_unique<discreteOperators>(settings_.nCells, meshWidth_, settings_.alpha);
	}
	else {
		std::cout << "Using Central Differences." << std::endl;
		discOps_ = std::make_unique<discreteOperators>(settings_.nCells, meshWidth_, 0.0);
	}

	if (settings_.pressureSolver == "SOR") {
		std::cout << "Using SOR solver." << std::endl;
		pressureSolver_ = std::make_unique<PressureSolver>(discOps_, settings_.epsilon,
		                                                   settings_.maximumNumberOfIterations, settings_.omega);
	}
	else if (settings_.pressureSolver == "GaussSeidel") {
		pressureSolver_ = std::make_unique<PressureSolver>(discOps_, settings_.epsilon,
		                                                   settings_.maximumNumberOfIterations, 1);
	}
	else { throw std::runtime_error("Unknown pressure solver."); }

	outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discOps_);
	outputWriterText_ = std::make_unique<OutputWriterText>(discOps_);
}

void Simulation::run() {
#ifndef NDEBUG
	unsigned int counter = 0;
#endif
	double currentTime = 0.0;

	while (currentTime < settings_.endTime) {
		setBoundaryValues();

		computeTimeStepWidth();

#ifndef NDEBUG
		std::cout << "**************************************" << std::endl;
		std::cout << "Computing Timestep " << counter++;
		const unsigned int expectedTimesteps = static_cast<int>(settings_.endTime/timeStepWidth_+1);
		std::cout << " of " << expectedTimesteps << " estimated Timesteps." << std::endl;
#endif

		if (currentTime + timeStepWidth_ > settings_.endTime) { timeStepWidth_ = settings_.endTime - currentTime; }
		currentTime += timeStepWidth_;

#ifndef NDEBUG
		std::cout << "Current Time: " << currentTime << "s";
		std::cout << "/" << settings_.endTime << "s";
		const double progress = (currentTime / settings_.endTime) * 100.0;
		std::cout << " (" << std::fixed << std::setprecision(2) << progress << "%)" << std::endl;
#endif

		setPreliminaryVelocities();

		setRightHandSide();

		pressureSolver_->solve();

		setVelocities();

		outputWriterParaview_->writeFile(currentTime);
		outputWriterText_->writeFile(currentTime);
	}

	std::cout << "Simulation finished." << std::endl;
}

void Simulation::setBoundaryValues() {
	const auto uBottom = settings_.dirichletBcBottom[0];
	const auto uTop = settings_.dirichletBcTop[0];
	const auto uLeft = settings_.dirichletBcLeft[0];
	const auto uRight = settings_.dirichletBcRight[0];

	const auto vBottom = settings_.dirichletBcBottom[1];
	const auto vTop = settings_.dirichletBcTop[1];
	const auto vLeft = settings_.dirichletBcLeft[1];
	const auto vRight = settings_.dirichletBcRight[1];

	auto& u = discOps_->u();
	auto& v = discOps_->v();
	auto& f = discOps_->f();
	auto& g = discOps_->g();

	// We expect u and f having the same size
	for (int i = u.beginI(); i < u.endI(); ++i) {
		const double valBeginJ = 2.0 * uBottom - u(i, u.beginJ() + 1);
		u(i, u.beginJ()) = valBeginJ;
		f(i, f.beginJ()) = valBeginJ;

		const double valEndJ = 2.0 * uTop - u(i, u.endJ() - 2);
		u(i, u.endJ() - 1) = valEndJ;
		f(i, f.endJ() - 1) = valEndJ;
	}

	// We expect v and g having the same size
	for (int i = v.beginI(); i < v.endI(); ++i) {
		v(i, v.beginJ()) = vBottom;
		g(i, g.beginJ()) = vBottom;
		v(i, v.endJ() - 1) = vTop;
		g(i, g.endJ() - 1) = vTop;
	}

	for (int j = u.beginJ(); j < u.endJ(); ++j) {
		u(u.beginI(), j) = uLeft;
		f(f.beginI(), j) = uLeft;
		u(u.endI() - 1, j) = uRight;
		f(f.endI() - 1, j) = uRight;
	}

	for (int j = v.beginJ(); j < v.endJ(); ++j) {
		const double valBeginI = 2.0 * vLeft - v(v.beginI() + 1, j);
		v(v.beginI(), j) = valBeginI;
		g(g.beginI(), j) = valBeginI;

		const double valEndI = 2.0 * vRight - v(v.endI() - 2, j);
		v(v.endI() - 1, j) = valEndI;
		g(g.endI() - 1, j) = valEndI;
	}
}

void Simulation::computeTimeStepWidth() {
	const double uMax = discOps_->u().absMax();
	const double vMax = discOps_->v().absMax();

	const double hx = discOps_->dx();
	const double hy = discOps_->dy();

	const double diffusive = (((hx * hx) + (hy * hy)) / 2.0) * (settings_.re / 2.0);

	// we need to handle velocities being 0
	const double convectiveU = (uMax > 0) ? hx / uMax : std::numeric_limits<double>::max();
	const double convectiveV = (vMax > 0) ? hy / vMax : std::numeric_limits<double>::max();

	const double dt = std::min({diffusive, convectiveU, convectiveV}) * settings_.tau;

	timeStepWidth_ = std::min(dt, settings_.maximumDt);
}

void Simulation::setPreliminaryVelocities() {
	const double invRe = 1.0 / settings_.re;

	auto& u = discOps_->u();
	auto& v = discOps_->v();
	auto& f = discOps_->f();
	auto& g = discOps_->g();

	for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
		for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
			const double uDxx = discOps_->computeD2uDx2(i, j);
			const double u2Dx = discOps_->computeDu2Dx(i, j);
			const double uDyy = discOps_->computeD2uDy2(i, j);
			const double uvDy = discOps_->computeDuvDy(i, j);

			const double convection = u2Dx + uvDy;
			const double diffusion = uDxx + uDyy;

			f(i, j) = u(i, j) + timeStepWidth_ * (invRe * diffusion - convection + settings_.g[0]);
		}
	}

	for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
		for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
			const double vDxx = discOps_->computeD2vDx2(i, j);
			const double v2Dy = discOps_->computeDv2Dy(i, j);
			const double vDyy = discOps_->computeD2vDy2(i, j);
			const double uvDx = discOps_->computeDuvDx(i, j);

			const double convection = v2Dy + uvDx;
			const double diffusion = vDxx + vDyy;

			g(i, j) = v(i, j) + timeStepWidth_ * (invRe * diffusion - convection + settings_.g[1]);
		}
	}
}

void Simulation::setRightHandSide() {
	const auto& dx = discOps_->dx();
	const auto& dy = discOps_->dy();

	const auto& f = discOps_->f();
	const auto& g = discOps_->g();
	auto& rhs = discOps_->rhs();

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
	auto& u = discOps_->u();
	auto& f = discOps_->f();

	for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
		for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
			u(i, j) = f(i, j) - timeStepWidth_ * discOps_->computeDpDx(i, j);
		}
	}

	auto& v = discOps_->v();
	auto& g = discOps_->g();

	for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
		for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
			v(i, j) = g(i, j) - timeStepWidth_ * discOps_->computeDpDy(i, j);
		}
	}
}
