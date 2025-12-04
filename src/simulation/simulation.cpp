#include "simulation.h"
#include "macros.h"
#include "outputWriter/outputWriterParaviewParallel.h"
#include "outputWriter/outputWriterTextParallel.h"
#include "settings.h"
#include "simulation/pressureSolver/conjugateGradientSolver.h"
#include "simulation/pressureSolver/redBlackSolver.h"
#include <chrono>
#include <iostream>
#include <ostream>

Simulation::Simulation(const Settings &settings) {
    settings_ = settings;
    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);

    if (partitioning_->onPrimaryRank())
        std::cout << "Initializing parallel Simulation on" << partitioning_->nRanks() << " ranks." << std::endl;

    for (int i = 0; i < 2; ++i) {
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];
    }

    if (settings_.useDonorCell) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Donor Cell." << std::endl;
        discOps_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, *partitioning_, settings_.alpha);
    } else {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Central Differences." << std::endl;
        discOps_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, *partitioning_, 0.0);
    }

    if (settings_.pressureSolver == IterSolverType::SOR) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using (Red-Black) SOR solver." << std::endl;
        pressureSolver_ =
            std::make_unique<RedBlackSolver>(discOps_, partitioning_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else if (settings_.pressureSolver == IterSolverType::GaussSeidel) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using (Red-Black) Gauss-Seidel solver." << std::endl;
        pressureSolver_ = std::make_unique<RedBlackSolver>(discOps_, partitioning_, settings_.epsilon, settings_.maximumNumberOfIterations, 1);
    } else {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Conjugate Gradient solver." << std::endl;
        pressureSolver_ = std::make_unique<ConjugateGradientSolver>(discOps_, partitioning_, settings_.epsilon, settings_.maximumNumberOfIterations);
    }

    partitioning_->barrier();
    partitioning_->printPartitioningInfo();

    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discOps_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discOps_, *partitioning_);
}

void Simulation::run() {

    auto start = std::chrono::high_resolution_clock::now();

    double currentTime = 0.0;

    std::vector uv = {&discOps_->u(), &discOps_->v()};
    std::vector fg = {&discOps_->f(), &discOps_->g()};

    setBoundaryUV();
    setBoundaryFG();

    while (currentTime < settings_.endTime) {
        setPreliminaryVelocities();
        partitioning_->nonBlockingExchange(fg);
        setBoundaryUV();

        TimeSteppingInfo timeSteppingInfo = computeTimeStepWidth(currentTime);
        timeStepWidth_ = timeSteppingInfo.timeStepWidth;

        partitioning_->waitForAllMPIRequests();
        setRightHandSide();
        partitioning_->exchange(discOps_->rhs());

        pressureSolver_->solve();

        setVelocities();
        partitioning_->nonBlockingExchange(uv);

        const int lastSec = static_cast<int>(currentTime);
        currentTime += timeStepWidth_;
        const int currentSec = static_cast<int>(currentTime);
        const bool writeOutput = (currentSec > lastSec);

        printConsoleInfo(currentTime, timeSteppingInfo);

        partitioning_->waitForAllMPIRequests();
        DEBUG(outputWriterText_->writeFile(currentTime));
        if (writeOutput) [[unlikely]] {
            outputWriterParaview_->writeFile(currentTime);
        }
    }

    partitioning_->barrier();

    if (partitioning_->onPrimaryRank()) {
        auto end = std::chrono::high_resolution_clock::now();
        uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Simulation finished in " << ms << "ms" << std::endl;
    }
}

void Simulation::printConsoleInfo(double currentTime, const TimeSteppingInfo &timeSteppingInfo) const {
    static double lastProgress10th = 0;
    if (partitioning_->onPrimaryRank()) {
        double progress = currentTime / settings_.endTime * 100;
        if (progress >= lastProgress10th) {
            std::cout << "Progress: ";
            std::cout << std::fixed << std::setw(5) << std::setprecision(1) << progress << "%: ";
            std::cout << std::fixed << std::setw(5) << std::setprecision(2) << currentTime << "s / ";
            std::cout << std::fixed << std::setprecision(2) << settings_.endTime << "s\n";
            DEBUG(std::cout << " -- max(u,v) = " << std::fixed << std::setprecision(2) << timeSteppingInfo.maxVelocity << ", ")
            DEBUG(std::cout << "dt = " << std::fixed << std::setprecision(4) << timeSteppingInfo.timeStepWidth << "\n")
            DEBUG(std::cout << " -- div constraint = " << std::fixed << std::setprecision(2) << timeSteppingInfo.convectiveConstraint << ", ")
            DEBUG(std::cout << "conectivity constraint = " << std::fixed << std::setprecision(2) << timeSteppingInfo.convectiveConstraint << "\n");
            lastProgress10th += 10;
        }
    }
}

void Simulation::setBoundaryUV() {
    auto &u = discOps_->u();
    auto &v = discOps_->v();

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        const auto uBottom = settings_.dirichletBcBottom[0];
        const auto vBottom = settings_.dirichletBcBottom[1];

        for (int i = u.beginI(); i < u.endI(); ++i) {
            u(i, u.beginJ()) = 2.0 * uBottom - u(i, u.beginJ() + 1);
        }

        for (int i = v.beginI(); i < v.endI(); ++i) {
            v(i, v.beginJ()) = vBottom;
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        const auto uTop = settings_.dirichletBcTop[0];
        const auto vTop = settings_.dirichletBcTop[1];

        for (int i = u.beginI(); i < u.endI(); ++i) {
            u(i, u.endJ() - 1) = 2.0 * uTop - u(i, u.endJ() - 2);
        }

        for (int i = v.beginI(); i < v.endI(); ++i) {
            v(i, v.endJ() - 1) = vTop;
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        const auto uLeft = settings_.dirichletBcLeft[0];
        const auto vLeft = settings_.dirichletBcLeft[1];

        for (int j = u.beginJ(); j < u.endJ(); ++j) {
            u(u.beginI(), j) = uLeft;
        }

        for (int j = v.beginJ(); j < v.endJ(); ++j) {
            v(v.beginI(), j) = 2.0 * vLeft - v(v.beginI() + 1, j);
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        const auto uRight = settings_.dirichletBcRight[0];
        const auto vRight = settings_.dirichletBcRight[1];

        for (int j = u.beginJ(); j < u.endJ(); ++j) {
            u(u.endI() - 1, j) = uRight;
        }

        for (int j = v.beginJ(); j < v.endJ(); ++j) {
            v(v.endI() - 1, j) = 2.0 * vRight - v(v.endI() - 2, j);
        }
    }
}

void Simulation::setBoundaryFG() {
    auto &f = discOps_->f();
    auto &g = discOps_->g();

    auto &u = discOps_->u();
    auto &v = discOps_->v();

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {

        for (int i = f.beginI(); i < f.endI(); ++i) {
            f(i, f.beginJ()) = u(i, f.beginJ());
        }

        for (int i = g.beginI(); i < g.endI(); ++i) {
            g(i, g.beginJ()) = v(i, g.beginJ());
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {

        for (int i = f.beginI(); i < f.endI(); ++i) {
            f(i, f.endJ() - 1) = u(i, f.endJ() - 1);
        }

        for (int i = g.beginI(); i < g.endI(); ++i) {
            g(i, g.endJ() - 1) = v(i, g.endJ() - 1);
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {

        for (int j = f.beginJ(); j < f.endJ(); ++j) {
            f(f.beginI(), j) = u(f.beginI(), j);
        }

        for (int j = g.beginJ(); j < g.endJ(); ++j) {
            g(g.beginI(), j) = v(g.beginI(), j);
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {

        for (int j = f.beginJ(); j < f.endJ(); ++j) {
            f(f.endI() - 1, j) = u(f.endI() - 1, j);
        }

        for (int j = g.beginJ(); j < g.endJ(); ++j) {
            g(g.endI() - 1, j) = v(g.endI() - 1, j);
        }
    }
}

TimeSteppingInfo Simulation::computeTimeStepWidth(double currentTime) {
    const double uMaxLocal = discOps_->u().absMax();
    const double vMaxLocal = discOps_->v().absMax();

    double uMaxGlobal{};
    double vMaxGlobal{};
    MPI_Allreduce(&uMaxLocal, &uMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&vMaxLocal, &vMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double hx = discOps_->dx();
    const double hy = discOps_->dy();

    const double dx2 = hx * hx;
    const double dy2 = hy * hy;
    const double diffusive = (settings_.re / 2.0) * ((dx2 * dy2) / (dx2 + dy2));

    // we need to handle velocities being 0
    const double convectiveU = (uMaxGlobal > 0) ? hx / uMaxGlobal : std::numeric_limits<double>::max();
    const double convectiveV = (vMaxGlobal > 0) ? hy / vMaxGlobal : std::numeric_limits<double>::max();

    double dt = std::min({diffusive, convectiveU, convectiveV}) * settings_.tau;

    dt = std::min(dt, settings_.maximumDt);

    if (currentTime + dt > settings_.endTime) [[unlikely]] {
        dt = settings_.endTime - currentTime;
    }

    TimeSteppingInfo info{.convectiveConstraint = std::min(convectiveU, convectiveV),
                          .diffusiveConstraint = diffusive,
                          .maxVelocity = std::max(uMaxGlobal, vMaxGlobal),
                          .timeStepWidth = dt};

    return info;
}

void Simulation::setPreliminaryVelocities() {
    const double invRe = 1.0 / settings_.re;

    auto &u = discOps_->u();
    auto &v = discOps_->v();
    auto &f = discOps_->f();
    auto &g = discOps_->g();

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
    const auto &dx = discOps_->dx();
    const auto &dy = discOps_->dy();

    const auto &f = discOps_->f();
    const auto &g = discOps_->g();
    auto &rhs = discOps_->rhs();

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
    auto &u = discOps_->u();
    auto &f = discOps_->f();

    for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
        for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
            u(i, j) = f(i, j) - timeStepWidth_ * discOps_->computeDpDx(i, j);
        }
    }

    auto &v = discOps_->v();
    auto &g = discOps_->g();

    for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
        for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
            v(i, j) = g(i, j) - timeStepWidth_ * discOps_->computeDpDy(i, j);
        }
    }
}
