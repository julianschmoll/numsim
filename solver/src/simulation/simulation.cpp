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
#include <cmath>

Simulation::Simulation(const Settings &settings, const std::string &folderName) {
    settings_ = settings;
    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);

    if (partitioning_->onPrimaryRank())
        std::cout << "Initializing parallel Simulation on " << partitioning_->nRanks() << " ranks." << std::endl;

    for (int i = 0; i < 2; ++i) {
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];
    }

    if (settings_.useDonorCell) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Donor Cell." << std::endl;
        discr_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, *partitioning_,
                                                       settings_.alpha);
    } else {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Central Differences." << std::endl;
        discr_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, *partitioning_, 0.0);
    }

    if (settings_.pressureSolver == IterSolverType::SOR) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using (Red-Black) SOR solver." << std::endl;
        pressureSolver_ =
                std::make_unique<RedBlackSolver>(discr_, partitioning_, settings_);
    } else if (settings_.pressureSolver == IterSolverType::GaussSeidel) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using (Red-Black) Gauss-Seidel solver." << std::endl;
        settings_.omega = 1;
        pressureSolver_ = std::make_unique<RedBlackSolver>(discr_, partitioning_, settings_);
    } else {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Conjugate Gradient solver." << std::endl;
        pressureSolver_ = std::make_unique<ConjugateGradientSolver>(discr_, partitioning_, settings_);
    }

    partitioning_->barrier();
    DEBUG(partitioning_->printPartitioningInfo());

    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discr_, *partitioning_, folderName);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discr_, *partitioning_, folderName);
}

void Simulation::writeOutput(const double currentTime, const int currentSec, const int lastSec) const {
    if (currentSec > lastSec && !settings_.generateTrainingData) [[unlikely]] {
        outputWriterParaview_->writeFile(currentTime);
    }
}

void Simulation::run() {

    const auto start = std::chrono::high_resolution_clock::now();

    double currentTime = 0.0;

    std::vector uv = {&discr_->u(), &discr_->v()};
    std::vector fg = {&discr_->f(), &discr_->g()};

    setBoundaryUV(currentTime);
    setBoundaryFG();

    while (currentTime < settings_.endTime) {
        TimeSteppingInfo timeSteppingInfo = computeTimeStepWidth(currentTime);
        timeStepWidth_ = timeSteppingInfo.timeStepWidth;

        setPreliminaryVelocities();
        setBoundaryFG();
        partitioning_->exchange(fg);

        setRightHandSide();
        pressureSolver_->solve();

        setVelocities();
        partitioning_->exchange(uv);
        setBoundaryUV(currentTime);

        const int lastSec = static_cast<int>(currentTime);
        currentTime += timeStepWidth_;
        const int currentSec = static_cast<int>(currentTime);

        printConsoleInfo(currentTime, timeSteppingInfo);
        DEBUG(outputWriterText_->writeFile(currentTime));

        writeOutput(currentTime, currentSec, lastSec);
    }

    if (settings_.generateTrainingData)
        outputWriterParaview_->writeFile(currentTime);

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
            DEBUG(std::cout << " -- max(u,v) = " << std::fixed << std::setprecision(2) << timeSteppingInfo.maxVelocity
                            << ", ")
            DEBUG(std::cout << "dt = " << std::fixed << std::setprecision(4) << timeSteppingInfo.timeStepWidth << "\n")
            DEBUG(std::cout << " -- div constraint = " << std::fixed << std::setprecision(2)
                            << timeSteppingInfo.convectiveConstraint << ", ")
            DEBUG(std::cout << "conectivity constraint = " << std::fixed << std::setprecision(2)
                            << timeSteppingInfo.convectiveConstraint << "\n");
            lastProgress10th += 10;
        }
    }
}

void Simulation::setBoundaryUV(double currentTime) {
    auto &u = discr_->u();
    auto &v = discr_->v();

    double speedVariance = settings_.dirichletAmplitude * sin(settings_.dirichletFrequency * (settings_.dirichletTimeShift + currentTime) * M_PI);

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        switch (settings_.boundaryBottom) {
            case BoundaryType::InflowNoSlip: {
                const auto uBottom = settings_.dirichletBcBottom[0] + speedVariance * settings_.dirichletBcBottom[0];
                const auto vBottom = settings_.dirichletBcBottom[1] + speedVariance * settings_.dirichletBcBottom[1];

                for (int i = discr_->beginI(u); i < discr_->endI(u); ++i) {
                    u(i, discr_->beginJ(u)) = 2.0 * uBottom - u(i, discr_->beginJ(u) + 1);
                }

                for (int i = discr_->beginI(v); i < discr_->endI(v); ++i) {
                    v(i, discr_->beginJ(v)) = vBottom;
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = discr_->beginI(u); i < discr_->endI(u); ++i) {
                    u(i, discr_->beginJ(u)) = u(i, discr_->beginJ(u) + 1);
                }

                for (int i = discr_->beginI(v); i < discr_->endI(v); ++i) {
                    v(i, discr_->beginJ(v)) = v(i, discr_->beginJ(v) + 1);
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        switch (settings_.boundaryTop) {
            case BoundaryType::InflowNoSlip: {
                const auto uTop = settings_.dirichletBcTop[0] + speedVariance * settings_.dirichletBcTop[0];
                const auto vTop = settings_.dirichletBcTop[1] + speedVariance * settings_.dirichletBcTop[1];

                for (int i = discr_->beginI(u); i < discr_->endI(u); ++i) {
                    u(i, discr_->endJ(u) - 1) = 2.0 * uTop - u(i, discr_->endJ(u) - 2);
                }

                for (int i = discr_->beginI(v); i < discr_->endI(v); ++i) {
                    v(i, discr_->endJ(v) - 1) = vTop;
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = discr_->beginI(u); i < discr_->endI(u); ++i) {
                    u(i, discr_->endJ(u) - 1) = u(i, discr_->endJ(u) - 2);
                }

                for (int i = discr_->beginI(v); i < discr_->endI(v); ++i) {
                    v(i, discr_->endJ(v) - 1) = v(i, discr_->endJ(v) - 2);
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        switch (settings_.boundaryLeft) {
            case BoundaryType::InflowNoSlip: {
                const auto uLeft = settings_.dirichletBcLeft[0] + speedVariance * settings_.dirichletBcLeft[0];
                const auto vLeft = settings_.dirichletBcLeft[1] + speedVariance * settings_.dirichletBcLeft[1];

                for (int j = discr_->beginJ(u); j < discr_->endJ(u); ++j) {
                    u(discr_->beginI(u), j) = uLeft;
                }

                for (int j = discr_->beginJ(v); j < discr_->endJ(v); ++j) {
                    v(discr_->beginI(v), j) = 2.0 * vLeft - v(discr_->beginI(v) + 1, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = discr_->beginJ(u); j < discr_->endJ(u); ++j) {
                    u(discr_->beginI(u), j) = u(discr_->beginI(u) + 1, j);
                }

                for (int j = discr_->beginJ(v); j < discr_->endJ(v); ++j) {
                    v(discr_->beginI(v), j) = v(discr_->beginI(v) + 1, j);
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        switch (settings_.boundaryRight) {
            case BoundaryType::InflowNoSlip: {
                const auto uRight = settings_.dirichletBcRight[0] + speedVariance * settings_.dirichletBcRight[0];
                const auto vRight = settings_.dirichletBcRight[1] + speedVariance * settings_.dirichletBcRight[1];

                for (int j = discr_->beginJ(u); j < discr_->endJ(u); ++j) {
                    u(discr_->endI(u) - 1, j) = uRight;
                }

                for (int j = discr_->beginJ(v); j < discr_->endJ(v); ++j) {
                    v(discr_->endI(v) - 1, j) = 2.0 * vRight - v(discr_->endI(v) - 2, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = discr_->beginJ(u); j < discr_->endJ(u); ++j) {
                    u(discr_->endI(u) - 1, j) = u(discr_->endI(u) - 2, j);
                }

                for (int j = discr_->beginJ(v); j < discr_->endJ(v); ++j) {
                    v(discr_->endI(v) - 1, j) = v(discr_->endI(v) - 2, j);
                }
                break;
            }
        }
    }
}

void Simulation::setBoundaryFG() {
    auto &f = discr_->f();
    auto &g = discr_->g();

    auto &u = discr_->u();
    auto &v = discr_->v();

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        switch (settings_.boundaryBottom) {
            case BoundaryType::InflowNoSlip: {
                for (int i = discr_->beginI(f); i < discr_->endI(f); ++i) {
                    f(i, discr_->beginJ(f)) = u(i, discr_->beginJ(f));
                }

                for (int i = discr_->beginI(g); i < discr_->endI(g); ++i) {
                    g(i, discr_->beginJ(g)) = v(i, discr_->beginJ(g));
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = discr_->beginI(f); i < discr_->endI(f); ++i) {
                    f(i, discr_->beginJ(f)) = f(i, discr_->beginJ(f) + 1);
                }

                for (int i = discr_->beginI(g); i < discr_->endI(g); ++i) {
                    g(i, discr_->beginJ(g)) = g(i, discr_->beginJ(g) + 1);
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        switch (settings_.boundaryTop) {
            case BoundaryType::InflowNoSlip: {
                for (int i = discr_->beginI(f); i < discr_->endI(f); ++i) {
                    f(i, discr_->endJ(f) - 1) = u(i, discr_->endJ(f) - 1);
                }

                for (int i = discr_->beginI(g); i < discr_->endI(g); ++i) {
                    g(i, discr_->endJ(g) - 1) = v(i, discr_->endJ(g) - 1);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = discr_->beginI(f); i < discr_->endI(f); ++i) {
                    f(i, discr_->endJ(f) - 1) = f(i, discr_->endJ(f) - 2);
                }

                for (int i = discr_->beginI(g); i < discr_->endI(g); ++i) {
                    g(i, discr_->endJ(g) - 1) = g(i, discr_->endJ(g) - 2);
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        switch (settings_.boundaryLeft) {
            case BoundaryType::InflowNoSlip: {
                for (int j = discr_->beginJ(f); j < discr_->endJ(f); ++j) {
                    f(discr_->beginI(f), j) = u(discr_->beginI(f), j);
                }

                for (int j = discr_->beginJ(g); j < discr_->endJ(g); ++j) {
                    g(discr_->beginI(g), j) = v(discr_->beginI(g), j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = discr_->beginJ(f); j < discr_->endJ(f); ++j) {
                    f(discr_->beginI(f), j) = f(discr_->beginI(f) + 1, j);
                }

                for (int j = discr_->beginJ(g); j < discr_->endJ(g); ++j) {
                    g(discr_->beginI(g), j) = g(discr_->beginI(g) + 1, j);
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        switch (settings_.boundaryRight) {
            case BoundaryType::InflowNoSlip: {
                for (int j = discr_->beginJ(f); j < discr_->endJ(f); ++j) {
                    f(discr_->endI(f) - 1, j) = u(discr_->endI(f) - 1, j);
                }

                for (int j = discr_->beginJ(g); j < discr_->endJ(g); ++j) {
                    g(discr_->endI(g) - 1, j) = v(discr_->endI(g) - 1, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = discr_->beginJ(f); j < discr_->endJ(f); ++j) {
                    f(discr_->endI(f) - 1, j) = f(discr_->endI(f) - 2, j);
                }

                for (int j = discr_->beginJ(g); j < discr_->endJ(g); ++j) {
                    g(discr_->endI(g) - 1, j) = g(discr_->endI(g) - 2, j);
                }
                break;
            }
        }
    }
}

TimeSteppingInfo Simulation::computeTimeStepWidth(double currentTime) {
    const double uMaxLocal = discr_->u().absMax();
    const double vMaxLocal = discr_->v().absMax();

    double uMaxGlobal{};
    double vMaxGlobal{};
    MPI_Allreduce(&uMaxLocal, &uMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&vMaxLocal, &vMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double hx = discr_->dx();
    const double hy = discr_->dy();

    const double dx2 = hx * hx;
    const double dy2 = hy * hy;
    const double diffusive = (settings_.re / 2.0) * ((dx2 * dy2) / (dx2 + dy2));

    // we need to handle velocities being 0
    const double convectiveU = (uMaxGlobal > 0) ? hx / uMaxGlobal : std::numeric_limits<double>::max();
    const double convectiveV = (vMaxGlobal > 0) ? hy / vMaxGlobal : std::numeric_limits<double>::max();

    double dt = std::min({diffusive, convectiveU, convectiveV}) * settings_.tau;

    dt = std::min(dt, settings_.maximumDt);

    if (currentTime + dt > settings_.endTime) {
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

    auto &u = discr_->u();
    auto &v = discr_->v();
    auto &f = discr_->f();
    auto &g = discr_->g();

    for (int j = discr_->beginJ(u) + 1; j < discr_->endJ(u) - 1; j++) {
        for (int i = discr_->beginI(u) + 1; i < discr_->endI(u) - 1; i++) {
            const double uDxx = discr_->computeD2uDx2(i, j);
            const double u2Dx = discr_->computeDu2Dx(i, j);
            const double uDyy = discr_->computeD2uDy2(i, j);
            const double uvDy = discr_->computeDuvDy(i, j);

            const double convection = u2Dx + uvDy;
            const double diffusion = uDxx + uDyy;

            f(i, j) = u(i, j) + timeStepWidth_ * (invRe * diffusion - convection + settings_.g[0]);
        }
    }

    for (int j = discr_->beginJ(v) + 1; j < discr_->endJ(v) - 1; j++) {
        for (int i = discr_->beginI(v) + 1; i < discr_->endI(v) - 1; i++) {
            const double vDxx = discr_->computeD2vDx2(i, j);
            const double v2Dy = discr_->computeDv2Dy(i, j);
            const double vDyy = discr_->computeD2vDy2(i, j);
            const double uvDx = discr_->computeDuvDx(i, j);

            const double convection = v2Dy + uvDx;
            const double diffusion = vDxx + vDyy;

            g(i, j) = v(i, j) + timeStepWidth_ * (invRe * diffusion - convection + settings_.g[1]);
        }
    }
}

void Simulation::setRightHandSide() {
    const auto &dx = discr_->dx();
    const auto &dy = discr_->dy();

    const auto &f = discr_->f();
    const auto &g = discr_->g();
    auto &rhs = discr_->rhs();

    const double invTimeStep = 1.0 / timeStepWidth_;
    const double invDx = 1.0 / dx;
    const double invDy = 1.0 / dy;

    for (int j = discr_->beginJ(rhs) + 1; j < discr_->endJ(rhs) - 1; j++) {
        for (int i = discr_->beginI(rhs) + 1; i < discr_->endI(rhs) - 1; i++) {
            const double diffF = (f(i, j) - f(i - 1, j)) * invDx;
            const double diffG = (g(i, j) - g(i, j - 1)) * invDy;
            rhs(i, j) = invTimeStep * (diffF + diffG);
        }
    }
}

void Simulation::setVelocities() {
    auto &u = discr_->u();
    auto &f = discr_->f();

    for (int j = discr_->beginJ(u) + 1; j < discr_->endJ(u) - 1; j++) {
        for (int i = discr_->beginI(u) + 1; i < discr_->endI(u) - 1; i++) {
            u(i, j) = f(i, j) - timeStepWidth_ * discr_->computeDpDx(i, j);
        }
    }

    auto &v = discr_->v();
    auto &g = discr_->g();

    for (int j = discr_->beginJ(v) + 1; j < discr_->endJ(v) - 1; j++) {
        for (int i = discr_->beginI(v) + 1; i < discr_->endI(v) - 1; i++) {
            v(i, j) = g(i, j) - timeStepWidth_ * discr_->computeDpDy(i, j);
        }
    }
}

void Simulation::saveState() {
    DEBUG(std::cout << "Simulation::saveState" << std::endl);
    uCheckpoint_ = discr_->u();
    vCheckpoint_ = discr_->v();
    pCheckpoint_ = discr_->p();
}

void Simulation::reloadLastState() {
    DEBUG(std::cout << "Simulation::reloadLastState" << std::endl);
    discr_->u() = uCheckpoint_;
    discr_->v() = vCheckpoint_;
    discr_->p() = pCheckpoint_;
}
