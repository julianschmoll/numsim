#include "simulation.h"
#include "macros.h"
#include "outputWriter/outputWriterParaviewParallel.h"
#include "outputWriter/outputWriterTextParallel.h"
#include "settings.h"
#include "simulation/pressureSolver/conjugateGradientSolver.h"
#include "simulation/pressureSolver/redBlackSolver.h"
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>

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
        discOps_ = std::make_unique<DiscreteOperators>(settings_, *partitioning_, settings_.alpha);
    } else {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Central Differences." << std::endl;
        discOps_ = std::make_unique<DiscreteOperators>(settings_, *partitioning_, 0.0);
    }

    if (settings_.pressureSolver == IterSolverType::SOR) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using (Red-Black) SOR solver." << std::endl;
        pressureSolver_ = std::make_unique<RedBlackSolver>(discOps_, partitioning_, settings_);
    } else if (settings_.pressureSolver == IterSolverType::GaussSeidel) {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using (Red-Black) Gauss-Seidel solver." << std::endl;
        settings_.omega = 1;
        pressureSolver_ = std::make_unique<RedBlackSolver>(discOps_, partitioning_, settings_);
    } else {
        if (partitioning_->onPrimaryRank())
            std::cout << " -- Using Conjugate Gradient solver." << std::endl;
        pressureSolver_ = std::make_unique<ConjugateGradientSolver>(discOps_, partitioning_, settings_);
    }

    partitioning_->barrier();
    DEBUG(partitioning_->printPartitioningInfo());

    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discOps_, *partitioning_, folderName);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discOps_, *partitioning_, folderName);
}

void Simulation::writeOutput(const int currentSec, const int lastSec, bool always) const {
    if (always || currentSec > lastSec) {
        outputWriterParaview_->writeFile(currentTime_);
        DEBUG(outputWriterText_->writeFile(currentTime_));
    }
}

void Simulation::setDisplacements(const std::vector<double> &topDisplacements, const std::vector<double> &bottomDisplacements) {
    // TODO: Unsere Displacements haben 2 Einträge mehr als die innere Zellenanzahl.. Das wäre für Parallelisierung wichtig.
    int n = discOps_->displacementsTop_.size();
    assert(n == static_cast<int>(topDisplacements.size()));
    assert(n == static_cast<int>(bottomDisplacements.size()));
    discOps_->displacementsTop_ = topDisplacements;
    discOps_->displacementsBottom_ = bottomDisplacements;

    double domainHeight = settings_.physicalSize[1];

    // TODO: Wir könnten die Geschwindigkeitsränder eventuell sogar hier aktualisieren
    // dann wäre kopieren und speichern der displacements unnötig.
    // Zu tatsächlichem Rand hinzufügen (lokal, da wir unsere globale Position kennen)
    for (int i = 0; i < n; i++) {
        discOps_->topBoundaryPosition_[i] += discOps_->displacementsTop_[i];
        discOps_->bottomBoundaryPosition_[i] += discOps_->displacementsBottom_[i];

        discOps_->topBoundaryPosition_[i] = std::min(domainHeight, discOps_->topBoundaryPosition_[i]);
        discOps_->bottomBoundaryPosition_[i] = std::max(0.0, discOps_->bottomBoundaryPosition_[i]);
    }
}

void Simulation::run() {

    const auto start = std::chrono::high_resolution_clock::now();

    currentTime_ = 0.0;

    std::vector uv = {&discOps_->u(), &discOps_->v()};
    std::vector fg = {&discOps_->f(), &discOps_->g()};

    DataField &p = discOps_->p();

    size_t n = settings_.nCells[0] + 2;
    std::vector<double> topDisplacements(n, 0);
    std::vector<double> bottomDisplacements(n, 0);
    for (size_t i = 0; i < discOps_->bottomBoundaryPosition_.size(); i++) {
        bottomDisplacements[i] = (i <= n / 2) ? 0.04 * i : 0.04 * (n - i);
    }
    initializeDisplacements(topDisplacements, bottomDisplacements);

    while (currentTime_ < settings_.endTime) {
        TimeSteppingInfo timeSteppingInfo = computeTimeStepWidth();
        timeStepWidth_ = timeSteppingInfo.timeStepWidth;

        const int lastSec = static_cast<int>(currentTime_);
        advanceFluidSolver(timeStepWidth_);
        const int currentSec = static_cast<int>(currentTime_);

        printConsoleInfo(timeSteppingInfo);
        writeOutput(currentSec, lastSec, false);
    }

    if (settings_.generateTrainingData)
        outputWriterParaview_->writeFile(currentTime_);

    partitioning_->barrier();

    if (partitioning_->onPrimaryRank()) {
        auto end = std::chrono::high_resolution_clock::now();
        uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Simulation finished in " << ms << "ms" << std::endl;
    }
}

void Simulation::initializeDisplacements(std::vector<double> &displacements) {
    setOuterVelocityBoundaries();

    constexpr int meshDim = 3;
    const int n = settings_.nCells[0];

    assert(displacements.size() == 2 * n * meshDim);

    // +1 because it's a flat array with x,y,z,x,y,z...
    const int topOffset = n * meshDim + 1;
    constexpr int bottomOffset = 1;

    for (int i = 0, idxTop = topOffset, idxBottom = bottomOffset; i < n; ++i, idxTop += meshDim, idxBottom += meshDim) {
        discOps_->topBoundaryPosition_[i + 1] = settings_.physicalSize[1] + displacements[idxTop];
        discOps_->bottomBoundaryPosition_[i + 1] = displacements[idxBottom];
    }
    discOps_->updateStructureCells(0); // 0 weil Wert egal
    setStructureBoundaries();
}

void Simulation::initializeDisplacements(const std::vector<double> &topDisplacements, const std::vector<double> &bottomDisplacements) {
    setOuterVelocityBoundaries();
    for (size_t i = 0; i < settings_.nCells[0] + 2; i++) {
        discOps_->topBoundaryPosition_[i] = settings_.physicalSize[1] + topDisplacements[i];
        discOps_->bottomBoundaryPosition_[i] = bottomDisplacements[i];
    }
    discOps_->updateStructureCells(0); // 0 weil Wert egal
    setStructureBoundaries();
}

void Simulation::advanceFluidSolver(double dt) {
    std::vector uv = {&discOps_->u(), &discOps_->v()};
    std::vector fg = {&discOps_->f(), &discOps_->g()};

    setOuterVelocityBoundaries();
    setStructureBoundaries();
    setPreliminaryVelocities();

    partitioning_->exchange(fg);

    setRightHandSide();
    pressureSolver_->solve(discOps_->p());
    setVelocities();
    partitioning_->exchange(uv);

    calculateForces();

    currentTime_ += dt;
}

void Simulation::updateSolid() {
    auto &q = discOps_->q();

    bool correctionRequired = discOps_->updateStructureCells(timeStepWidth_);

    if (correctionRequired) {
        setOuterVelocityBoundaries();
        setStructureBoundaries();
        setPreliminaryVelocities();

        setRightHandSide();
        pressureSolver_->solve(q);

        setVelocities();
    }

    DEBUG(std::cout << discOps_->structure_ << "\n");
}

void Simulation::printConsoleInfo(const TimeSteppingInfo &timeSteppingInfo) const {
    static double lastProgressStep = 0;
    if (partitioning_->onPrimaryRank()) {
        double progress = currentTime_ / settings_.endTime * 100;
        if (progress >= lastProgressStep) {
            std::cout << "Progress: ";
            std::cout << std::fixed << std::setw(5) << std::setprecision(1) << progress << "%: ";
            std::cout << std::fixed << std::setw(5) << std::setprecision(2) << currentTime_ << "s / ";
            std::cout << std::fixed << std::setprecision(2) << settings_.endTime << "s\n";
            DEBUG(std::cout << " -- max(u,v) = " << std::fixed << std::setprecision(2) << timeSteppingInfo.maxVelocity << ", ")
            DEBUG(std::cout << "dt = " << std::fixed << std::setprecision(4) << timeSteppingInfo.timeStepWidth << "\n")
            DEBUG(std::cout << " -- div constraint = " << std::fixed << std::setprecision(2) << timeSteppingInfo.convectiveConstraint << ", ")
            DEBUG(std::cout << "conectivity constraint = " << std::fixed << std::setprecision(2) << timeSteppingInfo.convectiveConstraint << "\n");
            lastProgressStep += 5;
        }
    }
}

// must be called after updating structure
void Simulation::setStructureBoundaries() {
    auto &u = discOps_->u();
    auto &v = discOps_->v();
    auto &g = discOps_->g();
    auto &f = discOps_->f();

    // TODO: we do not consider partition boundaries here

    if (settings_.boundaryBottom == BoundaryType::Coupled) {
        // u bottom border
        for (int i = u.minI() + 1; i <= u.maxI() - 1; ++i) {
            for (int j = u.minJ(); j <= u.maxJ() - 1; ++j) {
                if (discOps_->isFluid(i, j)) {
                    break;
                }
                const bool leftFluid = discOps_->isFluid(i - 1, j);
                const bool rightFluid = discOps_->isFluid(i + 1, j);
                const bool topFluid = discOps_->isFluid(i, j + 1);

                if (topFluid && leftFluid && rightFluid) {
                    f(i, j) = u(i, j) = 0;
                    f(i - 1, j) = u(i - 1, j) = 0;
                } else if (topFluid && leftFluid) {
                    f(i, j) = u(i, j) = -u(i, j + 1);
                    f(i - 1, j) = u(i - 1, j) = 0;
                } else if (topFluid && rightFluid) {
                    f(i, j) = u(i, j) = 0;
                } else if (topFluid) {
                    f(i, j) = u(i, j) = -u(i, j + 1);
                } else if (leftFluid) {
                    f(i - 1, j) = u(i - 1, j) = 0;
                } else if (rightFluid) {
                    f(i, j) = u(i, j) = 0;
                }
            }
        }
        // v bottom border
        for (int i = v.minI() + 1; i <= v.maxI() - 1; ++i) {
            const double dCenter = discOps_->bottomDisplacement(i) / timeStepWidth_;
            const double dLeft = discOps_->bottomDisplacement(i - 1) / timeStepWidth_;
            const double dRight = discOps_->bottomDisplacement(i + 1) / timeStepWidth_;

            for (int j = v.minJ(); j <= v.maxJ() - 1; ++j) {
                if (discOps_->isFluid(i, j)) {
                    break;
                }
                const bool leftFluid = discOps_->isFluid(i - 1, j);
                const bool rightFluid = discOps_->isFluid(i + 1, j);
                const bool topFluid = discOps_->isFluid(i, j + 1);

                if (topFluid && leftFluid && rightFluid) {
                    g(i, j) = v(i, j) = dCenter;
                } else if (topFluid && leftFluid) {
                    g(i, j) = v(i, j) = dCenter;
                } else if (topFluid && rightFluid) {
                    g(i, j) = v(i, j) = dCenter;
                } else if (topFluid) {
                    g(i, j) = v(i, j) = dCenter;
                } else if (leftFluid) {
                    g(i, j) = v(i, j) = dCenter + dRight - v(i - 1, j);
                } else if (rightFluid) {
                    g(i, j) = v(i, j) = dLeft + dCenter - v(i + i, j);
                }
            }
        }
    }

    if (settings_.boundaryTop == BoundaryType::Coupled) {
        // u top border
        for (int i = u.minI() + 1; i <= u.maxI() - 1; ++i) {
            for (int j = u.maxJ(); j <= u.minJ(); --j) {
                if (discOps_->isFluid(i, j)) {
                    break;
                }
                const bool leftFluid = discOps_->isFluid(i - 1, j);
                const bool rightFluid = discOps_->isFluid(i + 1, j);
                const bool bottomFluid = discOps_->isFluid(i, j - 1);

                if (bottomFluid && leftFluid && rightFluid) {
                    f(i, j) = u(i, j) = 0;
                    f(i - 1, j) = u(i - 1, j) = 0;
                } else if (bottomFluid && leftFluid) {
                    f(i, j) = u(i, j) = -u(i, j - 1);
                    f(i - 1, j) = u(i - 1, j) = 0;
                } else if (bottomFluid && rightFluid) {
                    f(i, j) = u(i, j) = 0;
                } else if (bottomFluid) {
                    f(i, j) = u(i, j) = -u(i, j - 1);
                } else if (leftFluid) {
                    f(i - 1, j) = u(i - 1, j) = 0;
                } else if (rightFluid) {
                    f(i, j) = u(i, j) = 0;
                }
            }
        }
        // v top border
        for (int i = v.minI() + 1; i <= v.maxI() - 1; ++i) {
            const double dCenter = discOps_->topDisplacement(i) / timeStepWidth_;
            const double dLeft = discOps_->topDisplacement(i - 1) / timeStepWidth_;
            const double dRight = discOps_->topDisplacement(i + 1) / timeStepWidth_;

            for (int j = v.maxJ() + 1; j >= v.minJ() + 1; --j) {
                if (discOps_->isFluid(i, j)) {
                    break;
                }
                const bool leftFluid = discOps_->isFluid(i - 1, j);
                const bool rightFluid = discOps_->isFluid(i + 1, j);
                const bool bottomFluid = discOps_->isFluid(i, j - 1);

                if (bottomFluid && leftFluid && rightFluid) {
                    g(i, j - 1) = v(i, j - 1) = dCenter;
                } else if (bottomFluid && leftFluid) {
                    g(i, j - 1) = v(i, j - 1) = dCenter;
                } else if (bottomFluid && rightFluid) {
                    g(i, j - 1) = v(i, j - 1) = dCenter;
                } else if (bottomFluid) {
                    g(i, j - 1) = v(i, j - 1) = dCenter;
                } else if (leftFluid) {
                    g(i, j - 1) = v(i, j - 1) = dCenter + dRight - v(i - 1, j - 1);
                } else if (rightFluid) {
                    g(i, j - 1) = v(i, j - 1) = dLeft + dCenter - v(i + i, j - 1);
                }
            }
        }
    }
}

void Simulation::setOuterVelocityBoundaries() {
    auto &u = discOps_->u();
    auto &v = discOps_->v();
    auto &f = discOps_->f();
    auto &g = discOps_->g();

    double speedVariance = settings_.dirichletAmplitude * sin(settings_.dirichletFrequency * (settings_.dirichletTimeShift + currentTime_) * M_PI);
    bool inflowActive = currentTime_ >= settings_.startBurst_ && currentTime_ <= settings_.endBurst_;

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        switch (settings_.boundaryBottom) {
        case BoundaryType::InflowNoSlip:
        case BoundaryType::Coupled: {
            const auto uBottom = inflowActive ? settings_.dirichletBcBottom[0] + speedVariance * settings_.dirichletBcBottom[0] : 0.0;
            const auto vBottom = inflowActive ? settings_.dirichletBcBottom[1] + speedVariance * settings_.dirichletBcBottom[1] : 0.0;

            for (int i = u.beginI(); i < u.endI(); ++i) {
                f(i, u.beginJ()) = u(i, u.beginJ()) = 2.0 * uBottom - u(i, u.beginJ() + 1);
            }

            for (int i = v.beginI(); i < v.endI(); ++i) {
                g(i, v.beginJ()) = v(i, v.beginJ()) = vBottom;
            }
            break;
        }

        case BoundaryType::Outflow: {
            for (int i = u.beginI(); i < u.endI(); ++i) {
                f(i, u.beginJ()) = u(i, u.beginJ()) = u(i, u.beginJ() + 1);
            }

            for (int i = v.beginI(); i < v.endI(); ++i) {
                g(i, v.beginJ()) = v(i, v.beginJ()) = v(i, v.beginJ() + 1);
            }
            break;
        }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        switch (settings_.boundaryTop) {
        case BoundaryType::InflowNoSlip:
        case BoundaryType::Coupled: {
            const auto uTop = inflowActive ? settings_.dirichletBcTop[0] + speedVariance * settings_.dirichletBcTop[0] : 0.0;
            const auto vTop = inflowActive ? settings_.dirichletBcTop[1] + speedVariance * settings_.dirichletBcTop[1] : 0.0;

            for (int i = u.beginI(); i < u.endI(); ++i) {
                f(i, u.maxJ()) = u(i, u.maxJ()) = 2.0 * uTop - u(i, u.endJ() - 2);
            }

            for (int i = v.beginI(); i < v.endI(); ++i) {
                g(i, v.maxJ()) = v(i, v.maxJ()) = vTop;
            }
            break;
        }

        case BoundaryType::Outflow: {
            for (int i = u.beginI(); i < u.endI(); ++i) {
                f(i, u.maxJ()) = u(i, u.maxJ()) = u(i, u.maxJ() - 1);
            }

            for (int i = v.beginI(); i < v.endI(); ++i) {
                g(i, v.maxJ()) = v(i, v.maxJ()) = v(i, v.maxJ() - 1);
            }
            break;
        }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        switch (settings_.boundaryLeft) {
        case BoundaryType::InflowNoSlip: {
            const auto uLeft = inflowActive ? settings_.dirichletBcLeft[0] + speedVariance * settings_.dirichletBcLeft[0] : 0.0;
            const auto vLeft = inflowActive ? settings_.dirichletBcLeft[1] + speedVariance * settings_.dirichletBcLeft[1] : 0.0;

            for (int j = u.beginJ(); j < u.endJ(); ++j) {
                f(u.minI(), j) = u(u.minI(), j) = uLeft;
            }

            for (int j = v.beginJ(); j < v.endJ(); ++j) {
                g(v.minI(), j) = v(v.minI(), j) = 2.0 * vLeft - v(v.minI() + 1, j);
            }
            break;
        }

        case BoundaryType::Outflow: {
            for (int j = u.beginJ(); j < u.endJ(); ++j) {
                f(u.minI(), j) = u(u.minI(), j) = u(u.minI() + 1, j);
            }

            for (int j = v.beginJ(); j < v.endJ(); ++j) {
                g(u.minI(), j) = v(v.minI(), j) = v(v.minI() + 1, j);
            }
            break;
        }

        case BoundaryType::Coupled:
            assert(false);
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        switch (settings_.boundaryRight) {
        case BoundaryType::InflowNoSlip: {
            const auto uRight = inflowActive ? settings_.dirichletBcRight[0] + speedVariance * settings_.dirichletBcRight[0] : 0.0;
            const auto vRight = inflowActive ? settings_.dirichletBcRight[1] + speedVariance * settings_.dirichletBcRight[1] : 0.0;

            for (int j = u.beginJ(); j < u.endJ(); ++j) {
                f(u.maxI(), j) = u(u.maxI(), j) = uRight;
            }

            for (int j = v.beginJ(); j < v.endJ(); ++j) {
                g(v.maxI(), j) = v(v.maxI(), j) = 2.0 * vRight - v(v.endI() - 2, j);
            }
            break;
        }

        case BoundaryType::Outflow: {
            for (int j = u.beginJ(); j < u.endJ(); ++j) {
                f(u.maxI(), j) = u(u.maxI(), j) = u(u.maxI() - 1, j);
            }

            for (int j = v.beginJ(); j < v.endJ(); ++j) {
                g(v.maxI(), j) = v(v.maxI(), j) = v(v.maxI() - 1, j);
            }
            break;
        }

        case BoundaryType::Coupled:
            assert(false);
        }
    }
}

TimeSteppingInfo Simulation::computeTimeStepWidth() {
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

    if (currentTime_ + dt > settings_.endTime) {
        dt = settings_.endTime - currentTime_;
    }

    TimeSteppingInfo info{.convectiveConstraint = std::min(convectiveU, convectiveV),
                          .diffusiveConstraint = diffusive,
                          .maxVelocity = std::max(uMaxGlobal, vMaxGlobal),
                          .timeStepWidth = dt};

    return info;
}

void Simulation::setTimeStepWidth(double dt) {
    timeStepWidth_ = dt;
}

void Simulation::setPreliminaryVelocities() {
    const double invRe = 1.0 / settings_.re;

    auto &u = discOps_->u();
    auto &v = discOps_->v();
    auto &f = discOps_->f();
    auto &g = discOps_->g();

    for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
        for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
            if (discOps_->isSolid(i, j) || discOps_->isSolid(i + 1, j))
                continue;
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
            if (discOps_->isSolid(i, j) || discOps_->isSolid(i, j + 1))
                continue;
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
            if (discOps_->isSolid(i, j))
                continue;
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
            if (discOps_->isSolid(i, j) || discOps_->isSolid(i + 1, j))
                continue;
            u(i, j) = f(i, j) - timeStepWidth_ * discOps_->computeDpDx(i, j);
        }
    }

    auto &v = discOps_->v();
    auto &g = discOps_->g();

    for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
        for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
            if (discOps_->isSolid(i, j) || discOps_->isSolid(i, j + 1))
                continue;
            v(i, j) = g(i, j) - timeStepWidth_ * discOps_->computeDpDy(i, j);
        }
    }
}

void Simulation::correctVelocities() {
    auto &u = discOps_->u();
    auto &f = discOps_->f();

    for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) {
        for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
            if (discOps_->isSolid(i, j) || discOps_->isSolid(i + 1, j))
                continue;
            u(i, j) = f(i, j) - timeStepWidth_ * discOps_->computeDqDx(i, j);
        }
    }

    auto &v = discOps_->v();
    auto &g = discOps_->g();

    for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) {
        for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
            if (discOps_->isSolid(i, j) || discOps_->isSolid(i, j + 1))
                continue;
            v(i, j) = g(i, j) - timeStepWidth_ * discOps_->computeDqDy(i, j);
        }
    }
}

void Simulation::calculateForces() {
    auto &v = discOps_->v();
    auto &p = discOps_->p();

    const double invRe = 1.0 / settings_.re;
    const double dx = discOps_->dx();
    const double dy = discOps_->dy();

    // collect forces orthogonal to top boundary
    for (int i = v.beginI() + 1; i < v.endI() - 1; ++i) {
        for (int j = v.endJ() - 1; j > v.beginJ(); --j) {
            if (discOps_->isFluid(i, j)) {
                const double fBoundary = invRe * discOps_->computeDvDy(i, j) - v(i, j) * v(i, j) - (p(i, j + 1) + p(i, j)) / 2;
                const double fBefore = invRe * discOps_->computeDvDy(i, j - 1) - v(i, j - 1) * v(i, j - 1) - (p(i, j) + p(i, j - 1)) / 2;

                const double lowerCellEdge = discOps_->globalDomainPosJ(j);
                const double boundary = discOps_->topBoundaryPosition(i);
                const double a = (boundary - lowerCellEdge) / dy;

                discOps_->topF(i) = -dx * ((1 - a) * fBefore + a * fBoundary);
                break;
            }
        }
    }

    // collect forces orthogonal to bottom boundary
    for (int i = v.beginI() + 1; i < v.endI() - 1; ++i) {
        for (int j = v.beginJ(); j < v.endJ() - 1; ++j) {
            if (discOps_->isFluid(i, j)) {
                const double fBoundary = invRe * discOps_->computeDvDy(i, j) - v(i, j - 1) * v(i, j - 1) - (p(i, j - 1) + p(i, j)) / 2;
                const double fBefore = invRe * discOps_->computeDvDy(i, j + 1) - v(i, j) * v(i, j) - (p(i, j) + p(i, j + 1)) / 2;

                const double lowerCellEdge = discOps_->globalDomainPosJ(j);
                const double boundary = discOps_->bottomBoundaryPosition(i);
                const double a = (boundary - lowerCellEdge) / dy;

                discOps_->bottomF(i) = dx * ((1 - a) * fBoundary + a * fBefore);
                break;
            }
        }
    }
}

void Simulation::test() {
    std::cout << "Simulation::test: Start" << std::endl;
    discOps_->test(settings_);
    std::cout << "Simulation::test: End" << std::endl;
}

// -----------------------------------------------//
// --------- preCICE Adapter Medthods ------------//
// -----------------------------------------------//

// ToDo: States are overkill, need to figure out what,
// apparently uv and dt was not enough
void Simulation::saveState() {
    DEBUG(std::cout << "\nSimulation::saveState()\n" << std::endl);
    uCheckpoint_ = discOps_->u();
    vCheckpoint_ = discOps_->v();
    pCheckpoint_ = discOps_->p();
    qCheckpoint_ = discOps_->q();
    fCheckpoint_ = discOps_->f();
    gCheckpoint_ = discOps_->g();
    timeStepWidthCheckpoint_ = timeStepWidth_;
    checkpointTime_ = currentTime_;
    displacementsTopCheckpoint_ = discOps_->displacementsTop_;
    displacementsBottomCheckpoint_ = discOps_->displacementsBottom_;
    topBoundaryPositionCheckpoint_ = discOps_->topBoundaryPosition_;
    bottomBoundaryPositionCheckpoint_ = discOps_->bottomBoundaryPosition_;
    structureCheckpoint_ = discOps_->structure_;
}

void Simulation::reloadLastState() {
    DEBUG(std::cout << "\nSimulation::reloadLastState()\n" << std::endl);
    discOps_->u() = uCheckpoint_;
    discOps_->v() = vCheckpoint_;
    discOps_->p() = pCheckpoint_;
    discOps_->q() = qCheckpoint_;
    discOps_->f() = fCheckpoint_;
    discOps_->g() = gCheckpoint_;
    currentTime_ = checkpointTime_;
    timeStepWidth_ = timeStepWidthCheckpoint_;
    discOps_->topBoundaryPosition_ = topBoundaryPositionCheckpoint_;
    discOps_->bottomBoundaryPosition_ = bottomBoundaryPositionCheckpoint_;
    discOps_->displacementsTop_ = displacementsTopCheckpoint_;
    discOps_->displacementsBottom_ = displacementsBottomCheckpoint_;
    discOps_->structure_ = structureCheckpoint_;
}

std::shared_ptr<Partitioning> Simulation::getPartitioning() const noexcept {
    return partitioning_;
}

void Simulation::getForces(std::vector<double> &forces) {
    const int meshDim = 3;
    const int n = settings_.nCells[0];

    constexpr int coordOffset = 1; // x: 0, y: 1, z: 2
    const int topOffset = coordOffset + n * meshDim;
    const int bottomOffset = coordOffset;

    forces.assign(forces.size(), 0.0);

    for (int i = 0, idxTop = topOffset, idxBottom = bottomOffset; i < n; ++i, idxTop += meshDim, idxBottom += meshDim) {
        forces[idxTop] = discOps_->topF(i);
        forces[idxBottom] = discOps_->bottomF(i);
    }
}

void Simulation::setDisplacements(std::vector<double> &displacements) {

    constexpr int meshDim = 3;
    const int n = settings_.nCells[0];

    std::vector top(n + 2, 0.0);
    std::vector bottom(n + 2, 0.0);
    assert(displacements.size() == 2 * n * meshDim);

    constexpr int bottomOffset = 1;

    // +1 because it's a flat array with x,y,z,x,y,z...
    const int topOffset = n * meshDim + 1;

    for (int i = 0, idxTop = topOffset, idxBottom = bottomOffset; i < n; ++i, idxTop += meshDim, idxBottom += meshDim) {
        top[i + 1] = displacements[idxTop];
        bottom[i + 1] = displacements[idxBottom];
    }
    setDisplacements(top, bottom);
}

double Simulation::getCurrentTime() {
    return currentTime_;
}

void Simulation::writeLineToFile(const std::string &filePath, const std::string &line) {
    std::ofstream file(filePath, std::ios::app); // std::ios::app = append mode

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return;
    }

    file << line << '\n';
    file.close();
}
