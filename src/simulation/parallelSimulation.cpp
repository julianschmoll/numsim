#include "parallelSimulation.h"

#include "grid/dataField.h"
#include "settings.h"
#include "simulation/discreteOperators.h"
#include "macros.h"
#include "outputWriter/outputWriterParaviewParallel.h"
#include "outputWriter/outputWriterTextParallel.h"
#include "pressureSolver/redBlack.h"
#include "simulation/partitioning.h"

#include <memory>
#include <chrono>
#include <ostream>

void ParallelSimulation::initialize(const Settings &settings) {
    std::cout << "Initializing parallel Simulation..." << std::endl;
    settings_ = settings;

    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);

    for (int i = 0; i < 2; ++i) {
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];
    }

    if (settings_.useDonorCell) {
        std::cout << "Using Donor Cell." << std::endl;
        discOps_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, *partitioning_, settings_.alpha);
    } else {
        std::cout << "Using Central Differences." << std::endl;
        discOps_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, *partitioning_, 0.0);
    }

    if (settings_.pressureSolver == IterSolverType::SOR) {
        std::cout << "Using SOR solver." << std::endl;
        pressureSolver_ = std::make_unique<RedBlack>(discOps_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    } else {
        pressureSolver_ = std::make_unique<RedBlack>(discOps_, settings_.epsilon, settings_.maximumNumberOfIterations, 1, partitioning_);
    }

    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discOps_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discOps_, *partitioning_);
}

void ParallelSimulation::run() {

    auto start = std::chrono::high_resolution_clock::now();

    double currentTime = 0.0;

    DEBUG(unsigned int counter = 0);

    //ToDo: Exchange velocity halo cells (necessary only if we do not always start with 0-velocities for inner, non-boundary cells)

    partitioning_->printPartitioningInfo();

    std::vector uv = {&discOps_->u(), &discOps_->v()};
    std::vector fg = {&discOps_->f(), &discOps_->g()};

    setBoundaryFG();

    partitioning_->exchange(fg);

    while (currentTime < settings_.endTime) {
        setBoundaryUV();
        setBoundaryFG(); // TODO: not every iteration but after u, v

        partitioning_->exchange(uv);

        computeTimeStepWidth();

        DEBUG(std::cout << "**************************************" << std::endl);
        DEBUG(std::cout << "Computing Timestep " << counter++);
        DEBUG(const auto expectedTimesteps = static_cast<unsigned int>(counter + (settings_.endTime - currentTime) / timeStepWidth_ - 1));
        DEBUG(std::cout << " of " << expectedTimesteps << " estimated Timesteps." << std::endl);

        if (currentTime + timeStepWidth_ > settings_.endTime) {
            timeStepWidth_ = settings_.endTime - currentTime;
        }

        // ToDo: Do we want to snap to integer values or just write when we crossed one
        const int lastSec = static_cast<int>(currentTime);
        currentTime += timeStepWidth_;
        const int currentSec = static_cast<int>(currentTime);
        const bool writeOutput = (currentSec > lastSec);

        DEBUG(std::cout << "Current Time: " << currentTime << "s");
        DEBUG(std::cout << "/" << settings_.endTime << "s");
        DEBUG(const double progress = (currentTime / settings_.endTime) * 100.0);
        DEBUG(std::cout << " (" << std::fixed << std::setprecision(2) << progress << "%)" << std::endl);

        setPreliminaryVelocities();

        partitioning_->exchange(fg);

        setRightHandSide();

        pressureSolver_->solve();

        setVelocities();
        DEBUG(std::cout << "Set velocities" << std::endl);

        partitioning_->exchange(uv);

        if (true) {
            //outputWriterParaview_->writeFile(currentTime);
            outputWriterText_->writeFile(currentTime);
        }
        break;
    }

    auto end = std::chrono::high_resolution_clock::now();
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Simulation finished in " << ms << "ms" << std::endl;
}

void ParallelSimulation::setBoundaryUV() {
    auto &u = discOps_->u();
    auto &v = discOps_->v();

    if (partitioning_->ownPartitionContainsBoundary(Direction::Bottom)) {
        const auto uBottom = settings_.dirichletBcBottom[0];
        const auto vBottom = settings_.dirichletBcBottom[1];

        for (int i = u.beginI(); i < u.endI(); ++i) {
            u(i, u.beginJ()) = 2.0 * uBottom - u(i, u.beginJ() + 1);
        }
        for (int i = v.beginI(); i < v.endI(); ++i) {
            v(i, v.beginJ()) = vBottom;
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Top)) {
        const auto uTop = settings_.dirichletBcTop[0];
        const auto vTop = settings_.dirichletBcTop[1];

        for (int i = u.beginI(); i < u.endI(); ++i) {
            u(i, u.endJ() - 1) = 2.0 * uTop - u(i, u.endJ() - 2);
        }
        for (int i = v.beginI(); i < v.endI(); ++i) {
            v(i, v.endJ() - 1) = vTop;
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Left)) {
        const auto uLeft = settings_.dirichletBcLeft[0];
        const auto vLeft = settings_.dirichletBcLeft[1];

        for (int j = u.beginJ(); j < u.endJ(); ++j) {
            u(u.beginI(), j) = uLeft;
        }
        for (int j = v.beginJ(); j < v.endJ(); ++j) {
            v(v.beginI(), j) = 2.0 * vLeft - v(v.beginI() + 1, j);
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Right)) {
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

void ParallelSimulation::setBoundaryFG() {
    auto &f = discOps_->f();
    auto &g = discOps_->g();

    auto &u = discOps_->u();
    auto &v = discOps_->v();

    if (partitioning_->ownPartitionContainsBoundary(Direction::Bottom)) {
        const auto fBottom = settings_.dirichletBcBottom[0];
        const auto gBottom = settings_.dirichletBcBottom[1];

        for (int i = f.beginI(); i< f.endI(); ++i) {
            f(i, f.beginJ()) = u(i, f.beginJ());
        }
        for (int i = g.beginI(); i < g.endI(); ++i) {
            g(i, g.beginJ()) = v(i, g.beginJ());
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Top)) {
        const auto fTop = settings_.dirichletBcTop[0];
        const auto gTop = settings_.dirichletBcTop[1];

        for (int i = f.beginI(); i< f.endI(); ++i) {
            f(i, f.endJ()-1) = u(i, f.endJ()-1);
        }
        for (int i = g.beginI(); i < g.endI(); ++i) {
            g(i, g.endJ()-1) = v(i, g.endJ()-1);
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Left)) {
        const auto fLeft = settings_.dirichletBcLeft[0];
        const auto gLeft = settings_.dirichletBcLeft[1];

        for (int j = f.beginJ(); j < f.endJ(); ++j) {
            f(f.beginI(), j) = u(f.beginI(), j);
        }
        for (int j = g.beginJ(); j < g.endJ(); ++j) {
            g(g.beginI(), j) = v(g.beginI(), j);
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Right)) {
        const auto fRight = settings_.dirichletBcRight[0];
        const auto gRight = settings_.dirichletBcRight[1];

        for (int j = f.beginJ(); j < f.endJ(); ++j) {
            f(f.endI() - 1, j) = u(f.endI() - 1, j);
        }
        for (int j = g.beginJ(); j < g.endJ(); ++j) {
            g(g.endI() - 1, j) = v(g.endI() - 1, j);
        }
    }
}

void ParallelSimulation::computeTimeStepWidth() {
    const double uMaxLocal = discOps_->u().absMax();
    const double vMaxLocal = discOps_->v().absMax();

    double uMaxGlobal;
    double vMaxGlobal;
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

    const double dt = std::min({diffusive, convectiveU, convectiveV}) * settings_.tau;

    timeStepWidth_ = std::min(dt, settings_.maximumDt);
}
