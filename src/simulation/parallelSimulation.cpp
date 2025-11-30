#include "parallelSimulation.h"

#include "grid/staggeredGrid.h"
#include "settings.h"
#include "simulation/discreteOperators.h"
#include "pressureSolver/pressureSolverSerial.h"
#include "macros.h"
#include "outputWriter/outputWriterParaviewParallel.h"
#include "outputWriter/outputWriterTextParallel.h"
#include "pressureSolver/redBlack.h"

#include <memory>
#include <chrono>

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
        discOps_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, settings_.alpha);
    } else {
        std::cout << "Using Central Differences." << std::endl;
        discOps_ = std::make_unique<DiscreteOperators>(partitioning_->nCellsLocal(), meshWidth_, 0.0);
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

    setBoundaryFG();

    while (currentTime < settings_.endTime) {
        setBoundaryUV();

        computeTimeStepWidth();

        DEBUG(std::cout << "**************************************" << std::endl);
        DEBUG(std::cout << "Computing Timestep " << counter++);
        DEBUG(const auto expectedTimesteps = static_cast<unsigned int>(counter + (settings_.endTime - currentTime) / timeStepWidth_ - 1));
        DEBUG(std::cout << " of " << expectedTimesteps << " estimated Timesteps." << std::endl);

        if (currentTime + timeStepWidth_ > settings_.endTime) {
            timeStepWidth_ = settings_.endTime - currentTime;
        }
        currentTime += timeStepWidth_;

        DEBUG(std::cout << "Current Time: " << currentTime << "s");
        DEBUG(std::cout << "/" << settings_.endTime << "s");
        DEBUG(const double progress = (currentTime / settings_.endTime) * 100.0);
        DEBUG(std::cout << " (" << std::fixed << std::setprecision(2) << progress << "%)" << std::endl);

        setPreliminaryVelocities();
        setRightHandSide();

        pressureSolver_->solve(); //ToDo: add communication for pressure halo cells after every iteration

        DEBUG(std::cout << "Solved pressure" << std::endl);

        setVelocities();

        DEBUG(std::cout << "Set velocities" << std::endl);

        exchangeVelocities();
        //ToDo: Exchange velocity halo cells

        if (partitioning_->ownRankNo() == ROOT_RANK) {
            outputWriterParaview_->writeFile(currentTime);
            outputWriterText_->writeFile(currentTime);
        } else {
            //ToDo: Put local values in global grid
            // is that actually a todo or isn't that done by the writer?
            //ToDo: Send global grid to ROOT_RANK
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Simulation finished in " << ms << "ms" << std::endl;
}

void ParallelSimulation::setBoundaryUV() {
    auto &u = discOps_->u();
    auto &v = discOps_->v();

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
}

void ParallelSimulation::setBoundaryFG() {
    auto &f = discOps_->f();
    auto &g = discOps_->g();

    if (partitioning_->ownPartitionContainsBoundary(Direction::Left)) {
        const auto fLeft = settings_.dirichletBcLeft[0];
        const auto gLeft = settings_.dirichletBcLeft[1];

        for (int j = f.beginJ(); j < f.endJ(); ++j) {
            f(f.beginI(), j) = fLeft;
        }
        for (int j = g.beginJ(); j < g.endJ(); ++j) {
            g(g.beginI(), j) = gLeft;
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Right)) {
        const auto fRight = settings_.dirichletBcRight[0];
        const auto gRight = settings_.dirichletBcRight[1];

        for (int j = f.beginJ(); j < f.endJ(); ++j) {
            f(f.endI() -1, j) = fRight;
        }
        for (int j = g.beginJ(); j < g.endJ(); ++j) {
            g(g.endI() - 1, j) = gRight;
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Top)) {
        const auto fTop = settings_.dirichletBcTop[0];
        const auto gTop = settings_.dirichletBcTop[1];

        for (int i = f.beginI(); i< f.endI(); ++i) {
            f(i, f.endJ()-1) = fTop;
        }
        for (int i = g.beginI(); i < g.endI(); ++i) {
            g(i, g.endJ()-1) = gTop;
        }
    }

    if (partitioning_->ownPartitionContainsBoundary(Direction::Bottom)) {
        const auto fBottom = settings_.dirichletBcBottom[0];
        const auto gBottom = settings_.dirichletBcBottom[1];

        for (int i = f.beginI(); i< f.endI(); ++i) {
            f(i, f.beginJ()) = fBottom;
        }
        for (int i = g.beginI(); i < g.endI(); ++i) {
            g(i, g.beginJ()) = gBottom;
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

void ParallelSimulation::exchangeVelocities() {
    int ownPartWidth = partitioning_->nCellsLocal()[0];
    int ownPartHeight = partitioning_->nCellsLocal()[1];

    int leftRankNo = partitioning_->neighborRankNo(Direction::Left);
    int rightRankNo = partitioning_->neighborRankNo(Direction::Right);
    int bottomRankNo = partitioning_->neighborRankNo(Direction::Bottom);
    int topRankNo = partitioning_->neighborRankNo(Direction::Top);

    DataField &u = discOps_->u();
    DataField &v = discOps_->v();

    MPI_Request requestULeft, requestVLeft, requestURight, requestVRight, requestUTop, requestVTop, requestUBottom, requestVBottom;

    //Receive requests should cancel automatically if source rank is MPI_PROC_NULL

    //LEFT SIDE
    std::vector<double> leftURecvBuffer(ownPartHeight), leftVRecvBuffer(ownPartHeight);
    MPI_Irecv(leftURecvBuffer.data(), ownPartHeight, MPI_DOUBLE, leftRankNo, U_TAG, MPI_COMM_WORLD, &requestULeft);
    MPI_Irecv(leftVRecvBuffer.data(), ownPartHeight, MPI_DOUBLE, leftRankNo, V_TAG, MPI_COMM_WORLD, &requestVLeft);

    if (leftRankNo != MPI_PROC_NULL) {
        std::vector<double> leftUSendBuffer(ownPartHeight), leftVSendBuffer(ownPartHeight);
        for (int j = 0; j < ownPartHeight; j++) {
            leftUSendBuffer[j] = u(0, j);
            leftVSendBuffer[j] = v(0, j);
        }
        MPI_Send(leftUSendBuffer.data(), ownPartHeight, MPI_DOUBLE, leftRankNo, U_TAG, MPI_COMM_WORLD);
        MPI_Send(leftVSendBuffer.data(), ownPartHeight, MPI_DOUBLE, leftRankNo, V_TAG, MPI_COMM_WORLD);
    }


    //RIGHT SIDE
    std::vector<double> rightURecvBuffer(ownPartHeight), rightVRecvBuffer(ownPartHeight);
    MPI_Irecv(rightURecvBuffer.data(), ownPartHeight, MPI_DOUBLE, rightRankNo, U_TAG, MPI_COMM_WORLD, &requestURight);
    MPI_Irecv(rightVRecvBuffer.data(), ownPartHeight, MPI_DOUBLE, rightRankNo, V_TAG, MPI_COMM_WORLD, &requestVRight);

    if (rightRankNo != MPI_PROC_NULL) {
        std::vector<double> rightUSendBuffer(ownPartHeight), rightVSendBuffer(ownPartHeight);
        for (int j = 0; j < ownPartHeight; j++) {
            rightUSendBuffer[j] = u(ownPartWidth - 1, j);
            rightVSendBuffer[j] = v(ownPartWidth - 1, j);
        }
        MPI_Send(rightUSendBuffer.data(), ownPartHeight, MPI_DOUBLE, rightRankNo, U_TAG, MPI_COMM_WORLD);
        MPI_Send(rightVSendBuffer.data(), ownPartHeight, MPI_DOUBLE, rightRankNo, V_TAG, MPI_COMM_WORLD);
    }

    MPI_Request requestUBottomSend, requestUTopSend, requestUBottomRecv, requestUTopRecv;
    MPI_Isend(bottomUSendBuffer.data(), ownPartWidth, MPI_DOUBLE, bottomRankNo, U_TAG, MPI_COMM_WORLD, &requestUBottomSend);
    MPI_Irecv(bottomURecvBuffer.data(), ownPartWidth, MPI_DOUBLE, bottomRankNo, U_TAG, MPI_COMM_WORLD, &requestUBottomRecv);

    MPI_Isend(topUSendBuffer.data(), ownPartWidth, MPI_DOUBLE, topRankNo, U_TAG, MPI_COMM_WORLD, &requestUTopSend);
    MPI_Irecv(topURecvBuffer.data(), ownPartWidth, MPI_DOUBLE, topRankNo, U_TAG, MPI_COMM_WORLD, &requestUTopRecv);


    //u: right --> left
    std::vector<double> rightUSendBuffer(ownPartHeight), leftURecvBuffer(ownPartHeight);

    for (int j = 0; j < ownPartHeight; j++) {
        rightUSendBuffer[j] = u(ownPartWidth - 1, j);
    }

    MPI_Request requestURightSend, requestULeftRecv;
    MPI_Isend(rightUSendBuffer.data(), ownPartHeight, MPI_DOUBLE, rightRankNo, U_TAG, MPI_COMM_WORLD, &requestURightSend);
    MPI_Irecv(leftURecvBuffer.data(), ownPartHeight, MPI_DOUBLE, leftRankNo, U_TAG, MPI_COMM_WORLD, &requestULeftRecv);

    //We need to waitall before pasting the values, since the left side of the current partition is the right side of the left  neighbor partition.
    //Hence, we can't do the directions sequentially like "complete left, then go on with right" etc.
    std::array<MPI_Request, 8> requests = {
        requestULeft, requestVLeft,
        requestURight, requestVRight,
        requestUBottom, requestVBottom,
        requestUTop, requestVTop
    };
    MPI_Waitall(8, requests.data(), MPI_STATUS_IGNORE);

    //Put received values into grid
    if (leftRankNo != MPI_PROC_NULL) {
        for (int j = 0; j < ownPartHeight; j++) {
            u(-1, j) = leftURecvBuffer[j];
            v(-1, j) = leftVRecvBuffer[j];
        }
    }
    if (rightRankNo != MPI_PROC_NULL) {
        for (int j = 0; j < ownPartHeight; j++) {
            u(ownPartWidth, j) = rightURecvBuffer[j];
            v(ownPartWidth, j) = rightVRecvBuffer[j];
        }
    }
    if (bottomRankNo != MPI_PROC_NULL) {
        for (int i = 0; i < ownPartHeight; i++) {
            u(i, -1) = bottomURecvBuffer[i];
            v(i, -1) = bottomVRecvBuffer[i];
        }
    }
    if (topRankNo != MPI_PROC_NULL) {
        for (int i = 0; i < ownPartHeight; i++) {
            u(i, ownPartHeight) = topURecvBuffer[i];
            v(i, ownPartHeight) = topVRecvBuffer[i];
        }
    }

    //ToDo (Bene): Fix errors and refactor this huge duplicating piece of sh*t so it becomes less ugly
}
