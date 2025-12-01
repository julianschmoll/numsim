#include "parallelSimulation.h"

#include "settings.h"
#include "simulation/discreteOperators.h"
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
        setRightHandSide();

        pressureSolver_->solve();

        setVelocities();

        DEBUG(std::cout << "Set velocities" << std::endl);

        exchangeVelocities();
        //ToDo: Exchange velocity halo cells

        if (writeOutput) {
            outputWriterParaview_->writeFile(currentTime);
            outputWriterText_->writeFile(currentTime);
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
            f(f.endI() - 1, j) = fRight;
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
    //ToDo: Use MPI_comm from partitioning

    int leftRankNo = partitioning_->neighborRankNo(Direction::Left);
    int rightRankNo = partitioning_->neighborRankNo(Direction::Right);
    int bottomRankNo = partitioning_->neighborRankNo(Direction::Bottom);
    int topRankNo = partitioning_->neighborRankNo(Direction::Top);

    DataField &u = discOps_->u();
    DataField &v = discOps_->v();

    int uWidth = u.endI() - u.beginI();
    int uHeight = u.endJ() - u.beginJ();

    int vWidth = v.endI() - v.beginI();
    int vHeight = v.endJ() - v.beginJ();

    int ownRankNumber = partitioning_->ownRankNo();

    /*
    for (int i = u.beginI() + 1; i < u.endI() - 1; i++) {
        for (int j = u.beginJ() + 1; j < u.endJ() - 1; j++) u(i,j) = ownRankNumber;
    }
    for (int i = v.beginI() + 1; i < v.endI() - 1; i++) {
        for (int j = v.beginJ() + 1; j < v.endJ() - 1; j++) v(i,j) = ownRankNumber;
    }*/

    //Receive requests should cancel automatically if source rank is MPI_PROC_NULL

    std::vector<double> leftVSendBuffer(vHeight), leftVRecvBuffer(vHeight); //TODO Eckrandwerte?
    std::vector<double> rightVSendBuffer(vHeight),  rightVRecvBuffer(vHeight);
    std::vector<double> topVSendBuffer(vWidth), topVRecvBuffer(vWidth);
    std::vector<double> bottomVSendBuffer(vWidth), bottomVRecvBuffer(vWidth);


    // v: left <--> right

    for (int j = v.beginJ(); j < v.endJ(); j++) { // TODO: begin, end?
        leftVSendBuffer[j + 1] = v(v.beginI() + 1, j);
        rightVSendBuffer[j + 1] = v(v.endI() - 2, j);
    }

    MPI_Request requestVLeftSend, requestVRightSend, requestVLeftRecv, requestVRightRecv;
    MPI_Isend(leftVSendBuffer.data(), leftVSendBuffer.size(), MPI_DOUBLE, leftRankNo, V_TAG, MPI_COMM_WORLD, &requestVLeftSend);
    MPI_Irecv(leftVRecvBuffer.data(), leftVRecvBuffer.size(), MPI_DOUBLE, leftRankNo, V_TAG, MPI_COMM_WORLD, &requestVLeftRecv);

    MPI_Isend(rightVSendBuffer.data(), rightVSendBuffer.size(), MPI_DOUBLE, rightRankNo, V_TAG, MPI_COMM_WORLD, &requestVRightSend);
    MPI_Irecv(rightVRecvBuffer.data(), rightVRecvBuffer.size(), MPI_DOUBLE, rightRankNo, V_TAG, MPI_COMM_WORLD, &requestVRightRecv);


    // v: top <--> bottom
    for (int i = v.beginI(); i < v.endI(); i++) {
        topVSendBuffer[i + 1] = v(i, v.endJ() - 2);
        bottomVSendBuffer[i + 1] = v(i, v.beginJ() + 1);
    }
    MPI_Request requestVTopSend, requestVBottomSend, requestVBottomRecv, requestVTopRecv;
    MPI_Isend(topVSendBuffer.data(), topVSendBuffer.size(), MPI_DOUBLE, topRankNo, V_TAG, MPI_COMM_WORLD, &requestVTopSend);
    MPI_Irecv(bottomVRecvBuffer.data(), bottomVRecvBuffer.size(), MPI_DOUBLE, bottomRankNo, V_TAG, MPI_COMM_WORLD, &requestVBottomRecv);
    MPI_Isend(bottomVSendBuffer.data(), bottomVSendBuffer.size(), MPI_DOUBLE, bottomRankNo, V_TAG, MPI_COMM_WORLD, &requestVBottomSend);
    MPI_Irecv(topVRecvBuffer.data(), topVRecvBuffer.size(), MPI_DOUBLE, topRankNo, V_TAG, MPI_COMM_WORLD, &requestVTopRecv);


    //u: bottom <--> top
    std::vector<double> bottomUSendBuffer(uWidth), bottomURecvBuffer(uWidth);
    std::vector<double> topUSendBuffer(uWidth), topURecvBuffer(uWidth);

    for (int i = u.beginI(); i < u.endI(); i++) {
        bottomUSendBuffer[i + 1] = u(i, u.beginJ() + 1);
        topUSendBuffer[i + 1] = u(i, u.endJ() - 2);
    }

    MPI_Request requestUBottomSend, requestUTopSend, requestUBottomRecv, requestUTopRecv;
    MPI_Isend(bottomUSendBuffer.data(), bottomUSendBuffer.size(), MPI_DOUBLE, bottomRankNo, U_TAG, MPI_COMM_WORLD, &requestUBottomSend);
    MPI_Irecv(bottomURecvBuffer.data(), bottomURecvBuffer.size(), MPI_DOUBLE, bottomRankNo, U_TAG, MPI_COMM_WORLD, &requestUBottomRecv);

    MPI_Isend(topUSendBuffer.data(), topUSendBuffer.size(), MPI_DOUBLE, topRankNo, U_TAG, MPI_COMM_WORLD, &requestUTopSend);
    MPI_Irecv(topURecvBuffer.data(), topURecvBuffer.size(), MPI_DOUBLE, topRankNo, U_TAG, MPI_COMM_WORLD, &requestUTopRecv);


    //u: right <--> left
    std::vector<double> rightUSendBuffer(uHeight), leftUSendBuffer(uHeight), leftURecvBuffer(uHeight), rightURecvBuffer(uHeight);

    for (int j = u.beginJ(); j < u.endJ(); j++) {
        rightUSendBuffer[j + 1] = u(u.endI() - 2, j);
        leftUSendBuffer[j + 1] = u(u.beginI() + 1 , j);
    }

    MPI_Request requestURightSend, requestULeftRecv, requestULeftSend, requestURightRecv;
    MPI_Isend(rightUSendBuffer.data(), rightUSendBuffer.size(), MPI_DOUBLE, rightRankNo, U_TAG, MPI_COMM_WORLD, &requestURightSend);
    MPI_Irecv(leftURecvBuffer.data(), leftURecvBuffer.size(), MPI_DOUBLE, leftRankNo, U_TAG, MPI_COMM_WORLD, &requestULeftRecv);
    MPI_Isend(leftUSendBuffer.data(), leftUSendBuffer.size(), MPI_DOUBLE, leftRankNo, U_TAG, MPI_COMM_WORLD, &requestULeftSend);
    MPI_Irecv(rightURecvBuffer.data(), rightURecvBuffer.size(), MPI_DOUBLE, leftRankNo, U_TAG, MPI_COMM_WORLD, &requestURightRecv);

    //We need to waitall before pasting the values, since the left side of the current partition is the right side of the left  neighbor partition.
    //Hence, we can't do the directions sequentially like "complete left, then go on with right" etc.
    std::array<MPI_Request, 12> requests = {
        requestVLeftSend, requestVRightSend, requestVLeftRecv, requestVRightRecv,
        requestVTopSend, requestVBottomRecv,
        requestUBottomSend, requestUTopSend, requestUBottomRecv, requestUTopRecv,
        requestURightSend, requestULeftRecv
    };
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

    //Put received values into grid
    if (leftRankNo != MPI_PROC_NULL) {
        for (int j = u.beginJ(); j < u.endJ(); j++) {
            u(-1, j) = leftURecvBuffer[j + 1];
        }
        for (int j = v.beginJ(); j < v.endJ(); j++) {
            v(-1, j) = leftVRecvBuffer[j + 1];
        }
    }
    if (rightRankNo != MPI_PROC_NULL) {
        for (int j = v.beginJ(); j < v.endJ(); j++) {
            v(v.endI() - 1, j) = rightVRecvBuffer[j + 1];
        }
    }
    if (bottomRankNo != MPI_PROC_NULL) {
        for (int i = u.beginI(); i < u.endI(); i++) {
            u(i, -1) = bottomURecvBuffer[i + 1];
        }
        for (int i = v.beginI(); i < v.endI(); i++) {
            v(i, -1) = bottomVRecvBuffer[i + 1];
        }
    }
    if (topRankNo != MPI_PROC_NULL) {
        for (int i = u.beginI(); i < u.endI(); i++) {
            u(i, u.endJ() - 1) = topURecvBuffer[i + 1];
        }
    }

    //ToDo (Bene): Fix errors and refactor this huge duplicating piece of sh*t so it becomes less ugly
}
