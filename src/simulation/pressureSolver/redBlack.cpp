#include "redBlack.h"
#include "macros.h"
#include <limits>
#include <iostream>

RedBlack::RedBlack(const std::shared_ptr<StaggeredGrid> &grid,
                   const double epsilon,
                   const int maximumNumberOfIterations,
                   const double omega,
                   const std::shared_ptr<Partitioning> &partitioning)
    : partitioning_(partitioning), grid_(grid), epsilon_(epsilon), maxNumberOfIterations_(maximumNumberOfIterations), omega_(omega) {}


// ToDo: fix solvererror not being local
void RedBlack::solve() {
    int it = 0;

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1 / dx2;
    const double invDy2 = 1 / dy2;
    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

    const int residualStartThreshold = static_cast<int>(0.7 * lastIterationCount_);

    localPressureError_ = std::numeric_limits<double>::max();
    globalPressureError_ = std::numeric_limits<double>::max();

    const int beginI = p.beginI() + 1;
    const int endI = p.endI() - 1;

    const int beginJ = p.beginJ() + 1;
    const int endJ = p.endJ() - 1;

    // ToDo: Do we need to reset error after every iteration?
    while (it < maxNumberOfIterations_) {
        // rb = 0 is red, rb=1 is black pass
        for (int rb = 0; rb < 2; ++rb) {
            for (int j = beginJ; j < endJ; ++j) {
                // ((beginI + j + rb) & 1) calculates offset 0 or 1
                for (int i = beginI + ((beginI + j + rb) & 1); i < endI; i += 2) {
                    const double pDxx = (p(i - 1, j) + p(i + 1, j)) * invDx2;
                    const double pDyy = (p(i, j - 1) + p(i, j + 1)) * invDy2;
                    const double pNew = (1 - omega_) * p(i, j) + omega_ * scalingFactor * (pDxx + pDyy - rhs(i, j));
                    localPressureError_ += (p(i, j) - pNew) * (p(i, j) - pNew);
                    p(i, j) = pNew;
                }
            }
            updatePressureBoundaries();
        }

        // ToDo: Global is always inf so this is wrong
        MPI_Allreduce(&localPressureError_, &globalPressureError_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // DEBUG(std::cout << "Local pressure Error: " << localPressureError_ << ", global pressure error: " << globalPressureError_ << std::endl;)

        if (it > residualStartThreshold && globalPressureError_ < epsilon_ * epsilon_) {
            break;
        }
        it++;
    }
    lastIterationCount_ = it;
}

void RedBlack::updatePressureBoundaries() {
    DataField &p = grid_->p();
    requests_.clear();

    const int numberX = p.endI() - p.beginI();
    const int numberY = p.endJ() - p.beginJ();

    // ToDo: This should probably be global object on the class so we only initialize once
    MPI_Datatype mpiColumn;
    MPI_Type_vector(numberY, 1, numberX, MPI_DOUBLE, &mpiColumn);
    MPI_Type_commit(&mpiColumn);

    for (const Direction direction : directions()) {
        if (!partitioning_->ownPartitionContainsBoundary(direction)) {
            const int neighborRank = partitioning_->neighborRankNo(direction);
            MPI_Request sendReq{}, recvReq{};
            // ToDo: I'm not sure about the indices here! (e.g. should right/top be -1 and 0 or -2 and -1)
            switch (direction) {
            case Direction::Left: {
                MPI_Isend(&p(p.beginI() + 1, p.beginJ()), 1, mpiColumn, neighborRank, P_TAG, MPI_COMM_WORLD, &sendReq);
                MPI_Irecv(&p(p.beginI(), p.beginJ()), 1, mpiColumn, neighborRank, P_TAG, MPI_COMM_WORLD, &recvReq);
                break;
            }
            case Direction::Right: {
                MPI_Isend(&p(p.endI() - 2, p.beginJ()), 1, mpiColumn, neighborRank, P_TAG, MPI_COMM_WORLD, &sendReq);
                MPI_Irecv(&p(p.endI() - 1, p.beginJ()), 1, mpiColumn, neighborRank, P_TAG, MPI_COMM_WORLD, &recvReq);
                break;
            }
            case Direction::Bottom: {
                MPI_Isend(&p(p.beginI(), p.beginJ() + 1), numberX, MPI_DOUBLE, neighborRank, P_TAG, MPI_COMM_WORLD, &sendReq);
                MPI_Irecv(&p(p.beginI(), p.beginJ()), numberX, MPI_DOUBLE, neighborRank, P_TAG, MPI_COMM_WORLD, &recvReq);
                break;
            }
            case Direction::Top: {
                MPI_Isend(&p(p.beginI(), p.endJ() - 2), numberX, MPI_DOUBLE, neighborRank, P_TAG, MPI_COMM_WORLD, &sendReq);
                MPI_Irecv(&p(p.beginI(), p.endJ() - 1), numberX, MPI_DOUBLE, neighborRank, P_TAG, MPI_COMM_WORLD, &recvReq);
                break;
            }
            }
            requests_.push_back(sendReq);
            requests_.push_back(recvReq);
        } else {
            setBoundaryValues(direction);
        }
    }

    if (!requests_.empty()) {
        MPI_Waitall(requests_.size(), requests_.data(), MPI_STATUSES_IGNORE);
    }
}

// I did this with a template set method before, but it's probably faster this way
void RedBlack::setBoundaryValues(const Direction direction) {
    DataField &p = grid_->p();

    switch (direction) {
    case Direction::Left: {
        for (int j = p.beginJ(); j < p.endJ(); j++) {
            p(p.beginI(), j) = p(p.beginI() + 1, j);
        }
        break;
    }
    case Direction::Right: {
        for (int j = p.beginJ(); j < p.endJ(); j++) {
            p(p.endI() - 1, j) = p(p.endI() - 2, j);
        }
        break;
    }
    case Direction::Bottom: {
        for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
            p(i, p.beginJ()) = p(i, p.beginJ() +1 );
        }
        break;
    }
    case Direction::Top: {
        for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
            p(i, p.endJ() - 1) = p(i, p.endJ() -2);
        }
        break;
    }
    }
}
