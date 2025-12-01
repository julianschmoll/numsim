#include "redBlack.h"
#include "macros.h"
#include <limits>
#include <iostream>

RedBlack::RedBlack(const std::shared_ptr<StaggeredGrid> &grid,
                   const double epsilon,
                   const int maximumNumberOfIterations,
                   const double omega,
                   const std::shared_ptr<Partitioning> &partitioning)
    : partitioning_(partitioning), grid_(grid), epsilon_(epsilon), maxNumberOfIterations_(maximumNumberOfIterations), omega_(omega) {
    const DataField &p = grid_->p();
    const int numberX = p.endI() - p.beginI();
    const int numberY = p.endJ() - p.beginJ();

    MPI_Type_vector(numberY, 1, numberX, MPI_DOUBLE, &mpiColumn_);
    MPI_Type_commit(&mpiColumn_);
}

RedBlack::~RedBlack() {
    MPI_Type_free(&mpiColumn_);
}


// ToDo: I think the residual is still fucked up a little
void RedBlack::solve() {
    int it = 0;

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1 / dx2;
    const double invDy2 = 1 / dy2;
    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

    globalPressureError_ = std::numeric_limits<double>::max();

    const int beginI = p.beginI() + 1;
    const int endI = p.endI() - 1;

    const int beginJ = p.beginJ() + 1;
    const int endJ = p.endJ() - 1;

    while (it < maxNumberOfIterations_) {
        localPressureError_ = 0.0;
        // rb = 0 is red, rb=1 is black pass
        // #pragma omp parallel for reduction(+:localPressureError_)
        for (int rb = 0; rb < 2; ++rb) {
            // #pragma omp simd reduction(+:localPressureError_)
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

        MPI_Allreduce(&localPressureError_, &globalPressureError_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // DEBUG(std::cout << "Local pressure Error: " << localPressureError_ << ", global pressure error: " << globalPressureError_ << std::endl;)
        if (globalPressureError_ < epsilon_ * epsilon_) {
            break;
        }
        it++;
    }
}

void RedBlack::updatePressureBoundaries() {
    DataField &p = grid_->p();
    requests_.clear();
    requests_.reserve(8);

    const int numberX = p.endI() - p.beginI();

    for (const Direction direction : directions()) {
        if (partitioning_->ownPartitionContainsBoundary(direction)) {
            setBoundaryValues(direction);
        } else {
            const int neighborRank = partitioning_->neighborRankNo(direction);
            MPI_Request sendReq{}, recvReq{};
            switch (direction) {
            case Direction::Left: {
                MPI_Isend(&p(p.beginI() + 1, p.beginJ()), 1, mpiColumn_, neighborRank, P_TAG, MPI_COMM_WORLD, &sendReq);
                MPI_Irecv(&p(p.beginI(), p.beginJ()), 1, mpiColumn_, neighborRank, P_TAG, MPI_COMM_WORLD, &recvReq);
                break;
            }
            case Direction::Right: {
                MPI_Isend(&p(p.endI() - 2, p.beginJ()), 1, mpiColumn_, neighborRank, P_TAG, MPI_COMM_WORLD, &sendReq);
                MPI_Irecv(&p(p.endI() - 1, p.beginJ()), 1, mpiColumn_, neighborRank, P_TAG, MPI_COMM_WORLD, &recvReq);
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
            p(i, p.beginJ()) = p(i, p.beginJ() + 1);
        }
        break;
    }
    case Direction::Top: {
        for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
            p(i, p.endJ() - 1) = p(i, p.endJ() - 2);
        }
        break;
    }
    }
}
