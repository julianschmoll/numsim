#include "partitioning.h"
#include <mpi.h>
#include <iostream>

int Partitioning::ownRankNo() const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

void Partitioning::initialize(std::array<int, 2> nCellsGlobal) {
    ownRankNo_ = ownRankNo();
    int numberRanks = nRanks();
    nCellsGlobal_ = nCellsGlobal;

    constexpr int dimensions = 2;
    int part[2] = {0,0};

    // partitions stores the desired dimensions of the grid (e.g. 2x3 for 6 cells)
    MPI_Dims_create(numberRanks, dimensions, part);

    constexpr int periodic[dimensions] = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, part, periodic, false, &cartComm_);

    partitions = {part[0], part[1]};

    MPI_Cart_coords(cartComm_, ownRankNo_, dimensions, rankCoordinates_.data());
    const int extraCellsX = rankCoordinates_[0] < nCellsGlobal_[0] % partitions[0] ? 1 : 0;
    const int extraCellsY = rankCoordinates_[1] < nCellsGlobal_[1] % partitions[1] ? 1 : 0;

    nCellsLocal_[0] = nCellsGlobal[0] / partitions[0] + extraCellsX;
    nCellsLocal_[1] = nCellsGlobal[1] / partitions[1] + extraCellsY;

    // offset vector of the partition
    const int partWidth = nCellsGlobal_[0] / partitions[0];
    const int partHeight = nCellsGlobal_[1] / partitions[1];

    offsetX_ = partWidth * rankCoordinates_[0] + std::min(rankCoordinates_[0], nCellsGlobal[0] % partitions[0]);
    offsetY_ = partHeight * rankCoordinates_[1] + std::min(rankCoordinates_[1], nCellsGlobal[1] % partitions[1]);
}

int Partitioning::nRanks() const {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

std::array<int, 2> Partitioning::nCellsLocal() const {
    return nCellsLocal_;
}

std::array<int, 2> Partitioning::nCellsGlobal() const {
    return nCellsGlobal_;
}

std::array<int, 2> Partitioning::getCurrentRankCoords() const {
    std::array<int, dimensions_> coordinates;
    MPI_Cart_coords(cartComm_, ownRankNo_, dimensions_, coordinates.data());
    return coordinates;
}

bool Partitioning::ownPartitionContainsBottomBoundary() const {
    return rankCoordinates_[1] == 0;
}

bool Partitioning::ownPartitionContainsTopBoundary() const {
    /*int coordinates[dimensions_] = {0, 0};
    *MPI_Cart_coords(cartComm_, ownRankNo_, dimensions_, coordinates);
    return*/
    return rankCoordinates_[1] == partitions[1] - 1;
}

bool Partitioning::ownPartitionContainsLeftBoundary() const {
    return rankCoordinates_[0] == 0;
}

bool Partitioning::ownPartitionContainsRightBoundary() const {
    return rankCoordinates_[0] == partitions[0] - 1;
}

int Partitioning::leftNeighbourRankNo() const {
    int destRank;
    int srcRank = ownRankNo_;
    MPI_Cart_shift(cartComm_, 0, -1, &srcRank, &destRank);
    return destRank;
    // const int rowStart = ownRankNo_ - (ownRankNo_ % partitions[0]);
    // return ownRankNo_ - 1 < rowStart ? MPI_PROC_NULL : ownRankNo_ - 1;
}

int Partitioning::rightNeighbourRankNo() const {
    int destRank;
    int srcRank = ownRankNo_;
    MPI_Cart_shift(cartComm_, 0, 1, &srcRank, &destRank);
    return destRank;
    // const int rowEnd = ownRankNo_ - (ownRankNo_ % partitions[0]) + partitions[0] - 1;
    // return ownRankNo_ + 1 > rowEnd ? MPI_PROC_NULL : ownRankNo_ + 1;
}

int Partitioning::topNeighbourRankNo() const {
    int destRank;
    int srcRank = ownRankNo_;
    MPI_Cart_shift(cartComm_, 1, 1, &srcRank, &destRank);
    return destRank;
    //return ownRankNo_ + partitions[0] >= nRanks() ? MPI_PROC_NULL : ownRankNo_ + partitions[0];
}

int Partitioning::bottomNeighbourRankNo() const {
    int destRank;
    int srcRank = ownRankNo_;
    MPI_Cart_shift(cartComm_, 1, -1, &srcRank, &destRank);
    return destRank;
    //return ownRankNo_ - partitions[0] < 0 ? MPI_PROC_NULL : ownRankNo_ - partitions[0];
}

std::array<int, 2> Partitioning::nodeOffset() const {
    return {offsetX_, offsetY_};
}
