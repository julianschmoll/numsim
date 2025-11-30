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
    const int numberRanks = nRanks();
    nCellsGlobal_ = nCellsGlobal;

    constexpr int dimensions = 2;
    int part[2] = {0, 0};

    // partitions stores the desired dimensions of the grid (e.g. 2x3 for 6 cells)
    MPI_Dims_create(numberRanks, dimensions, part);

    constexpr int periodic[dimensions] = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, part, periodic, false, &cartComm_);

    partitions_ = {part[0], part[1]};

    MPI_Cart_coords(cartComm_, ownRankNo_, dimensions, rankCoordinates_.data());
    const int extraCellsX = rankCoordinates_[0] < nCellsGlobal_[0] % partitions_[0] ? 1 : 0;
    const int extraCellsY = rankCoordinates_[1] < nCellsGlobal_[1] % partitions_[1] ? 1 : 0;

    nCellsLocal_[0] = nCellsGlobal[0] / partitions_[0] + extraCellsX;
    nCellsLocal_[1] = nCellsGlobal[1] / partitions_[1] + extraCellsY;

    // offset vector of the partition
    const int partWidth = nCellsGlobal_[0] / partitions_[0];
    const int partHeight = nCellsGlobal_[1] / partitions_[1];

    offsetX_ = partWidth * rankCoordinates_[0] + std::min(rankCoordinates_[0], nCellsGlobal[0] % partitions_[0]);
    offsetY_ = partHeight * rankCoordinates_[1] + std::min(rankCoordinates_[1], nCellsGlobal[1] % partitions_[1]);

    std::cout << "I'm partition " << ownRankNo_ << " with cells "
    << offsetX_ << "-" << offsetX_ + nCellsLocal_[0] -1 << ", " << offsetY_ << "-" << offsetY_ + nCellsLocal_[1] -1
    << ". Left neighbor: " << neighborRankNo(Direction::Left)
    << " right neighbor: " << neighborRankNo(Direction::Right)
    << " top neighbor: " << neighborRankNo(Direction::Top)
    << " bottom neighbor: " << neighborRankNo(Direction::Bottom) << std::endl;
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

bool Partitioning::ownPartitionContainsBoundary(const Direction direction) const {
    switch (direction) {
    case Direction::Left:
        return rankCoordinates_[0] == 0;
    case Direction::Right:
        return rankCoordinates_[0] == partitions_[0] - 1;
    case Direction::Bottom:
        return rankCoordinates_[1] == 0;
    case Direction::Top:
        return rankCoordinates_[1] == partitions_[1] - 1;
    }
    return false;
}

int Partitioning::neighborRankNo(const Direction direction) const {
    int destRank;
    int srcRank = ownRankNo_;
    switch (direction) {
    case Direction::Left:
        MPI_Cart_shift(cartComm_, 0, -1, &srcRank, &destRank);
        break;
    case Direction::Right:
        MPI_Cart_shift(cartComm_, 0, 1, &srcRank, &destRank);
        break;
    case Direction::Top:
        MPI_Cart_shift(cartComm_, 1, 1, &srcRank, &destRank);
        break;
    case Direction::Bottom:
        MPI_Cart_shift(cartComm_, 1, -1, &srcRank, &destRank);
        break;
    }
    return destRank;
}

std::array<int, 2> Partitioning::nodeOffset() const {
    return {offsetX_, offsetY_};
}
