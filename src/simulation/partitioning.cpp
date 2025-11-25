#include "partitioning.h"
#include <mpi.h>

int Partitioning::ownRankNo() const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

void Partitioning::initialize(std::array<int, 2> nCellsGlobal) {
    int rank = ownRankNo();
    int nRanks = nRanks();
    nCellsGlobal_ = nCellsGlobal;

    constexpr int dimensions = 2;

    // partitions stores the desired dimensions of the grid (e.g. 2x3 for 6 cells)
    MPI_Dims_create(nRanks, dimensions, partitions.data());

    constexpr int periodic[dimensions] = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, partitions.data(), periodic, false, &cartComm_);

    MPI_Cart_coords(cartComm_, rank, dimensions, rankCoordinates_.data());
    const int extraCellsX = rankCoordinates_[0] < nCellsGlobal_[0] % partitions[0] ? 1 : 0;
    const int extraCellsY = rankCoordinates_[1] < nCellsGlobal_[1] % partitions[1] ? 1 : 0;

    nCellsLocal_[0] = nCellsGlobal[0] / partitions[0] + extraCellsX;
    nCellsLocal_[1] = nCellsGlobal[1] / partitions[1] + extraCellsY;

    // offset vector of the partition
    const int partWidth = nCellsGlobal_[0] / partitions[0]; //ToDo: Replace later with nCellsLocal_[0]
    const int partHeight = nCellsGlobal_[1] / partitions[1]; //ToDo: Replace later with nCellsLocal_[1]

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
    MPI_Cart_coords(cartComm_, ownRankNo(), dimensions_, coordinates.data());
    return coordinates;
}

bool Partitioning::ownPartitionContainsBottomBoundary() const {
    return rankCoordinates_[1] == 0;
}

bool Partitioning::ownPartitionContainsTopBoundary() const {
    /*int coordinates[dimensions_] = {0, 0};
    *MPI_Cart_coords(cartComm_, ownRankNo(), dimensions_, coordinates);
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
    /*int rank_source, rank_dest;
    // direction=0 (x-axis), disp=-1 (shift to the left)
    MPI_Cart_shift(comm_cart_, 0, -1, &rank_source, &rank_dest);
    // rank_dest is the rank to send data TO
    return rank_dest;
     *
     */

    const int rank = ownRankNo();
    const int rowStart = rank - (rank % partitions[0]);
    return rank - 1 < rowStart ? MPI_PROC_NULL : rank - 1;
}

int Partitioning::rightNeighbourRankNo() const {
    const int rank = ownRankNo();
    const int rowEnd = rank - (rank % partitions[0]) + partitions[0] - 1;
    return rank + 1 > rowEnd ? MPI_PROC_NULL : rank + 1;
}

int Partitioning::topNeighbourRankNo() const {
    const int rank = ownRankNo();
    return rank + partitions[0] >= nRanks() ? MPI_PROC_NULL : rank + partitions[0];
}

int Partitioning::bottomNeighbourRankNo() const {
    const int rank = ownRankNo();
    return rank - partitions[0] < 0 ? MPI_PROC_NULL : rank - partitions[0];
}

std::array<int, 2> Partitioning::nodeOffset() const {
    return {offsetX_, offsetY_};
}
