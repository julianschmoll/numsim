#include "partitioning.h"
#include <iostream>
#include <mpi.h>

Partitioning::Partitioning(std::array<int, 2> nCellsGlobal) : nCellsGlobal_(nCellsGlobal), rankCoordinates_({0, 0}) {
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);

    const int numberRanks = nRanks();

    constexpr int dimensions = 2;
    int part[2] = {0, 0}; // std::array::data() doesn't work for some reason.

    // partitions stores the desired dimensions of the grid (e.g. 2x3 for 6 cells)
    MPI_Dims_create(numberRanks, dimensions, part);

    constexpr int periodic[dimensions] = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, part, periodic, false, &cartComm_);

    partitions_ = {part[0], part[1]};

    MPI_Cart_coords(cartComm_, ownRank_, dimensions, rankCoordinates_.data());
    const int extraCellsX = rankCoordinates_[0] < nCellsGlobal_[0] % partitions_[0] ? 1 : 0;
    const int extraCellsY = rankCoordinates_[1] < nCellsGlobal_[1] % partitions_[1] ? 1 : 0;

    nCellsLocal_[0] = nCellsGlobal[0] / partitions_[0] + extraCellsX;
    nCellsLocal_[1] = nCellsGlobal[1] / partitions_[1] + extraCellsY;

    // offset vector of the partition
    const int partWidth = nCellsGlobal_[0] / partitions_[0];
    const int partHeight = nCellsGlobal_[1] / partitions_[1];

    offsetX_ = partWidth * rankCoordinates_[0] + std::min(rankCoordinates_[0], nCellsGlobal[0] % partitions_[0]);
    offsetY_ = partHeight * rankCoordinates_[1] + std::min(rankCoordinates_[1], nCellsGlobal[1] % partitions_[1]);
}

int Partitioning::nRanks() const {
    int size{};
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
    std::array<int, 2> coordinates{};
    MPI_Cart_coords(cartComm_, ownRank_, 2, coordinates.data());
    return coordinates;
}

void Partitioning::barrier() const {
    MPI_Barrier(cartComm_);
}

bool Partitioning::onPrimaryRank() const {
    return ownRank_ == 0;
}

std::array<int, 2> Partitioning::nodeOffset() const {
    return {offsetX_, offsetY_};
}

int Partitioning::ownRank() const {
    return ownRank_;
}

void Partitioning::printPartitioningInfo() const {
    std::cout << "Partitioning: rank " << ownRank_ << ": (" << rankCoordinates_[0] << ", " << rankCoordinates_[1] << ")\n";
    if (onPrimaryRank()) {
        std::cout << " -- nCellsGlobal_ = [" << nCellsLocal_[0] << ", " << nCellsLocal_[1] << "]\n";
    }
    std::cout << " -- nCellsLocal_ = [" << nCellsLocal_[0] << ", " << nCellsLocal_[1] << "]\n";
    std::cout << " -- Neighbours: top = " << neighborRank<Direction::Top>() << ", "
              << "right = " << neighborRank<Direction::Right>() << ", "
              << "bottom = " << neighborRank<Direction::Bottom>() << ", "
              << "left = " << neighborRank<Direction::Left>();
    std::cout << std::endl;
}

void Partitioning::exchange(const std::vector<DataField *> &fields) const {
    std::vector<MPI_Request> requests;

    for (DataField *field : fields) {
        if (!ownContainsBoundary<Direction::Top>()) {
            requests.emplace_back(sendBorder<Direction::Top>(*field));
            requests.emplace_back(recvBorder<Direction::Top>(*field));
        }
        if (!ownContainsBoundary<Direction::Bottom>()) {
            requests.emplace_back(sendBorder<Direction::Bottom>(*field));
            requests.emplace_back(recvBorder<Direction::Bottom>(*field));
        }
        if (!ownContainsBoundary<Direction::Left>()) {
            requests.emplace_back(sendBorder<Direction::Left>(*field));
            requests.emplace_back(recvBorder<Direction::Left>(*field));
        }
        if (!ownContainsBoundary<Direction::Right>()) {
            requests.emplace_back(sendBorder<Direction::Right>(*field));
            requests.emplace_back(recvBorder<Direction::Right>(*field));
        }
    }
    if (!requests.empty()) {
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUS_IGNORE);
    }
}