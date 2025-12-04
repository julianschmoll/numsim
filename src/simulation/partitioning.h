#pragma once

#include "grid/dataField.h"
#include <array>
#include <cassert>
#include <mpi.h>

using Rank = int;

// This avoids magic numbers in shift and can be used in direction as well
enum class Direction { Left, Right, Top, Bottom };

// This can be used to iterate over every direction and avoid redundancy in code :)
constexpr std::array<Direction, 4> directions() {
    return {Direction::Left, Direction::Right, Direction::Bottom, Direction::Top};
}

/**
 * @class Partitioning
 * @brief Handles partitioning the simulation domain into subdomains for multithreaded simulation.
 */
class Partitioning {
    /// Dimensions of the local domain
    std::array<int, 2> nCellsLocal_{};

    /// Dimensions of the global domain
    std::array<int, 2> nCellsGlobal_{};

    /// Number of ranks in each direction.
    std::array<int, 2> partitions_{};

    std::array<int, 2> rankCoordinates_;

    /// MPI communicator with cartesian coordinate information attached.
    MPI_Comm cartComm_ = nullptr;

    /// Node offset in x direction
    int offsetX_ = 0;

    /// Node offset in y direction
    int offsetY_ = 0;

    /// MPI Rank number of the current Rank
    Rank ownRank_ = 0;

public:
    Partitioning() = delete;

    /**
     * Handles partitioning the grid and computes all domain dimensions and offsets.
     *
     * @param nCellsGlobal: Size of the global domain
     */
    explicit Partitioning(std::array<int, 2> nCellsGlobal);

    /**
     * Gets the local number of cells in the own subdomain.
     *
     * @return Dimension of partition.
     */
    std::array<int, 2> nCellsLocal() const;

    /**
     * Gets the global number of cells in the whole computational domain
     *
     * @return Dimension of the global domain.
     */
    std::array<int, 2> nCellsGlobal() const;

    /**
     * Gets the own MPI Rank Number.
     *
     * @return Rank Number.
     */
    Rank ownRank() const;

    /**
     * Gets the number of MPI ranks.
     *
     * @return Number of ranks.
     */
    int nRanks() const;
    /**
     * Gets the offset values for counting local nodes in x and y direction.
     *
     * @return Offset to bottom left corner of the partition domain.
     */
    std::array<int, 2> nodeOffset() const;

    /**
     * Gets the coordinates of the own rank.
     *
     * @return MPI Cartesian coordinates of the own rank.
     */
    std::array<int, 2> getCurrentRankCoords() const;

    /**
     * Performs a blocking exchange of borders of the data fields.
     *
     * @param fields Vector holding data fields to exchange.
     */
    void exchange(const std::vector<DataField *> &fields);

    /**
     * Performs a blocking exchange of the borders of one data field.
     *
     * @param field Data field to exchange data.
     */
    void exchange(DataField &field);

    /**
     * Performs a non-blocking exchange of borders of the data fields.
     *
     * @param fields Vector holding data fields to exchange.
     */
    void nonBlockingExchange(const std::vector<DataField *> &fields);

    /**
     * Performs a non-blocking exchange of the borders of one data field.
     *
     * @param field Data field to exchange data.
     */
    void nonBlockingExchange(DataField &field);

    /**
     * MPI Barrier.
     *
     * Waits until all processes have reached the barrier.
     */
    void barrier() const;

    /**
     * Check whether the current rank is the primary rank.
     *
     * @return If the current rank is the primary rank.
     */
    bool onPrimaryRank() const;

    /**
     * Prints info about how the domain is partitioned.
     */
    void printPartitioningInfo() const;

    /**
     * Waits until all non-blocking exchanges are finished.
     */
    void waitForAllMPIRequests();

    /**
     * Uses MPI Reduce to calculate a sum from multiple subdomains.
     *
     * @param localSum Sum on one partition.
     * @return Summed up local sums.
     */
    double collectSum(double localSum) const;

    /**
     * Checks if the partition contains a certain boundary.
     *
     * @tparam direction Direction to check.
     * @return if the partition contains that boundary.
     */
    template <Direction direction> bool ownContainsBoundary() const {
        if constexpr (direction == Direction::Left) {
            return rankCoordinates_[0] == 0;
        }
        if constexpr (direction == Direction::Right) {
            return rankCoordinates_[0] == partitions_[0] - 1;
        }
        if constexpr (direction == Direction::Bottom) {
            return rankCoordinates_[1] == 0;
        }
        if constexpr (direction == Direction::Top) {
            return rankCoordinates_[1] == partitions_[1] - 1;
        }
        return false;
    }

    /**
     * Gets the rank number of neighboring rank.
     *
     * @tparam direction Direction to get rank number from.
     * @return Rank number.
     */
    template <Direction direction> Rank neighborRank() const {
        Rank destRank{};
        Rank srcRank = ownRank_;
        if constexpr (direction == Direction::Left) {
            MPI_Cart_shift(cartComm_, 0, -1, &srcRank, &destRank);
        }
        if constexpr (direction == Direction::Right) {
            MPI_Cart_shift(cartComm_, 0, 1, &srcRank, &destRank);
        }
        if constexpr (direction == Direction::Top) {
            MPI_Cart_shift(cartComm_, 1, 1, &srcRank, &destRank);
        }
        if constexpr (direction == Direction::Bottom) {
            MPI_Cart_shift(cartComm_, 1, -1, &srcRank, &destRank);
        }
        return destRank;
    }

    /**
     * Performs MPI Isend to exchange border values on a specified border.
     *
     * @tparam direction Direction of the border.
     * @param field Data field to send border data from.
     * @return Send request handle.
     */
    template <Direction direction> MPI_Request sendBorder(DataField &field) const {
        MPI_Request sendReq{};

        int i = 0, j = 0;
        int count = field.cols() - 2;
        auto type = MPI_DOUBLE;

        if constexpr (direction == Direction::Left) {
            i = field.beginI() + 1;
            j = field.beginJ();
            type = field.mpiColType();
            count = 1;
        } else if constexpr (direction == Direction::Right) {
            i = field.endI() - 2;
            j = field.beginJ();
            type = field.mpiColType();
            count = 1;
        } else if constexpr (direction == Direction::Bottom) {
            i = field.beginI() + 1;
            j = field.beginJ() + 1;
        } else if constexpr (direction == Direction::Top) {
            i = field.beginI() + 1;
            j = field.endJ() - 2;
        } else {
            assert(false);
        }

        MPI_Isend(&field(i, j), count, type, neighborRank<direction>(), field.getID(), cartComm_, &sendReq);

        return sendReq;
    }

    /**
     * Performs MPI IReceive to exchange border values on a specified border.
     *
     * @tparam direction Direction of the border.
     * @param field Data field to receive border data from.
     * @return Send request handle.
     */
    template <Direction direction> MPI_Request recvBorder(DataField &field) const {
        MPI_Request recvReq{};

        int i = 0, j = 0;
        int count = field.cols() - 2;
        auto type = MPI_DOUBLE;

        if constexpr (direction == Direction::Left) {
            i = field.beginI();
            j = field.beginJ();
            type = field.mpiColType();
            count = 1;
        } else if constexpr (direction == Direction::Right) {
            i = field.endI() - 1;
            j = field.beginJ();
            type = field.mpiColType();
            count = 1;
        } else if constexpr (direction == Direction::Bottom) {
            i = field.beginI() + 1;
            j = field.beginJ();
        } else if constexpr (direction == Direction::Top) {
            i = field.beginI() + 1;
            j = field.endJ() - 1;
        } else {
            assert(false);
        }

        MPI_Irecv(&field(i, j), count, type, neighborRank<direction>(), field.getID(), cartComm_, &recvReq);

        return recvReq;
    }

private:
    /// Holds all MPI requests of non-blocking exchanges.
    mutable std::vector<MPI_Request> requests_;
};
