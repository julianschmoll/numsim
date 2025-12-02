#pragma once

#include "../grid/dataField.h"
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

class Partitioning {
    std::array<int, 2> nCellsLocal_{};
    std::array<int, 2> nCellsGlobal_{};

    /// Number of ranks in each direction.
    std::array<int, 2> partitions_{};

    std::array<int, 2> rankCoordinates_;
    
    /// MPI communicator with cartesian coordinate information attached.
    MPI_Comm cartComm_ = nullptr;

    /// Node offsets
    int offsetX_ = 0;
    int offsetY_ = 0;

    Rank ownRank_ = 0;


public:
    Partitioning() = delete;

    //! compute partitioning, set internal variables
    explicit Partitioning(std::array<int, 2> nCellsGlobal);

    //! get the local number of cells in the own subdomain
    std::array<int, 2> nCellsLocal() const;

    //! get the global number of cells in the whole computational domain
    //! used in OutputWriterParaviewParallel
    std::array<int, 2> nCellsGlobal() const;

    //! get the own MPI rank no
    //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
    Rank ownRank() const;

    //! number of MPI ranks
    int nRanks() const;

    //! get the offset values for counting local nodes in x and y direction.
    //! (i_local,j_local) + nodeOffset = (i_global,j_global)
    //! used in OutputWriterParaviewParallel
    std::array<int, 2> nodeOffset() const;
    
    std::array<int, 2> getCurrentRankCoords() const;
    
    void exchange(const std::vector<DataField *> &fields) const;
    
    void barrier() const;
    
    bool onPrimaryRank() const;
    
    void printPartitioningInfo() const;

    //! if the own partition has part of the bottom boundary of the whole domain
    template<Direction direction>
    bool ownContainsBoundary() const {
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
        
    //! get the rank no of the left neighbouring rank
    template <Direction direction>
    Rank neighborRank() const {
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

    template <Direction direction>
    MPI_Request sendBorder(DataField &field) const {
        MPI_Request sendReq{};

        int i = 0, j = 0;
        int count = field.cols();
        MPI_Datatype type = MPI_DOUBLE;
        MPI_Datatype coltype{};

        // For some reason MPI does not let us store this type in DataField
        MPI_Type_vector(field.rows(), 1, field.cols(), type, &coltype);
        MPI_Type_commit(&coltype);

        if constexpr (direction == Direction::Left) {
            i = field.beginI() + 1;
            j = field.beginJ();
            type = coltype;
            count = 1;
        } else if constexpr (direction == Direction::Right) {
            i = field.endI() - 2;
            j = field.beginJ();
            type = coltype;
            count = 1;
        } else if constexpr (direction == Direction::Bottom) {
            i = field.beginI();
            j = field.beginJ() + 1;
        } else if constexpr (direction == Direction::Top) {
            i = field.beginI();
            j = field.endJ() - 2;
        } else {
            assert(false);
        }

        MPI_Isend(&field(i, j), count, type, neighborRank<direction>(), field.getID(), cartComm_, &sendReq);

        MPI_Type_free(&coltype);

        return sendReq;
    }

    template <Direction direction>
    MPI_Request recvBorder(DataField &field) const {
        MPI_Request recvReq{};

        int i = 0, j = 0;
        int count = field.cols();
        MPI_Datatype type = MPI_DOUBLE;
        MPI_Datatype coltype{};

        // For some reason MPI does not let us store this type in DataField
        MPI_Type_vector(field.rows(), 1, field.cols(), type, &coltype);
        MPI_Type_commit(&coltype);

        if constexpr (direction == Direction::Left) {
            i = field.beginI();
            j = field.beginJ();
            type = coltype;
            count = 1;
        } else if constexpr (direction == Direction::Right) {
            i = field.endI() - 1;
            j = field.beginJ();
            type = coltype;
            count = 1;
        } else if constexpr (direction == Direction::Bottom) {
            i = field.beginI();
            j = field.beginJ();
        } else if constexpr (direction == Direction::Top) {
            i = field.beginI();
            j = field.endJ() - 1;
        } else {
            assert(false);
        }

        MPI_Irecv(&field(i, j), count, type, neighborRank<direction>(), field.getID(), cartComm_, &recvReq);

        MPI_Type_free(&coltype);

        return recvReq;
    }
};
