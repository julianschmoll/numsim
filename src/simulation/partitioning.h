#pragma once

#include "../grid/dataField.h"
#include <array>
#include <cassert>
#include <iostream>
#include <mpi.h>
#include <string>

// This avoids magic numbers in shift and can be used in direction as well
enum class Direction { Left, Right, Top, Bottom };

// This can be used to iterate over every direction and avoid redundancy in code :)
constexpr std::array<Direction, 4> directions() {
    return {Direction::Left, Direction::Right, Direction::Bottom, Direction::Top};
}

class Partitioning {
    /// Dimensions of the problem: 2D grid with x,y directions.
    static constexpr int dimensions_ = 2; // TODO: move

    std::array<int, 2> nCellsLocal_ = {};
    std::array<int, 2> nCellsGlobal_ = {};

    /// Number of ranks in each direction.
    std::array<int, 2> partitions_ = {};

    /// MPI communicator with cartesian coordinate information attached.
    MPI_Comm cartComm_ = nullptr;

    /// Node offsets
    int offsetX_ = 0;
    int offsetY_ = 0;

    int ownRankNo_ = 0;

    std::array<int, dimensions_> rankCoordinates_ = {};

public:
    // TODO: tabs

    //! compute partitioning, set internal variables
    void initialize(std::array<int, 2> nCellsGlobal);

    //! get the local number of cells in the own subdomain
    std::array<int, 2> nCellsLocal() const;

    //! get the global number of cells in the whole computational domain
    //! used in OutputWriterParaviewParallel
    std::array<int, 2> nCellsGlobal() const;

    //! get the own MPI rank no
    //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
    int ownRankNo() const;

    //! number of MPI ranks
    int nRanks() const;

    //! if the own partition has part of the bottom boundary of the whole domain
    bool ownPartitionContainsBoundary(Direction direction) const;

    //! get the rank no of the left neighbouring rank
    int neighborRankNo(Direction direction) const;

    //! get the offset values for counting local nodes in x and y direction.
    //! (i_local,j_local) + nodeOffset = (i_global,j_global)
    //! used in OutputWriterParaviewParallel
    std::array<int, 2> nodeOffset() const;

    std::array<int, 2> getCurrentRankCoords() const;

    void exchange(const std::vector<DataField *> &fields) const;

    inline std::string dirToStr(Direction dir) const {
        if (dir == Direction::Left)
            return "Left";
        else if (dir == Direction::Right)
            return "Right";
        else if (dir == Direction::Bottom)
            return "Bottom";
        else
            return "Top";
    }

    void printPartitioningInfo() const;

    template <Direction direction> MPI_Request sendBorder(DataField &field) const {
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

        if (neighborRankNo(direction) != MPI_PROC_NULL) {
            // std::cout << "[" << ownRankNo_ << "] send " << dirToStr(direction) << " to [" << neighborRankNo(direction) << "]: ";
            // std::cout << "n=" << 1 << ", size=[" << field.cols() << ", " << field.rows() << "]" << ", start i=" << i << ",j=" << j << ", id=" <<
            // field.getID() << "\n";
        }

        MPI_Isend(&field(i, j), count, type, neighborRankNo(direction), field.getID(), cartComm_, &sendReq);

        MPI_Type_free(&coltype);

        return sendReq;
    }

    template <Direction direction> MPI_Request recvBorder(DataField &field) const {
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

        if (neighborRankNo(direction) != MPI_PROC_NULL) {
            // std::cout << "[" << ownRankNo_ << "] recv " << dirToStr(direction) << " to [" << neighborRankNo(direction) << "]: ";
            // std::cout << "n=" << count << ", size=[" << field.cols() << ", " << field.rows() << "]" << ", start i=" << i << ",j=" << j << ", id="
            // << field.getID() << "\n";
        }

        MPI_Irecv(&field(i, j), count, type, neighborRankNo(direction), field.getID(), cartComm_, &recvReq);

        MPI_Type_free(&coltype);

        return recvReq;
    }
};
