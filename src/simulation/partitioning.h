#pragma once

#include <array>
#include <mpi.h>

// This avoids magic numbers in shift and can be used in direction as well
enum class Direction {
    Left, Right, Top, Bottom
};

// This can be used to iterate over every direction and avoid redundancy in code :)
constexpr std::array<Direction, 4> directions() {
    return {Direction::Left, Direction::Right, Direction::Bottom, Direction::Top};
}

class Partitioning
{
    /// Dimensions of the problem: 2D grid with x,y directions.
    static constexpr int dimensions_ = 2;

    std::array<int,2> nCellsLocal_ = {};
    std::array<int,2> nCellsGlobal_ = {};

    /// Number of ranks in each direction.
    std::array<int,2> partitions_ = {};

    /// MPI communicator with cartesian coordinate information attached.
    MPI_Comm cartComm_ = nullptr;

    /// Node offsets
    int offsetX_ = 0;
    int offsetY_ = 0;

    int ownRankNo_ = 0;

    std::array<int, dimensions_> rankCoordinates_ = {};

public:

  //! compute partitioning, set internal variables
  void initialize(std::array<int,2> nCellsGlobal);

  //! get the local number of cells in the own subdomain
  std::array<int,2> nCellsLocal() const;

  //! get the global number of cells in the whole computational domain
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nCellsGlobal() const;

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
  std::array<int,2> nodeOffset() const;

  std::array<int,2> getCurrentRankCoords() const;
};
