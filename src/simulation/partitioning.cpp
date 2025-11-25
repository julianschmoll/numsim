#include "partitioning.h"

int Partitioning::ownRankNo() const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

//! number of MPI ranks
int nRanks() const;