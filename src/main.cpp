#include "settings.h"
#include "simulation/simulation.h"
#include "simulation/parallelSimulation.h"
#include "macros.h"
#include "simulation/partitioning.h"

#include <mpi.h>
#include <iostream>

// TODO: Fas funktioniert perfekt:
void testCode() {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
        // --- SENDER (P0): Column Vector ---
        int matrix_send[4 * 4] = {
            10, 11, 12, 13,
            14, 15, 16, 17,
            18, 19, 20, 21,
            22, 23, 24, 25
        };
        
        MPI_Datatype send_column_type;
        // Stride is 4 (4 integers)
        MPI_Type_vector(4, 1, 4, MPI_INT, &send_column_type);
        MPI_Type_commit(&send_column_type);

        printf("P0: Sending data as a COLUMN vector (stride 4).\n");
        // We send 1 instance of the column type, pointing to the start of the matrix
        MPI_Send(&(matrix_send[0]), 1, send_column_type, 1, 0, MPI_COMM_WORLD);

        MPI_Type_free(&send_column_type);

    } else if (world_rank == 1) {
        // --- RECEIVER (P1): Row Vector ---
        int matrix_recv[4 * 2] = {0, 0, 0, 0, 0, 0, 0, 0}; // Initialize to zero
        
        MPI_Datatype recv_row_type;
        // Stride is 1 (contiguous memory for the first row)
        MPI_Type_vector(4, 1, 2, MPI_INT, &recv_row_type);
        MPI_Type_commit(&recv_row_type);
        
        printf("P1: Receiving data as a ROW vector (stride 1).\n");
        // We receive 1 instance of the row type, pointing to the start of the receive matrix's first row
        MPI_Recv(&(matrix_recv[1]), 1, recv_row_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Type_free(&recv_row_type);

        printf("P1: Resulting first row of the matrix: { ");
        // Print only the first row where data was received
        for (int i = 0; i < 4; i++) {
            printf("%d%s", matrix_recv[2 * i + 1], (i == 4 - 1) ? "" : ", ");
        }
        printf(" }\n");
        // Expected output: { 10, 14, 18, 22 }
    }
}


int main(int argc, char *argv[]) {
    // we need an input file being specified
    if (argc == 1) {
        std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
        return EXIT_FAILURE;
    }
    const std::string filename = argv[1];

    Settings settings;
    settings.loadFromFile(filename);

    DEBUG(settings.printSettings());

    MPI_Init(&argc, &argv);

    int ownRankNo = 0;
    int nRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    //testCode();
    {
        ParallelSimulation simulation {};
        simulation.initialize(settings);
        simulation.run();
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
