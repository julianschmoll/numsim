#include "outputWriter.h"

OutputWriter::OutputWriter(std::shared_ptr<StaggeredGrid> grid, const Partitioning &partitioning)
    : grid_(grid), partitioning_(partitioning), fileNo_(0) {
    // create "out" subdirectory if it does not yet exist
    int returnValue = system("mkdir -p out");
    if (returnValue != 0)
        std::cout << "Could not create subdirectory \"out\"." << std::endl;
}
