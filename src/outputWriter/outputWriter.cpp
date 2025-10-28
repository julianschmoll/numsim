#include "outputWriter.h"

outputWriter::outputWriter(const std::shared_ptr<discretization>& discretization)
    : discretization_(discretization), fileNo_(0) {
    // create "out" subdirectory if it does not yet exist
    int returnValue = system("mkdir -p out");
    if (returnValue != 0)
        std::cout << "Could not create subdirectory \"out\"." << std::endl;
}
