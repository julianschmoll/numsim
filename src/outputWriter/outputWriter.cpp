#include "outputWriter.h"

OutputWriter::OutputWriter(std::shared_ptr<StaggeredGrid> grid, const Partitioning &partitioning, const std::string &folderName)
    : grid_(grid), partitioning_(partitioning), fileNo_(0), folderName_(folderName) {
    // create "out" subdirectory if it does not yet exist
    std::string command = "mkdir -p " + folderName_;
    int returnValue = system(command.c_str());
    if (returnValue != 0) {
        std::cout << "Could not create subdirectory \"out\"." << std::endl;
    }
}
