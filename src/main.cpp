#include "write_paraview_output.h"

#include <iostream>
#include <cstdlib>
#include <vector>

int main(int argc, char *argv[]) {
    // write 5 output files
    for (int i = 0; i < 5; i++) {
        writeParaviewOutput(i);
    }
    return EXIT_SUCCESS;
}
