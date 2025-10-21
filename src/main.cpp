// #include "outputWriter.h"
#include "settings.h"

//#include <iostream>

int main(int argc, char *argv[]) {
  // we need an input file being specified
  if (argc == 1) {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
    return EXIT_FAILURE;
  }

  // read in the first argument
  std::string filename = argv[1];

  Settings settings;
  settings.loadFromFile(filename);

#ifndef NDEBUG
  settings.printSettings();
#endif

  // then we call the solver here :)
  return EXIT_SUCCESS;
}
