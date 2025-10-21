#include "settings.h"

// #include <fstream>
// #include <iomanip>

void Settings::loadFromFile(std::string filename) { return; }

void Settings::printSettings() {
  std::cout << "Settings: " << std::endl
            
            << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1]
            << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl

            << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0]
            << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt
            << std::endl

            << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << ","
            << dirichletBcBottom[1] << ")"

            << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1]
            << ")"

            << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1]
            << ")"

            << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1]
            << ")" << std::endl

            << "  useDonorCell: " << std::boolalpha << useDonorCell
            << ", alpha: " << alpha << std::endl

            << "  pressureSolver: " << pressureSolver << ", omega: " << omega
            << ", epsilon: " << epsilon
            << ", maximumNumberOfIterations: " << maximumNumberOfIterations
            << std::endl;
}
