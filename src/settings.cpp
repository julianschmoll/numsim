#include "settings.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

void Settings::loadFromFile(const std::string &filename) {
    std::cout << "loading settings from " << filename << std::endl;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open: " + filename);
    }

    std::string line;
    std::unordered_map<std::string, std::string> settings;

    while (std::getline(file, line)) {
        // those are comments
        if (line.find('#') == 0)
            continue;

        // strip inline comments
        line = line.substr(0, line.find('#'));

        // strip whitespace
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

        if (line.empty())
            continue;

        auto equalPos = line.find('=');

        if (equalPos == std::string::npos)
            continue;

        std::string key = line.substr(0, equalPos);
        std::string value = line.substr(equalPos + 1);
        settings[key] = value;
    }

    // this is somewhat nasty but I figured it would be the quickest
    // IDK, which cpp standard we use, std::unordered_map only exists
    // in cpp 20 or above, that's why I used count...
    if (settings.count("nCellsX"))
        nCells[0] = std::stoi(settings["nCellsX"]);
    if (settings.count("nCellsY"))
        nCells[1] = std::stoi(settings["nCellsY"]);

    if (settings.count("physicalSizeX"))
        physicalSize[0] = std::stod(settings["physicalSizeX"]);
    if (settings.count("physicalSizeY"))
        physicalSize[1] = std::stod(settings["physicalSizeY"]);

    if (settings.count("re"))
        re = std::stod(settings["re"]);
    if (settings.count("endTime"))
        endTime = std::stod(settings["endTime"]);
    if (settings.count("tau"))
        tau = std::stod(settings["tau"]);
    if (settings.count("maximumDt"))
        maximumDt = std::stod(settings["maximumDt"]);

    if (settings.count("gX"))
        g[0] = std::stod(settings["gX"]);
    if (settings.count("gY"))
        g[1] = std::stod(settings["gY"]);

    if (settings.count("useDonorCell"))
        useDonorCell = (settings["useDonorCell"] == "true");

    if (settings.count("alpha"))
        alpha = std::stod(settings["alpha"]);

    if (settings.count("dirichletBottomX"))
        dirichletBcBottom[0] = std::stod(settings["dirichletBottomX"]);
    if (settings.count("dirichletBottomY"))
        dirichletBcBottom[1] = std::stod(settings["dirichletBottomY"]);

    if (settings.count("dirichletTopX"))
        dirichletBcTop[0] = std::stod(settings["dirichletTopX"]);
    if (settings.count("dirichletTopY"))
        dirichletBcTop[1] = std::stod(settings["dirichletTopY"]);

    if (settings.count("dirichletLeftX"))
        dirichletBcLeft[0] = std::stod(settings["dirichletLeftX"]);
    if (settings.count("dirichletLeftY"))
        dirichletBcLeft[1] = std::stod(settings["dirichletLeftY"]);

    if (settings.count("dirichletRightX"))
        dirichletBcRight[0] = std::stod(settings["dirichletRightX"]);
    if (settings.count("dirichletRightY"))
        dirichletBcRight[1] = std::stod(settings["dirichletRightY"]);

    if (settings.count("pressureSolver")) {
        if (settings["pressureSolver"] == "SOR") {
            pressureSolver = IterSolverType::SOR;
        } else if (settings["pressureSolver"] == "GaussSeidel") {
            pressureSolver = IterSolverType::GaussSeidel;
        } else if (settings["pressureSolver"] == "CG") {
            pressureSolver = IterSolverType::CG;
        }  else {
            throw std::runtime_error("Unnknown solver type for \"pressureSolver\": use \"SOR\" or \"GaussSeidel\".");
        }
    }
    if (settings.count("omega"))
        omega = std::stod(settings["omega"]);
    if (settings.count("epsilon"))
        epsilon = std::stod(settings["epsilon"]);
    if (settings.count("maximumNumberOfIterations"))
        maximumNumberOfIterations = static_cast<int>(std::stod(settings["maximumNumberOfIterations"]));
}

void Settings::printSettings() const {
    std::cout << "Settings: " << std::endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl

              << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt
              << std::endl

              << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"

              << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"

              << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"

              << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl

              << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl

              << "  pressureSolver: " << (pressureSolver == IterSolverType::SOR ? "SOR" : "GaussSeidel") << ", omega: " << omega
              << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}
