#pragma once

#include <array>
#include <iostream>

/**
 * All settings that parametrize a simulation run.
 */
struct Settings {
  // number of cells in x and y direction
  std::array<int, 2> nCells;

  // physical size of the domain
  std::array<double, 2> physicalSize;

  // reynolds number
  double re = 1000;

  // end time of the simulation
  double endTime = 10.0;

  // safety factor for time step width
  double tau = 0.5;

  // maximum time step width
  double maximumDt = 0.1;

  // external forces
  std::array<double, 2> g{0., 0.};

  // if the donor cell scheme schould be used
  bool useDonorCell = false;

  // factor for donor-cell scheme
  double alpha = 0.5;

  // prescribed values of u,v at bottom of domain
  std::array<double, 2> dirichletBcBottom;

  // prescribed values of u,v at top of domain
  std::array<double, 2> dirichletBcTop;

  // prescribed values of u,v at left of domain
  std::array<double, 2> dirichletBcLeft;

  // prescribed values of u,v at right of domain
  std::array<double, 2> dirichletBcRight;

  // which pressure solver to use, "GaussSeidel" or "SOR"
  std::string pressureSolver = "SOR";

  // overrelaxation factor
  double omega = 1.0;

  // tolerance for the residual in the pressure solver
  double epsilon = 1e-5;

  // maximum number of iterations in the solver
  int maximumNumberOfIterations = 1e5;

  /**
  * Parses a text file with settings. 
  *
  * Each line Contains settings tuple like "setting" = "value".
  *
  * @param filename The name of the file to parse
  */
  void loadFromFile(const std::string& filename);

  /**
  * Outputs Settings to console.
  */
  void printSettings() const;
};
