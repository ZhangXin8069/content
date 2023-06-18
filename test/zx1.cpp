#include <complex>
#include <iostream>
#include <vector>

#include "bicgstab_solver.h"

int main() {
  // Set up matrix and right-hand side vector
  std::vector<std::vector<std::complex<double>>> A = {{1.0, 2.0}, {2.0, 3.0}};
  std::vector<std::complex<double>> b = {1.0, 2.0};

  // Create solver and solve linear system
  BiCGStabSolver<std::complex<double>> solver(A, b);
  std::vector<std::complex<double>> x = solver.Solve();

  // Print solution vector
  std::cout << "Solution: ";
  for (const auto& element : x) {
    std::cout << element << " ";
  }
  std::cout << std::endl;

  return 0;
}