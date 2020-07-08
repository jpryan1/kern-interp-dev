// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <string.h>
#include <fstream>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include "kern_interp/ki_mat.h"
#include "kern_interp/skel_factorization/skel_factorization.h"
#include "kern_interp/quadtree/quadtree.h"
#include "kern_interp/kernel/kernel.h"
#include "kern_interp/linear_solve.h"
#include "kern_interp/boundaries/ex1boundary.h"

namespace kern_interp {


void run_experiment1(int num_boundary_points, int fact_threads) {
  double id_tol = 1e-6;
  Kernel::Pde pde = Kernel::Pde::STOKES;
  // int num_boundary_points = pow(2, 13);
  int domain_size = 350;
  // int domain_size = 10;
  int domain_dimension = 2;
  int solution_dimension = 2;

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex1Boundary());
  boundary->initialize(num_boundary_points, BoundaryCondition::DEFAULT);
  double start = omp_get_wtime();

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), solution_dimension,
                           domain_dimension);
  double x_min = boundary->points[0], x_max = boundary->points[0],
         y_min = boundary->points[1], y_max = boundary->points[1];
  for (int i = 0; i < boundary->points.size(); i += 2) {
    x_min = std::min(x_min, boundary->points[i]);
    x_max = std::max(x_max, boundary->points[i]);
    y_min = std::min(y_min, boundary->points[i + 1]);
    y_max = std::max(y_max, boundary->points[i + 1]);
  }
  std::vector<double> domain_points;
  get_domain_points(domain_size, &domain_points, x_min, x_max, y_min,
                    y_max);
  Kernel kernel(solution_dimension, domain_dimension,
                pde, boundary.get(), domain_points);

  ki_Mat solution = boundary_integral_solve(kernel, *(boundary.get()),
                    &quadtree, id_tol, fact_threads, domain_points);
  // boundary->perturbation_parameters[0] += 0.2;
  // boundary->initialize(num_boundary_points, BoundaryCondition::DEFAULT);
  // quadtree.perturb(*boundary.get());
  // kernel.update_data(boundary.get());
  // solution = boundary_integral_solve(kernel, *(boundary.get()),
  //                                    &quadtree, id_tol, fact_threads, domain_points);

  std::ofstream sol_out;
  sol_out.open("output/data/ie_solver_solution.txt");
  int points_index = 0;
  for (int i = 0; i < solution.height(); i += 2) {
    sol_out << domain_points[points_index] << "," <<
            domain_points[points_index + 1] << ",";
    points_index += 2;
    sol_out << solution.get(i, 0) << "," << solution.get(i + 1, 0)
            << std::endl;
  }
  sol_out.close();
  std::ofstream bound_out;
  bound_out.open("output/data/ie_solver_boundary.txt");
  for (int i = 0; i < boundary->points.size(); i += 2) {
    bound_out << boundary->points[i] << "," << boundary->points[i + 1]
              << std::endl;
  }
  bound_out.close();
}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  openblas_set_num_threads(1);
  // for (int num = 15; num <= 19; num++) {
  //   std::cout << "pow " << num << std::endl;
  //   for (int k = 0; k < 3; k++) {
  //     for (int thrpow = 0; thrpow <= 3; thrpow++) {
  // kern_interp::run_experiment1(pow(2, num), pow(2, thrpow));
  kern_interp::run_experiment1(pow(2, 12), pow(2, 2));
  //     }
  //   }
  // }
  return 0;
}
