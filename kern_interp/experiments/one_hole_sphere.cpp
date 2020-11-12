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
#include "kern_interp/boundaries/concentric_spheres.h"
#include "kern_interp/boundaries/donut.h"

namespace kern_interp {

double laplace_error3d(const ki_Mat& domain,
                       const std::vector<double>& domain_points,
                       Boundary * boundary, BoundaryCondition bc) {
  double diff_norm = 0;
  double norm_of_true = 0;
  for (int i = 0; i < domain_points.size(); i += 3) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    double x2 = domain_points[i + 2];
    PointVec x(x0, x1, x2);
    PointVec center(0.0, 0.0, 0.0);
    double r = (x - center).norm();
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    double potential;
    if (bc == BoundaryCondition::ELECTRON_3D) {
      potential = -1.0 / (4.0 * M_PI * sqrt(pow(x0 + 3, 2) + pow(x1+2 ,
                                            2) + pow(x2+2 , 2)));
    } else if (bc == BoundaryCondition::LAPLACE_CHECK_3D) {
      double c2 = (-2.0 / 9.0);
      double c1 = 3 - c2;
      potential = (c1) + (c2 / r);
    }
    diff_norm += pow(potential - domain.get(i / 3, 0), 2);
    norm_of_true += pow(potential, 2);
  }
  diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
  return diff_norm;
}

void run_one_hole_sphere() {
  srand(0);
  int fact_threads = 8;
  double id_tol = 1e-4;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Sphere());
  // Boundary condition is flow past noslip interior hole.
  boundary->initialize(pow(2, 7),  BoundaryCondition::ELECTRON_3D);


  std::vector<double> old_domain_points, domain_points;
  get_domain_points3d(6,  &old_domain_points, boundary.get(), 0.1, 1);
  for (int i = 0; i < old_domain_points.size(); i += 3) {
    if (boundary->is_in_domain(PointVec(old_domain_points[i],
                                        old_domain_points[i + 1],
                                        old_domain_points[i + 2]))) {
      domain_points.push_back(old_domain_points[i]);
      domain_points.push_back(old_domain_points[i + 1]);
      domain_points.push_back(old_domain_points[i + 2]);
    }
  }

  Kernel kernel(1, 3, Kernel::Pde::LAPLACE, boundary.get(), domain_points);
  double cstart = omp_get_wtime();
  kernel.compute_diag_entries_3dlaplace(boundary.get());
  double cend = omp_get_wtime();
  std::cout << "Time to compute diagonal " << (cend - cstart) << std::endl;
  double err = solve_err(kernel, boundary.get(), id_tol);
  // ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
  //                                      id_tol, fact_threads, domain_points);
  // double err = laplace_error3d(sol, domain_points, boundary.get(),
  //                              BoundaryCondition::ELECTRON_3D);
  std::cout<<"Err "<<err<<std::endl;
  // std::ofstream sol_out;
  // sol_out.open("output/data/ie_solver_solution.txt");
  // int points_index = 0;
  // for (int i = 0; i < sol.height(); i += 3) {
  //   sol_out << domain_points[points_index] << "," <<
  //           domain_points[points_index + 1] << ","
  //           << domain_points[points_index + 2] << ",";
  //   points_index += 3;
  //   sol_out << sol.get(i, 0) << "," << sol.get(i + 1, 0)
  //           << "," << sol.get(i + 2, 0)
  //           << std::endl;
  // }
  // sol_out.close();
}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  openblas_set_num_threads(1);
  kern_interp::run_one_hole_sphere();
  return 0;
}

