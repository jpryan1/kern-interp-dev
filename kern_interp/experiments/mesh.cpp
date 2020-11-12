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
#include "kern_interp/boundaries/mesh.h"
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
      potential = -1.0 / (4.0 * M_PI * sqrt(pow(x0, 2) + pow(x1 - 5,
                                            2) + pow(x2 - 5, 2)));
    } else if (bc == BoundaryCondition::ALL_ONES) {
      double c2 = (-2.0 / 9.0);
      double c1 = 3 - c2;
      potential = 1;//(c1) + (c2 / r);
    }
    diff_norm += pow(potential - domain.get(i / 3, 0), 2);
    norm_of_true += pow(potential, 2);
  }
  diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
  return diff_norm;
}


void get_mesh_domain_points(std::vector<double>* domain_points, double min,
                            double max) {
  // x val 0 everywhere, just yz plane
  double delt = (max - min) / 100.0;
  for (double y = min; y < max; y += delt) {
    for (double z = min; z < max; z += delt) {
      domain_points->push_back(0);
      domain_points->push_back(y);
      domain_points->push_back(z);
    }
  }
}


void run_mesh() {
  srand(0);
  int num_threads = 8;
  double id_tol = 1e-4;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Mesh());
  // Boundary condition is flow past noslip interior hole.
  boundary->initialize(pow(2, 7),  BoundaryCondition::ELECTRON_3D);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), 1, 3);

  std::vector<double> old_domain_points, domain_points;
  get_mesh_domain_points(&domain_points, quadtree.min, quadtree.max);

  Kernel kernel(1, 3, Kernel::Pde::LAPLACE, boundary.get(),
                domain_points);
  kernel.compute_diag_entries_3dlaplace(boundary.get());

  SkelFactorization skel_factorization(id_tol, num_threads, false);
  ki_Mat K_domain = kernel.forward();

  ki_Mat solution;
  if (boundary->holes.size() > 0) {
    ki_Mat U = initialize_U_mat(kernel.pde, boundary->holes, boundary->points,
                                kernel.domain_dimension);
    ki_Mat Psi = initialize_Psi_mat(kernel.pde, boundary->holes, *boundary,
                                    kernel.domain_dimension);
    ki_Mat U_forward = initialize_U_mat(kernel.pde, boundary->holes,
                                        kernel.domain_points, kernel.domain_dimension);
    quadtree.U = U;
    quadtree.Psi = Psi;
    skel_factorization.skeletonize(kernel, &quadtree);
    ki_Mat mu, alpha;
    linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu,
                 &alpha);
    solution = K_domain * mu + U_forward * alpha;
  } else {
    skel_factorization.skeletonize(kernel, &quadtree);
    ki_Mat mu;
    linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu,
                 nullptr);
    solution = K_domain * mu ;
  }

  for (int i = 0; i < kernel.domain_points.size(); i += kernel.domain_dimension) {
    std::vector<double> vec;

    for (int j = 0; j < kernel.domain_dimension; j++) {
      vec.push_back(kernel.domain_points[i + j]);
    }
    PointVec point(vec);
    if (!boundary->is_in_domain(point)) {
      int pt_idx = i / kernel.domain_dimension;
      for (int j = 0; j < kernel.solution_dimension; j++) {
        solution.set(kernel.solution_dimension * pt_idx + j, 0, 0.);
      }
    }
  }

  std::cout << "err " << laplace_error3d(solution,
                                         kernel.domain_points,
                                         boundary.get(),
                                         BoundaryCondition::ELECTRON_3D) << std::endl;
  std::ofstream sol_out;
  sol_out.open("output/data/cross_section_sol.txt");
  int points_index = 0;
  for (int i = 0; i < solution.height(); i++) {
    sol_out << domain_points[points_index + 1] << "," <<
            domain_points[points_index + 2] << ",";
    points_index += 3;
    sol_out << solution.get(i, 0)
            << std::endl;
  }
  sol_out.close();

}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  openblas_set_num_threads(1);
  kern_interp::run_mesh();
  return 0;
}

