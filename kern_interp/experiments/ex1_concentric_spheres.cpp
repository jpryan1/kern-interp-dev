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

namespace kern_interp {


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
      potential = (-1.0 / (4.0 * M_PI * sqrt(pow(x0, 2) + pow(x1 ,
                                             2) + pow(x2 + 1.5, 2))));

      // potential -= (-1.0 / (4.0 * M_PI * sqrt(pow(x0, 2) + pow(x1 ,
      //                                         2) + pow(x2 + 1.1, 2))));

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


void run_ex1_concentric_spheres(bool is_stokes) {
  srand(0);
  int fact_threads = 8;
  double id_tol = 1e-4;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Sphere());
  // Boundary condition is flow past noslip interior hole.
  BoundaryCondition bc =  BoundaryCondition::ELECTRON_3D;
  int solution_dim = 1;
  Kernel::Pde pde = Kernel::Pde::LAPLACE;

  if (is_stokes) {
    bc = BoundaryCondition::STOKES_3D_MIX;
    solution_dim = 3;
    pde = Kernel::Pde::STOKES;
  }

  boundary->initialize(0, bc);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), solution_dim, 3);
  std::vector<double> old_domain_points, domain_points;
  get_mesh_domain_points(&domain_points, quadtree.min, quadtree.max);

  Kernel kernel(solution_dim, 3, pde, boundary.get(),
                domain_points);
  if (is_stokes) {
    kernel.compute_diag_entries_3dstokes(boundary.get());
  } else {
    kernel.compute_diag_entries_3dlaplace(boundary.get());
  }
  ki_Mat solution = boundary_integral_solve(kernel, *(boundary.get()),
                    &quadtree, id_tol, fact_threads, domain_points);
  if (!is_stokes) {
    std::cout << "err " << laplace_error3d(solution, kernel.domain_points,
                                           boundary.get(),
                                           BoundaryCondition::ELECTRON_3D)
              << std::endl;
  }
  std::ofstream sol_out;
  sol_out.open("output/data/cross_section_sol.txt");
  for (int i = 0; i < solution.height(); i += 3) {
    sol_out << domain_points[i + 1] << "," <<
            domain_points[i + 2] << ",";
    if (is_stokes) {
      sol_out << solution.get(i + 1,
                              0) << "," << solution.get(i + 2, 0)
              << std::endl;

    } else {
      sol_out << solution.get(i, 0)
              << std::endl;
    }
  }
  sol_out.close();
}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  openblas_set_num_threads(1);
  kern_interp::run_ex1_concentric_spheres(true);
  return 0;
}
