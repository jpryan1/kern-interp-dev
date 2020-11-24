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
#include "kern_interp/boundaries/pipe_mesh.h"

namespace kern_interp {


void get_mesh_domain_points(std::vector<double>* domain_points, double min,
                            double max) {
  // x val 0 everywhere, just yz plane
  double delt = (max - min) / 50.0;
  for (double y = min; y < max; y += delt) {
    for (double z = min; z < max; z += delt) {
      domain_points->push_back(0);
      domain_points->push_back(y);
      domain_points->push_back(z);
    }
  }
}

void run_ex2_pipe_w_ellipsoid_holes(bool is_stokes) {
  srand(0);
  int fact_threads = 8;
  double id_tol = 1e-3;
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new PipeMesh());
  // Boundary condition is flow past noslip interior hole.
  BoundaryCondition bc =  BoundaryCondition::LAPLACE_PIPE_HOLES;
  int solution_dim = 1;
  Kernel::Pde pde = Kernel::Pde::LAPLACE;

  if (is_stokes) {
    bc = BoundaryCondition::STOKES_PIPE_HOLES;
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




  boundary->perturbation_parameters[0] = 0.05;
  boundary->initialize(2, BoundaryCondition::STOKES_PIPE_HOLES);
  kernel.update_data(boundary.get());
  kernel.compute_diag_entries_3dstokes(boundary.get());
  quadtree.perturb(*boundary.get());
  solution = boundary_integral_solve(kernel, *(boundary.get()),
                                     &quadtree, id_tol, fact_threads, domain_points);
  std::ofstream sol_out;
  sol_out.open("output/data/cross_section_sol.txt");
  for (int i = 0; i < domain_points.size(); i += 3) {
    sol_out << domain_points[i + 1] << "," <<
            domain_points[i + 2] << ",";

    if (is_stokes) {
      sol_out << solution.get(i + 1, 0) << "," << solution.get(i + 2, 0)
              << std::endl;
    } else {
      sol_out << solution.get(i / 3, 0)
              << std::endl;
    }
  }

  QuadTree fresh;
  fresh.initialize_tree(boundary.get(), 3,  3);
  ki_Mat new_sol = boundary_integral_solve(kernel, *(boundary.get()), &fresh,
                   id_tol, fact_threads, domain_points);
  std::cout << "Err: " << (new_sol - solution).vec_two_norm() /
            new_sol.vec_two_norm() << std::endl;
  sol_out.close();
}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  openblas_set_num_threads(1);
  kern_interp::run_ex2_pipe_w_ellipsoid_holes(true);
  return 0;
}
