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
#include "kern_interp/boundaries/sphere.h"
#include "kern_interp/boundaries/donut.h"

namespace kern_interp {


void run_one_hole_sphere() {

  int num_threads = 8;
  double id_tol = 1e-3;

  std::unique_ptr<Boundary> boundary3d =
    std::unique_ptr<Boundary>(new Sphere());
  Hole hole3d;
  hole3d.center = PointVec(0.5, 0.5, 0.5);
  hole3d.radius = 0.1;
  boundary3d->holes.push_back(hole3d);
  boundary3d->initialize(pow(2, 13),  BoundaryCondition::STOKES_3D_MIX);
  std::cout << "num nodes " << boundary3d->weights.size() << std::endl;
  QuadTree quadtree3d;
  quadtree3d.initialize_tree(boundary3d.get(), 3, 3);
  std::vector<double> domain_points3d;
  int domain_size = 10;

  get_domain_points3d(domain_size, &domain_points3d, boundary3d.get(), 0, 1);

  Kernel kernel3d(3, 3, Kernel::Pde::STOKES, boundary3d.get(), domain_points3d);
  // TODO(John) this should be part of kernel init
  double start = omp_get_wtime();
  kernel3d.compute_diag_entries_3dstokes(boundary3d.get());
  double end = omp_get_wtime();
  std::cout << "Compute diag takes " << (end - start) << std::endl;
  ki_Mat sol3d = boundary_integral_solve(kernel3d, *(boundary3d.get()),
                                         &quadtree3d,
                                         id_tol, num_threads, domain_points3d);

  double err = stokes_err_3d(sol3d, domain_points3d, boundary3d.get(),
                             hole3d.radius, STOKES_MIXER);
  std::cout << "Err: " << err << std::endl;
// double err = stokes_err_3d(sol, domain_points, boundary.get(),
//     hole.radius, STOKES_MIXER);
// //   std::cout<<"err of "<<err<<std::endl;
//   std::ofstream sol_out;
//   sol_out.open("output/data/ie_solver_solution.txt");
//   int points_index = 0;
//   for (int i = 0; i < sol3d.height(); i += 3) {
//     sol_out << domain_points3d[points_index] << "," <<
//             domain_points3d[points_index + 1] << "," <<
//             domain_points3d[points_index + 2] << ",";
//     points_index += 3;
//     sol_out << sol3d.get(i, 0) << "," <<
//             sol3d.get(i + 1, 0) << "," <<
//             sol3d.get(i + 2, 0) << std::endl;
//   }
//   sol_out.close();
}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  openblas_set_num_threads(1);
  kern_interp::run_one_hole_sphere();
  return 0;
}

