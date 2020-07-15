// Copyright John Paul Ryan 2019
#include <omp.h>
#include <string.h>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <cassert>
#include "kern_interp/ki_mat.h"
#include "kern_interp/boundaries/circle.h"
#include "kern_interp/boundaries/annulus.h"
#include "kern_interp/boundaries/donut.h"
#include "kern_interp/boundaries/ex1boundary.h"
#include "kern_interp/boundaries/ex2boundary.h"
#include "kern_interp/boundaries/ex3boundary.h"
#include "kern_interp/boundaries/sphere.h"
#include "kern_interp/boundaries/cubic_spline.h"
#include "kern_interp/skel_factorization/skel_factorization.h"
#include "kern_interp/quadtree/quadtree.h"
#include "kern_interp/kernel/kernel.h"
#include "kern_interp/linear_solve.h"
#include "gtest/gtest.h"

namespace kern_interp {


void check_solve_err(const Kernel& kernel, Boundary* boundary) {
  double id_tol = 1e-6;

  EXPECT_LE(solve_err(kernel, boundary, id_tol), 10 * id_tol);
}


TEST(IeSolverTest, LaplaceCircleBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(pow(2, 10),  BoundaryCondition::SINGLE_ELECTRON);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, LaplaceNeumannCircleBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(pow(2, 10),  BoundaryCondition::ALL_ONES);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE_NEUMANN, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, StokesCircleBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(pow(2, 10),  BoundaryCondition::TANGENT_VEC);
  Kernel kernel(2, 2, Kernel::Pde::STOKES, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, LaplaceStarfishBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(pow(2, 10),  BoundaryCondition::SINGLE_ELECTRON);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, LaplaceNeumannStarfishBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(pow(2, 10),  BoundaryCondition::ALL_ONES);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE_NEUMANN, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, StokesStarfishBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new CubicSpline());
  boundary->initialize(pow(2, 10),  BoundaryCondition::TANGENT_VEC);
  Kernel kernel(2, 2, Kernel::Pde::STOKES, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}



TEST(IeSolverTest, LaplaceAnnulusBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(pow(2, 10),  BoundaryCondition::SINGLE_ELECTRON);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, LaplaceNeumannAnnulusBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(pow(2, 10),  BoundaryCondition::DEFAULT);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE_NEUMANN, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


TEST(IeSolverTest, StokesAnnulusBackwardError) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(pow(2, 10),  BoundaryCondition::TANGENT_VEC);
  Kernel kernel(2, 2, Kernel::Pde::STOKES, boundary.get(),
                std::vector <double>());
  check_solve_err(kernel, boundary.get());
}


double laplace_error(const ki_Mat& domain,
                     const std::vector<double>& domain_points,
                     Boundary * boundary) {
  double diff_norm = 0;
  double norm_of_true = 0;
  for (int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    PointVec x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    double potential = log(sqrt(pow(x0 + 3, 2) + pow(x1 + 2, 2))) / (2 * M_PI);

    diff_norm += pow(potential - domain.get(i / 2, 0), 2);
    norm_of_true += pow(potential, 2);
  }
  diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
  return diff_norm;
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
    PointVec center(0.5, 0.5, 0.5);
    double r = (x - center).norm();
    if (!boundary->is_in_domain(x)) {
      continue;
    }
    double potential;
    if (bc == BoundaryCondition::ELECTRON_3D) {
      potential = -1.0 / (4.0 * M_PI * sqrt(pow(x0 + 3, 2) + pow(x1 + 2,
                                            2) + pow(x2 + 2, 2)));
    } else if (bc == BoundaryCondition::LAPLACE_CHECK_3D) {
      double c2 = (-2.0 / 9.0);
      double c1 = 3 - c2;
      potential = (c1) + (c2 / r);
    }
    // std::cout<<potential<< " vs "<<domain.get(i / 3, 0)<<std::endl;
    diff_norm += pow(potential - domain.get(i / 3, 0), 2);
    norm_of_true += pow(potential, 2);
  }
  diff_norm = sqrt(diff_norm) / sqrt(norm_of_true);
  return diff_norm;
}
double laplace_neumann_error(const ki_Mat& domain,
                             const std::vector<double>& domain_points,
                             Boundary * boundary) {
  std::vector<double> res;
  std::vector<double> truth;

  for (int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    PointVec x(x0, x1);
    if (boundary->is_in_domain(x)) {
      res.push_back(domain.get(i / 2, 0));
      double r = sqrt(pow(x0 - 0.5, 2) + pow(x1 - 0.5, 2));
      truth.push_back(log(r));
    }
  }
  double res_avg = 0.;
  double truth_avg = 0.;
  for (int i = 0; i < res.size(); i++) {
    res_avg += res[i];
    truth_avg += truth[i];
  }
  res_avg /= res.size();
  truth_avg /= truth.size();

  for (int i = 0; i < res.size(); i++) {
    res[i] -= res_avg;
    truth[i] -= truth_avg;
  }

  double diff_norm = 0.;
  double truth_norm = 0.;
  for (int i = 0; i < res.size(); i++) {
    diff_norm += pow(res[i] - truth[i], 2);
    truth_norm += pow(truth[i], 2);
  }
  return sqrt(diff_norm) / sqrt(truth_norm);
}


TEST(IeSolverTest, LaplaceCircleAnalyticAgreementElectron) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(pow(2, 12),  BoundaryCondition::SINGLE_ELECTRON);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), 1, 2);
  std::vector<double> domain_points;
  get_domain_points(20, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE, boundary.get(), domain_points);
  ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
                                       1e-13, 4, domain_points);
  double err = laplace_error(sol, domain_points, boundary.get());
  EXPECT_LE(err, 10 * 1e-13);
}

// TEST(IeSolverTest, LaplaceSphereAnalyticAgreementElectron) {
//   srand(0);
//   openblas_set_num_threads(1);
//   int num_threads = 8;
//   double id_tol = 1e-6;std::cout << "HERE" << std::endl;

//   std::unique_ptr<Boundary> boundary =
//     std::unique_ptr<Boundary>(new Sphere());
//   Hole hole;
//   hole.center = PointVec(0.5, 0.5, 0.5);
//   hole.radius = 0.1;
//   boundary->holes.push_back(hole);
//   int mode = 4;  std::cout << "HERE" << std::endl;

//   boundary->initialize(mode,  BoundaryCondition::LAPLACE_CHECK_3D);

//   QuadTree quadtree;
//   quadtree.initialize_tree(boundary.get(), 1, 3);
//   std::vector<double> old_domain_points, domain_points;

//   get_domain_points3d(10, &old_domain_points, boundary.get(), 0.1, 1);
//   for (int i = 0; i < old_domain_points.size(); i += 3) {

//     if (boundary->is_in_domain(PointVec(old_domain_points[i],
//                                         old_domain_points[i + 1],
//                                         old_domain_points[i + 2]))) {
//       domain_points.push_back(old_domain_points[i]);
//       domain_points.push_back(old_domain_points[i + 1]);
//       domain_points.push_back(old_domain_points[i + 2]);
//     }
//   }
//   Kernel kernel(1, 3, Kernel::Pde::LAPLACE, boundary.get(), domain_points);

//   kernel.compute_diag_entries_3dlaplace(boundary.get());
//   ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
//                                        id_tol, num_threads, domain_points);

//   double err = laplace_error3d(sol, domain_points, boundary.get(),
//                                BoundaryCondition::LAPLACE_CHECK_3D);
//   EXPECT_LE(err, 10 * 1e-6);
// }


// TEST(IeSolverTest, StokesSphereAnalyticAgreement) {
//   srand(0);
//   openblas_set_num_threads(1);

//   int num_threads = 8;
//   double id_tol = 1e-6;

//   std::unique_ptr<Boundary> boundary3d =
//     std::unique_ptr<Boundary>(new Sphere());
//   Hole hole3d;
//   hole3d.center = PointVec(0.5, 0.5, 0.5);
//   hole3d.radius = 0.1;
//   boundary3d->holes.push_back(hole3d);
//   int mode = 4;

//   boundary3d->initialize(mode,  BoundaryCondition::STOKES_3D_MIX);

//   QuadTree quadtree3d;
//   quadtree3d.initialize_tree(boundary3d.get(), 3, 3);
//   std::vector<double> domain_points3d;
//   int domain_size = 10;

//   get_domain_points3d(domain_size, &domain_points3d, boundary3d.get(), 0, 1);

//   Kernel kernel3d(3, 3, Kernel::Pde::STOKES, boundary3d.get(), domain_points3d);
//   // TODO(John) this should be part of kernel init
//   kernel3d.compute_diag_entries_3dstokes(boundary3d.get());

//   ki_Mat sol3d = boundary_integral_solve(kernel3d, *(boundary3d.get()),
//                                          &quadtree3d,
//                                          id_tol, num_threads, domain_points3d);

//   double err = stokes_err_3d(sol3d, domain_points3d, boundary3d.get(),
//                              hole3d.radius, STOKES_MIXER);
//   EXPECT_LE(err, 10 * id_tol);
// }


TEST(IeSolverTest, LaplaceAnnulusAnalyticAgreementElectron) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Annulus());
  boundary->initialize(pow(2, 12),  BoundaryCondition::SINGLE_ELECTRON);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), 1, 2);
  std::vector<double> domain_points;
  get_domain_points(20, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE, boundary.get(), domain_points);
  ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
                                       1e-13,
                                       4, domain_points);
  double err = laplace_error(sol, domain_points, boundary.get());
  EXPECT_LE(err, 10 * 1e-13);
}


TEST(IeSolverTest, LaplaceNeumannDonutAnalyticAgreement) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Donut());
  boundary->initialize(pow(2, 12),  BoundaryCondition::DEFAULT);

  // TODO(John) fix this hack
  int hole_nodes = pow(2, 12) / 5;
  int num_points = pow(2, 12) + hole_nodes * 1;

  boundary->boundary_values = ki_Mat(num_points, 1);
  for (int i = 0; i < pow(2, 12); i++) {
    boundary->boundary_values.set(i, 0, 1);
  }
  for (int i = pow(2, 12); i < num_points; i++) {
    boundary->boundary_values.set(i, 0, -2.);
  }

  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), 1, 2);
  std::vector<double> domain_points;
  get_domain_points(20, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  Kernel kernel(1, 2, Kernel::Pde::LAPLACE_NEUMANN, boundary.get(),
                domain_points);
  ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
                                       1e-13,
                                       4, domain_points);
  double err = laplace_neumann_error(sol, domain_points, boundary.get());
  EXPECT_LE(err, 10 * 1e-13);
}


TEST(IeSolverTest, StokesCircleAnalyticAgreementTangent) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Circle());
  boundary->initialize(pow(2, 12),  BoundaryCondition::TANGENT_VEC);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), 2, 2);
  std::vector<double> domain_points;
  get_domain_points(20, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  Kernel kernel(2, 2, Kernel::Pde::STOKES, boundary.get(), domain_points);
  ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
                                       1e-13,
                                       4, domain_points);
  ki_Mat stokes_sol = stokes_true_sol(domain_points, boundary.get(), -1, 2);

  double err = (stokes_sol - sol).vec_two_norm() / stokes_sol.vec_two_norm();
  EXPECT_LE(err, 10 * 1e-13);
}


TEST(IeSolverTest, StokesDonutAnalyticAgreementTangent) {
  srand(0);
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Donut());
  boundary->initialize(pow(2, 12),  BoundaryCondition::TANGENT_VEC);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), 2, 2);
  std::vector<double> domain_points;
  get_domain_points(20, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min, quadtree.max);
  Kernel kernel(2, 2, Kernel::Pde::STOKES, boundary.get(), domain_points);
  ki_Mat sol = boundary_integral_solve(kernel, *(boundary.get()), &quadtree,
                                       1e-13,
                                       4, domain_points);
  ki_Mat stokes_sol = stokes_true_sol(domain_points, boundary.get(), -1, 2);

  double err = (stokes_sol - sol).vec_two_norm() / stokes_sol.vec_two_norm();
  EXPECT_LE(err, 10 * 1e-13);
}


TEST(IeSolverTest, Ex2UpdateLosesNoAcc) {
  double id_tol = 1e-6;
  Kernel::Pde pde = Kernel::Pde::STOKES;
  int num_boundary_points = pow(2, 16);
  int domain_size = 20;
  int domain_dimension = 2;
  int solution_dimension = 2;
  int fact_threads = 4;

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex2Boundary());
  boundary->initialize(num_boundary_points,
                       BoundaryCondition::DEFAULT);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), solution_dimension,
                           domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(domain_size, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  Kernel kernel(solution_dimension, domain_dimension,
                pde, boundary.get(), domain_points);

  for (int frame = 0; frame < 10; frame++) {
    int rand_idx = floor(8 * (rand() / (0. + RAND_MAX)));
    boundary->perturbation_parameters[rand_idx] = 0.35
        + 0.3 * (frame / 10.0);
    boundary->initialize(num_boundary_points, BoundaryCondition::DEFAULT);
    quadtree.perturb(*boundary.get());
    kernel.update_data(boundary.get());
    ki_Mat solution = boundary_integral_solve(kernel, *(boundary.get()),
                      &quadtree, id_tol, fact_threads, domain_points);

    QuadTree fresh;
    fresh.initialize_tree(boundary.get(), 2,  2);
    ki_Mat new_sol = boundary_integral_solve(kernel, *(boundary.get()), &fresh,
                     id_tol, fact_threads, domain_points);
    std::cout << "err: " << (new_sol - solution).vec_two_norm() /
              new_sol.vec_two_norm() << std::endl;
    ASSERT_LE((new_sol - solution).vec_two_norm() / new_sol.vec_two_norm(),
              20 * id_tol);
  }
}


TEST(IeSolverTest, TreeCopyGivesSameAnswer) {
  double id_tol = 1e-6;
  Kernel::Pde pde = Kernel::Pde::STOKES;
  int num_boundary_points = pow(2, 12);
  int domain_size = 20;
  int domain_dimension = 2;
  int solution_dimension = 2;
  int fact_threads = 4;

  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex1Boundary());
  boundary->initialize(num_boundary_points,
                       BoundaryCondition::DEFAULT);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), solution_dimension,
                           domain_dimension);
  std::vector<double> domain_points;
  get_domain_points(domain_size, &domain_points, quadtree.min,
                    quadtree.max, quadtree.min,
                    quadtree.max);
  Kernel kernel(solution_dimension, domain_dimension,
                pde, boundary.get(), domain_points);
  ki_Mat solution = boundary_integral_solve(kernel, *(boundary.get()),
                    &quadtree, id_tol, fact_threads, domain_points);

  QuadTree fresh;
  quadtree.copy_into(&fresh);
  ki_Mat new_sol = boundary_integral_solve(kernel, *(boundary.get()), &fresh,
                   id_tol,
                   fact_threads, domain_points);
  ASSERT_LE((new_sol - solution).vec_two_norm(), 1e-15);
}

}  // namespace kern_interp
