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
#include "kern_interp/boundaries/ex3boundary.h"

#define THREE_B 1
#define DELTA_X 0.0001

namespace kern_interp {


void get_fin_diff_vals(double* samples, Boundary* boundary, int perturbed_param,
                       double id_tol, int fact_threads,
                       BoundaryCondition boundary_condition,
                       const QuadTree& quadtree, const Kernel& kernel,
                       const std::vector<double>& domain_points,
                       double* findiff) {
  QuadTree trees[4];
  std::unique_ptr<Boundary> boundaries[4];
  Kernel kernels[4];
  for (int i = 0; i < 4; i++) {
    quadtree.copy_into(&(trees[i]));
  }
  for (int i = 0; i < 4; i++) {
    boundary->perturbation_parameters[perturbed_param] = samples[i];
    boundary->initialize(boundary->weights.size(), boundary_condition);
    boundaries[i] = boundary->clone();
    kernels[i] = Kernel(kernel.solution_dimension, kernel.domain_dimension,
                        kernel.pde, boundaries[i].get(), domain_points);
    trees[i].perturb(*(boundaries[i].get()));
  }
  int outer_threads;
  if (fact_threads == 4) {
    outer_threads = 1;
  } else if (fact_threads == 2) {
    outer_threads = 2;
  } else {
    outer_threads = 4;
  }
  #pragma omp parallel for num_threads(outer_threads)
  for (int i = 0; i < 4; i++) {
    ki_Mat solution = boundary_integral_solve(kernels[i],
                      *(boundaries[i].get()), &trees[i], id_tol,
                      fact_threads, domain_points);
    double gradient;
    if (THREE_B) {
      gradient = -solution.get(0, 0);
    } else {
      gradient = (solution.get(1, 0) - solution.get(0, 0))
                 / (2 * DELTA_X);
    }
    findiff[i] = gradient;
  }
}


void enforce_separation(double* ang1, double* ang2) {
  double current_ang1 = *ang1;
  double current_ang2 = *ang2;

  while (current_ang1 > 2 * M_PI) current_ang1 -= 2 * M_PI;
  while (current_ang1 < 0) current_ang1 += 2 * M_PI;
  while (current_ang2 > 2 * M_PI) current_ang2 -= 2 * M_PI;
  while (current_ang2 < 0) current_ang2 += 2 * M_PI;

  double* lowerang;
  double* upperang;
  if (current_ang1 < current_ang2) {
    lowerang = &current_ang1;
    upperang = &current_ang2;
  } else {
    lowerang = &current_ang2;
    upperang = &current_ang1;
  }
  double dist = std::min(*upperang - *lowerang,
                         *lowerang + 2 * M_PI - *upperang);
  if (dist < M_PI / 4.) {
    double prob = (M_PI / 4.) - dist;
    if (*upperang - *lowerang < *lowerang + 2 * M_PI - *upperang) {
      *upperang += prob;
      *lowerang -= prob;
    } else {
      *upperang -= prob;
      *lowerang += prob;
    }
  }
  *ang1 = current_ang1;
  *ang2 = current_ang2;
}


void run_experiment3(int num_inner_threads) {
  double start_alpha = 1;
  double alpha_decay = 0.8;
  double h = 1e-4;
  double id_tol = 1e-6;
  int num_boundary_points = pow(2, 12);
  int domain_dimension = 2;
  int fact_threads = num_inner_threads;
  int FRAME_CAP = 30;
  BoundaryCondition boundary_condition;
  int solution_dimension;
  Kernel::Pde pde;
  if (THREE_B) {
    boundary_condition = BoundaryCondition::EX3B;
    solution_dimension = 2;
    pde = Kernel::Pde::STOKES;
  } else {
    boundary_condition = BoundaryCondition::EX3A;
    solution_dimension = 1;
    pde = Kernel::Pde::LAPLACE_NEUMANN;
  }
  std::unique_ptr<Boundary> boundary =
    std::unique_ptr<Boundary>(new Ex3Boundary());
  boundary->initialize(num_boundary_points, boundary_condition);

  double current_ang1 = -2.879, current_ang2 = 0.010;
  if (THREE_B) {
    current_ang1 = 1.4;
    current_ang2 = 3;
  }
  boundary->perturbation_parameters[0] = current_ang1;
  boundary->perturbation_parameters[1] = current_ang2;
  boundary->initialize(num_boundary_points, boundary_condition);
  QuadTree quadtree;
  quadtree.initialize_tree(boundary.get(), solution_dimension,
                           domain_dimension);

  std::vector<double> domain_points;
  domain_points.push_back(0.0 - DELTA_X);
  domain_points.push_back(0.0 + DELTA_X);
  if (!THREE_B) {
    domain_points.push_back(0.0 + DELTA_X);
    domain_points.push_back(0.0 + DELTA_X);
  }

  Kernel kernel(solution_dimension, domain_dimension,
                pde, boundary.get(), domain_points);
  ki_Mat solution = boundary_integral_solve(kernel,  *(boundary.get()),
                    &quadtree, id_tol, 4, domain_points);
  double prev_gradient;
  if (THREE_B) {
    prev_gradient = -solution.get(0, 0);
  } else {
    prev_gradient = (solution.get(1, 0) - solution.get(0, 0))
                    / (2 * DELTA_X);
  }


  for (int step = 0; step < FRAME_CAP; step++) {
    double starttime = omp_get_wtime();
    double findiff1[4];
    double samples1[4] = {current_ang1 - 2 * h, current_ang1 - h,
                          current_ang1 + h, current_ang1 + 2 * h
                         };
    double findiff2[4];
    double samples2[4] = {current_ang2 - 2 * h, current_ang2 - h,
                          current_ang2 + h, current_ang2 + 2 * h
                         };
    get_fin_diff_vals(samples1, boundary.get(), 0, id_tol, fact_threads,
                      boundary_condition, quadtree, kernel, domain_points,
                      findiff1);
    boundary->perturbation_parameters[0] = current_ang1;
    get_fin_diff_vals(samples2, boundary.get(), 1, id_tol, fact_threads,
                      boundary_condition, quadtree, kernel, domain_points,
                      findiff2);
    double grad1 = (findiff1[0] - 8 * findiff1[1] + 8 * findiff1[2]
                    - findiff1[3]) / (12 * h);
    double grad2 = (findiff2[0] - 8 * findiff2[1] + 8 * findiff2[2]
                    - findiff2[3]) / (12 * h);
    double alpha = start_alpha;
    double endtime = omp_get_wtime();
    std::cout << endtime - starttime << std::endl;
    while (alpha > 0.01) {
      double trial_ang1 = current_ang1 + alpha * grad1;
      double trial_ang2 = current_ang2 + alpha * grad2;
      enforce_separation(&trial_ang1, &trial_ang2);

      boundary->perturbation_parameters[0] = trial_ang1;
      boundary->perturbation_parameters[1] = trial_ang2;
      boundary->initialize(num_boundary_points, boundary_condition);
      quadtree.perturb(*boundary);
      kernel.update_data(boundary.get());
      ki_Mat solution = boundary_integral_solve(kernel, *(boundary.get()),
                        &quadtree, id_tol, fact_threads, domain_points);

      double gradient;
      if (THREE_B) {
        gradient = -solution.get(0, 0);
      } else {
        gradient = (solution.get(1, 0) - solution.get(0, 0))
                   / (2 * DELTA_X);
      }
      if (prev_gradient < gradient) {
        // std::cout << "theta1: " << trial_ang1 << " theta2: " << trial_ang2 <<
        //           " theta3: "
        //           << gradient << std::endl;
        prev_gradient = gradient;
        current_ang1 = trial_ang1;
        current_ang2 = trial_ang2;
        break;
      } else {
        alpha *= alpha_decay;
      }
    }
    if (alpha <= 0.01) {
      std::cout << "Line search did not terminate" << std::endl;
      exit(0);
    }
  }

  // std::ofstream sol_out;
  // sol_out.open("output/data/ie_solver_solution.txt");
  // int points_index = 0;
  // for (int i = 0; i < solution.height(); i += 2) {
  //   sol_out << domain_points[points_index] << "," <<
  //           domain_points[points_index + 1] << ",";
  //   points_index += 2;
  //   sol_out << solution.get(i, 0) << "," << solution.get(i + 1, 0) << std::endl;
  // }
  // sol_out.close();
  // std::ofstream bound_out;
  // bound_out.open("output/data/ie_solver_boundary.txt");
  // for (int i = 0; i < boundary->points.size(); i += 2) {
  //   bound_out << boundary->points[i] << "," << boundary->points[i + 1]
  //             << std::endl;
  // }
  // bound_out.close();
}

}  // namespace kern_interp


int main(int argc, char** argv) {
  srand(0);
  omp_set_nested(1);
  openblas_set_num_threads(1);
  std::cout << "inner_threads " << 1 << std::endl;
  for (int k = 0; k < 3; k++) {
    std::cout << "k " << k << std::endl;
    kern_interp::run_experiment3(1);
  }
  std::cout << "inner_threads " << 2 << std::endl;

  for (int k = 0; k < 3; k++) {
    std::cout << "k " << k << std::endl;

    kern_interp::run_experiment3(2);
  }
  std::cout << "inner_threads " << 4 << std::endl;

  for (int k = 0; k < 3; k++) {
    std::cout << "k " << k << std::endl;

    kern_interp::run_experiment3(4);
  }
  return 0;
}
