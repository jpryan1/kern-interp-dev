// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include "kern_interp/ki_mat.h"
#include "kern_interp/quadtree/quadtree.h"
#include "kern_interp/linear_solve.h"


namespace kern_interp {


ki_Mat initialize_U_mat(const Kernel::Pde pde,
                        const std::vector<Hole>& holes,
                        const std::vector<double>& tgt_points, int domain_dimension) {
  ki_Mat U;
  if (holes.size() == 0) return ki_Mat(0, 0);
  // Get rid of magic numbers
  switch (pde) {
    case Kernel::Pde::LAPLACE: {
      if (domain_dimension == 3) {
        U = ki_Mat(tgt_points.size() / 3, holes.size());
        for (int i = 0; i < tgt_points.size(); i += 3) {
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            PointVec r = PointVec(tgt_points[i], tgt_points[i + 1],
                                  tgt_points[i + 2]) - hole.center;
            U.set(i / 3, hole_idx, 1.0 / (r.norm()));
          }
        }
      } else {
        U = ki_Mat(tgt_points.size() / 2, holes.size());
        for (int i = 0; i < tgt_points.size(); i += 2) {
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            PointVec r = PointVec(tgt_points[i], tgt_points[i + 1]) - hole.center;
            U.set(i / 2, hole_idx, log(r.norm()));
          }
        }
      }
      break;
    }
    case Kernel::Pde::STOKES: {
      if (domain_dimension == 3) {
        U = ki_Mat(tgt_points.size(), 6 * holes.size());
        for (int i = 0; i < tgt_points.size(); i += 3) {
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            PointVec r = PointVec(tgt_points[i], tgt_points[i + 1],
                                  tgt_points[i + 2]) - hole.center;

            for (int p = 0; p < 3; p++) {
              for (int s = 0; s < 3; s++) {
                int kron = (p == s ? 1 : 0);
                U.set(i + p, 6 * hole_idx + s, (1.0 / (8.0 * M_PI))
                      * ((kron / r.norm())
                         + (r.a[p]*r.a[s] / pow(r.norm(), 3))));
              }
            }
            //Cross product matrix
            //neg because it's gxr originally
            double rotscale = (-1.0 / (8.0 * M_PI * pow(r.norm(), 3)));
            U.set(i, 6 * hole_idx + 4, rotscale * (-r.a[2]));
            U.set(i, 6 * hole_idx + 5, rotscale * (r.a[1]));
            U.set(i + 1, 6 * hole_idx + 3, rotscale * (r.a[2]));
            U.set(i + 1, 6 * hole_idx + 5, rotscale * (-r.a[0]));
            U.set(i + 2, 6 * hole_idx + 3, rotscale * (-r.a[1]));
            U.set(i + 2, 6 * hole_idx + 4, rotscale * (r.a[0]));
          }
        }
      } else {
        U = ki_Mat(tgt_points.size(), 3 * holes.size());
        for (int i = 0; i < tgt_points.size(); i += 2) {
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            PointVec r = PointVec(tgt_points[i], tgt_points[i + 1]) - hole.center;
            double scale = 1.0 / (4 * M_PI);
            U.set(i, 3 * hole_idx, scale *
                  (log(1 / r.norm()) +
                   (1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[0]));
            U.set(i + 1, 3 * hole_idx, scale *
                  ((1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[0]));
            U.set(i, 3 * hole_idx + 1, scale *
                  ((1.0 / pow(r.norm(), 2)) * r.a[0] * r.a[1]));
            U.set(i + 1, 3 * hole_idx + 1, scale *
                  (log(1 / r.norm()) +
                   (1.0 / pow(r.norm(), 2)) * r.a[1] * r.a[1]));
            U.set(i, 3 * hole_idx + 2, r.a[1] * (scale / pow(r.norm(), 2)));
            U.set(i + 1, 3 * hole_idx + 2, -r.a[0] * (scale / pow(r.norm(), 2)));
          }
        }
      }
      break;
    }
    case Kernel::Pde::LAPLACE_NEUMANN: {
      U = ki_Mat(0, 0);
      break;
    }
  }
  return U;
}


ki_Mat initialize_Psi_mat(const Kernel::Pde pde,
                          const std::vector<Hole>& holes,
                          const Boundary& boundary, int domain_dimension) {
  if (holes.size() == 0) return ki_Mat(0, 0);
  ki_Mat Psi;
  switch (pde) {
    case Kernel::Pde::LAPLACE: {
      if (domain_dimension == 3) {
        Psi = ki_Mat(holes.size(), boundary.points.size() / 3);
        for (int i = 0; i < boundary.points.size(); i += 3) {
          PointVec x = PointVec(boundary.points[i], boundary.points[i + 1],
                                boundary.points[i + 2]);
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            if ((x - hole.center).norm() < hole.radius + 1e-8) {
              Psi.set(hole_idx, i / 3, boundary.weights[i / 3]);
              break;
            }
          }
        }
      } else {
        Psi = ki_Mat(holes.size(), boundary.points.size() / 2);
        for (int i = 0; i < boundary.points.size(); i += 2) {
          PointVec x = PointVec(boundary.points[i], boundary.points[i + 1]);
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            if ((x - hole.center).norm() < hole.radius + 1e-8) {
              Psi.set(hole_idx, i / 2, boundary.weights[i / 2]);
              break;
            }
          }
        }
      }
      break;
    }
    case Kernel::Pde::STOKES: {
      if (domain_dimension == 3) {

        Psi = ki_Mat(6 * holes.size(), boundary.points.size());
        for (int i = 0; i < boundary.points.size(); i += 3) {
          PointVec x = PointVec(boundary.points[i], boundary.points[i + 1],
                                boundary.points[i + 2]);

          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            PointVec r = x - hole.center;
            if (r.norm() < hole.radius + 1e-8) {

              Psi.set(6 * hole_idx, i, boundary.weights[i / 3]);
              Psi.set(6 * hole_idx + 1, i + 1, boundary.weights[i / 3]);
              Psi.set(6 * hole_idx + 2, i + 2, boundary.weights[i / 3]);

              //negative because transpose
              Psi.set(6 * hole_idx + 4, i, -boundary.weights[i / 3] * (-r.a[2]));
              Psi.set(6 * hole_idx + 5, i, -boundary.weights[i / 3] * (r.a[1]));
              Psi.set(6 * hole_idx + 3, i + 1, -boundary.weights[i / 3] * (r.a[2]));
              Psi.set(6 * hole_idx + 5, i + 1, -boundary.weights[i / 3] * (-r.a[0]));
              Psi.set(6 * hole_idx + 3, i + 2, -boundary.weights[i / 3] * (-r.a[1]));
              Psi.set(6 * hole_idx + 4, i + 2, -boundary.weights[i / 3] * (r.a[0]));

              break;
            }
          }
        }

      } else {
        Psi = ki_Mat(3 * holes.size(), boundary.points.size());
        for (int i = 0; i < boundary.points.size(); i += 2) {
          PointVec x = PointVec(boundary.points[i], boundary.points[i + 1]);
          for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
            Hole hole = holes[hole_idx];
            if ((x - hole.center).norm() < hole.radius + 1e-8) {
              Psi.set(3 * hole_idx, i, boundary.weights[i / 2]);
              Psi.set(3 * hole_idx + 1, i + 1, boundary.weights[i / 2]);
              Psi.set(3 * hole_idx + 2, i, boundary.weights[i / 2]*x.a[1]);
              Psi.set(3 * hole_idx + 2, i + 1, -boundary.weights[i / 2]*x.a[0]);
              break;
            }
          }
        }
      }
      break;
    }
    case Kernel::Pde::LAPLACE_NEUMANN: {
      Psi = ki_Mat(0, 0);
      break;
    }
  }
  return Psi;
}


void linear_solve(const SkelFactorization& skel_factorization,
                  const QuadTree& quadtree, const ki_Mat& f, ki_Mat* mu,
                  ki_Mat* alpha) {
  *mu = ki_Mat(f.height(), 1);
  if (alpha == nullptr) {
    skel_factorization.solve(quadtree, mu, f);
  } else {
    *alpha = ki_Mat(quadtree.U.width(), 1);
    skel_factorization.multiply_connected_solve(quadtree, mu, alpha, f);
  }
}


void schur_solve(const SkelFactorization & skel_factorization,
                 const QuadTree & quadtree, const ki_Mat & U,
                 const ki_Mat & Psi,
                 const ki_Mat & f, const ki_Mat & K_domain,
                 const ki_Mat & U_forward,  ki_Mat * solution) {
  ki_Mat mu(f.height(), 1);
  if (U.width() == 0) {
    double start = omp_get_wtime();
    linear_solve(skel_factorization, quadtree, f, &mu);
    double end = omp_get_wtime();
    std::cout << "solve time " << end - start << std::endl;
    *solution = K_domain * mu;
  } else {
    ki_Mat alpha;
    double start = omp_get_wtime();
    linear_solve(skel_factorization, quadtree, f, &mu, &alpha);
    double end = omp_get_wtime();
    std::cout << "solve time " << end - start << std::endl;
    *solution = (K_domain * mu) + (U_forward * alpha);
    // std::cout<<"Alpha "<<alpha.get(0,0)<<" "<<alpha.get(1,0)<<" "<<alpha.get(2,0)<<std::endl;
  }
}


ki_Mat boundary_integral_solve(const Kernel& kernel, const Boundary& boundary,
                               QuadTree * quadtree, double id_tol,
                               int fact_threads,
                               const std::vector<double>& domain_points) {
  // Consider making init instead of constructor for readability
  SkelFactorization skel_factorization(id_tol, fact_threads);
  ki_Mat K_domain = kernel.forward();
  ki_Mat f = boundary.boundary_values;
  ki_Mat U = initialize_U_mat(kernel.pde, boundary.holes, boundary.points,
                              kernel.domain_dimension);
  ki_Mat Psi = initialize_Psi_mat(kernel.pde, boundary.holes,
                                  boundary, kernel.domain_dimension);
  ki_Mat U_forward = initialize_U_mat(kernel.pde, boundary.holes,
                                      kernel.domain_points, kernel.domain_dimension);

  ki_Mat domain_solution((kernel.domain_points.size() / kernel.domain_dimension)*
                         kernel.solution_dimension, 1);
  quadtree->U = U;
  quadtree->Psi = Psi;

  double start = omp_get_wtime();
  skel_factorization.skeletonize(kernel, quadtree);
  double end = omp_get_wtime();
  std::cout  << "skel time " << end - start << std::endl;

  schur_solve(skel_factorization, *quadtree, U, Psi, f, K_domain,
              U_forward, &domain_solution);


  for (int i = 0; i < kernel.domain_points.size(); i += kernel.domain_dimension) {
    std::vector<double> vec;
    for (int j = 0; j < kernel.domain_dimension; j++) {
      vec.push_back(kernel.domain_points[i + j]);
    }
    PointVec point(vec);

    if (!boundary.is_in_domain(point)) {
      int pt_idx = i / kernel.domain_dimension;
      for (int j = 0; j < kernel.solution_dimension; j++) {
        domain_solution.set(kernel.solution_dimension * pt_idx + j, 0, 0.);
      }
    }
  }

  return domain_solution;
}


void get_domain_points(int domain_size, std::vector<double>* points,
                       double x_min, double x_max, double y_min, double y_max) {
  for (int i = 0; i < domain_size; i++) {
    double x = x_min + ((i + 0.0) / (domain_size - 1)) * (x_max - x_min);
    for (int j = 0; j < domain_size; j++) {
      double y = y_min + ((j + 0.0) / (domain_size - 1)) * (y_max - y_min);
      points->push_back(x);
      points->push_back(y);
    }
  }
}

void get_domain_points3d(int domain_size, std::vector<double>* points,
                         Boundary* boundary, double min, double max) {
  double eps = 0.01;


  for (int i = 0; i < domain_size; i++) {
    double x = -0.5 + ((i + 0.0) / (domain_size - 1)) * (2);
    for (int j = 0; j < domain_size; j++) {
      double y = -0.5 + ((j + 0.0) / (domain_size - 1)) * (2);
      for (int k = 0; k < domain_size; k++) {
        double z = -0.5 + ((k + 0.0) / (domain_size - 1)) * (2);
        if (!boundary->is_in_domain(PointVec(x, y, z))) {
          continue;
        }
        points->push_back(x);
        points->push_back(y);
        points->push_back(z);
      }
    }
  }



  // for (int i = 0; i < domain_size; i++) {
  //   double r = eps + min + ((max - min - eps) * (i / (domain_size + 0.)));
  //   // double x = min + ((i + 0.0) / (domain_size - 1)) * (max - min);
  //   for (int j = 0; j < domain_size; j++) {
  //     double theta = 2 * M_PI * (j / (domain_size + 0.));
  //     // double y = min + ((j + 0.0) / (domain_size - 1)) * (max - min);
  //     for (int k = 0; k < domain_size; k++) {
  //       double phi = M_PI * (k / (domain_size + 0.));
  //       // double z = min + ((k + 0.0) / (domain_size - 1)) * (max - min);
  //       // double phi=0;
  //       // double theta=0;
  //       double x = 0.5 + r * sin(phi) * cos(theta);
  //       double y = 0.5 + r * sin(phi) * sin(theta);
  //       double z = 0.5 + r * cos(phi);
  //       points->push_back(x);
  //       points->push_back(y);
  //       points->push_back(z);
  //     }
  //   }
  // }
}


double solve_err(const Kernel& kernel, Boundary* boundary, double id_tol) {
  int num_threads = 4;

  QuadTree quadtree;
  quadtree.initialize_tree(boundary, kernel.solution_dimension,
                           kernel.domain_dimension);

  SkelFactorization skel_factorization(id_tol, num_threads);

  double max = 0.;
  ki_Mat dense;

  if (boundary->holes.size() > 0
      && kernel.pde != Kernel::Pde::LAPLACE_NEUMANN) {
    ki_Mat U = initialize_U_mat(kernel.pde, boundary->holes, boundary->points,
                                kernel.domain_dimension);
    ki_Mat Psi = initialize_Psi_mat(kernel.pde, boundary->holes, *boundary,
                                    kernel.domain_dimension);
    quadtree.U = U;
    quadtree.Psi = Psi;
    skel_factorization.skeletonize(kernel, &quadtree);

    ki_Mat mu, alpha;
    linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu,
                 &alpha);
    ki_Mat stacked(mu.height() + alpha.height(), 1);
    stacked.set_submatrix(0, mu.height(), 0, 1, mu);
    stacked.set_submatrix(mu.height(), stacked.height(), 0, 1, alpha);

    std::vector<int> all_dofs;
    for (int i = 0; i < kernel.solution_dimension * boundary->weights.size();
         i++) {
      all_dofs.push_back(i);
    }
    ki_Mat kern = kernel(all_dofs, all_dofs);
    int hole_factor = kernel.solution_dimension == 2 ? 3 : 1;
    int added = hole_factor * boundary->holes.size();
    dense = ki_Mat(all_dofs.size() + added, all_dofs.size() + added);
    ki_Mat ident(added, added);
    ident.eye(added);
    dense.set_submatrix(0, all_dofs.size(), 0, all_dofs.size(), kern);
    dense.set_submatrix(0, all_dofs.size(), all_dofs.size(), dense.width(), U);
    dense.set_submatrix(all_dofs.size(), dense.height(),
                        0, all_dofs.size(), Psi);
    dense.set_submatrix(all_dofs.size(), dense.height(),
                        all_dofs.size(), dense.width(), -ident);

    ki_Mat fzero_prime = dense * stacked;

    ki_Mat err1 = (fzero_prime(0, mu.height(), 0, 1)
                   - boundary->boundary_values);
    ki_Mat err2 = (fzero_prime(mu.height(), fzero_prime.height(), 0, 1));

    double err = sqrt(pow(err1.vec_two_norm(), 2) + pow(err1.vec_two_norm(), 2));

    return err / boundary->boundary_values.vec_two_norm();
  } else {
    skel_factorization.skeletonize(kernel, &quadtree);
    ki_Mat mu;
    linear_solve(skel_factorization, quadtree, boundary->boundary_values, &mu);
    std::vector<int> all_dofs;
    for (int i = 0; i < kernel.solution_dimension * boundary->weights.size();
         i++) {
      all_dofs.push_back(i);
    }

    double start = omp_get_wtime();
    ki_Mat bigK = kernel(all_dofs, all_dofs);
    double end = omp_get_wtime();
    std::cout << "Big kern call took " << (end - start) << std::endl;
    ki_Mat err = (bigK * mu) - boundary->boundary_values;
    return err.vec_two_norm() / boundary->boundary_values.vec_two_norm();
  }
}



ki_Mat stokes_true_sol(const std::vector<double>& domain_points,
                       Boundary * boundary, double c1, double c2) {
  double diff_norm;
  double norm_of_true;

  ki_Mat truth(domain_points.size(), 1);

  for (int i = 0; i < domain_points.size(); i += 2) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    PointVec x(x0, x1);
    if (!boundary->is_in_domain(x)) {
      truth.set(i, 0, 0.0);
      truth.set(i + 1, 0, 0.0);
      continue;
    }

    PointVec center(0.5, 0.5);
    PointVec r = x - center;
    PointVec true_vec = PointVec(-r.a[1], r.a[0]);
    switch (boundary->holes.size()) {
      case 0:  // circle
        true_vec = true_vec * (r.norm() / true_vec.norm()) * 4.;
        break;
      case 1:  // donut
      default:
        double p = (c1 / r.norm()) + c2 * r.norm();
        true_vec = true_vec * (1. / true_vec.norm()) * p;
        break;
    }
    truth.set(i, 0, true_vec.a[0]);
    truth.set(i + 1, 0, true_vec.a[1]);
  }

  return truth;
}


double stokes_err_3d(const ki_Mat& domain,
                     const std::vector<double>& domain_points,
                     Boundary* boundary, double rsmall, double w) {


  ki_Mat abcd(4, 4);
  abcd.set(0, 0, 1.);
  abcd.set(0, 1, 1.);
  abcd.set(0, 2, 1.);
  abcd.set(0, 3, 1.);

  abcd.set(1, 0, 1.);
  abcd.set(1, 1, -1.);
  abcd.set(1, 2, -2.);
  abcd.set(1, 3, -4.);

  abcd.set(2, 0, 1.0 / pow(rsmall, 3));
  abcd.set(2, 1, 1.0 / rsmall);
  abcd.set(2, 2, 1.0);
  abcd.set(2, 3, pow(rsmall, 2));

  abcd.set(3, 0, 1.0 / pow(rsmall, 3));
  abcd.set(3, 1, -1.0 / rsmall);
  abcd.set(3, 2, -2);
  abcd.set(3, 3, -4 * pow(rsmall, 2));

  ki_Mat b(4, 1);
  b.set(0, 0, 0.5);
  b.set(1, 0, -1.0);
  b.set(2, 0, w / 2.0);
  b.set(3, 0, w);

  ki_Mat x_vals(4, 1);
  abcd.left_multiply_inverse(b, &x_vals);

  double err = 0.;
  double tot = 0.;
  for (int i = 0; i < domain_points.size(); i += 3) {
    double x0 = domain_points[i];
    double x1 = domain_points[i + 1];
    double x2 = domain_points[i + 2];

    PointVec x(x0, x1, x2);
    if (!boundary->is_in_domain(x)) {
      continue;
    }

    PointVec center(0.5, 0.5, 0.5);
    PointVec point = x - center;

    double r, phi;
    r = point.norm();
    phi = acos(point.a[2] / point.norm());

    double ur, uphi;

    ur = 2 * cos(phi) * (
           (x_vals.get(0, 0) / pow(r, 3))
           + (x_vals.get(1, 0) / r)
           + (x_vals.get(2, 0))
           + (x_vals.get(3, 0) * pow(r, 2)));
    uphi = sin(phi) * (
             (x_vals.get(0, 0) / pow(r, 3))
             + (-x_vals.get(1, 0) / r)
             + (-2 * x_vals.get(2, 0))
             + (-4 * x_vals.get(3, 0) * pow(r, 2)));

    PointVec guessvel = PointVec(domain.get(i, 0), domain.get(i + 1, 0),
                                 domain.get(i + 2, 0));
    double guessmag = guessvel.norm();
    double predmag = sqrt(pow(ur, 2) + pow(uphi, 2));
    err += pow(guessmag - predmag, 2);
    tot += pow(predmag, 2);
  }
  return sqrt(err) / sqrt(tot);
}


}  // namespace kern_interp
