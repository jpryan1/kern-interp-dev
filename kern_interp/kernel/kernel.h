// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_KERNEL_KERNEL_H_
#define KERN_INTERP_KERNEL_KERNEL_H_

#include <vector>
#include <unordered_map>

#include "kern_interp/boundaries/boundary.h"
#include "kern_interp/ki_mat.h"
#include "kern_interp/legendre.h"
#include "kern_interp/quadtree/quadtree.h"

#define NUM_PROXY_POINTS 128
#define RADIUS_RATIO 1.5

namespace kern_interp {

struct Kernel {
  enum Pde {
    LAPLACE,
    LAPLACE_NEUMANN,
    GAUSS,
    STOKES
  };


  int solution_dimension, domain_dimension;
  Kernel::Pde pde;
  std::vector<double> boundary_points_, boundary_normals_,
      boundary_weights_, boundary_curvatures_, domain_points;

  std::vector<double> pxy_thetas, pxy_theta_weights, pxy_phis, pxy_phi_weights;
  std::vector<double> boundary_diags;
  std::vector<ki_Mat> boundary_diag_tensors;


  Kernel() {}
  Kernel(int solution_dimension_, int domain_dimension_,
         Kernel::Pde pde_, Boundary* boundary,
         std::vector<double> domain_points_) :
    solution_dimension(solution_dimension_),
    domain_dimension(domain_dimension_),
    pde(pde_),
    boundary_points_(boundary->points),
    boundary_normals_(boundary->normals),
    boundary_weights_(boundary->weights),
    boundary_curvatures_(boundary->curvatures),
    domain_points(domain_points_) {
    int num_circumf_points = 8;
    int num_phi_points = num_circumf_points / 2;
    double phis[num_phi_points];
    double phi_weights[num_phi_points];
    double phi_start = 0.;
    double phi_end = M_PI;
    cgqf(num_phi_points, 1, 0.0, 0.0, phi_start, phi_end, phis, phi_weights);
    for (int i = 0; i < num_circumf_points; i++) {
      double theta = 2 * M_PI * i * (1.0 / num_circumf_points);
      pxy_thetas.push_back(theta);
      pxy_theta_weights.push_back(2 * M_PI / (num_circumf_points));
    }
    for (int j = 0; j < num_phi_points; j++) {
      double phi = phis[j];
      pxy_phis.push_back(phi);
      pxy_phi_weights.push_back(phi_weights[j]*sin(phi));
    }
  }

  void compute_diag_entries_3dlaplace(Boundary* boundary);
  void compute_diag_entries_3dstokes(Boundary* boundary);

  void update_data(Boundary* boundary);

  void one_d_kern(int mat_idx, ki_Mat* ret, double r1,
                  double r2, double tn1, double tn2, double sn1, double sn2,
                  double sw, double sc, bool forward = false) const;
  void two_d_kern(int mat_idx, int tgt_parity, int src_parity,
                  ki_Mat* ret, double r1, double r2, double tn1, double tn2,
                  double sn1, double sn2, double sw, double sc,
                  bool forward = false) const;
  void three_d_laplace(int mat_idx, ki_Mat* ret, double r1, double r2, double r3,
                       double diag, double sn1, double sn2, double sn3,
                       double sw) const;
  void three_d_stokes(int mat_idx, int tgt_parity, int src_parity,
                      ki_Mat* ret, double r1, double r2, double r3, const ki_Mat& diag,
                      double tn1, double tn2, double tn3, double sn1, double sn2, double sn3,
                      double sw) const;

  ki_Mat operator()(const std::vector<int> & tgt_inds,
                    const std::vector<int> & src_inds,
                    bool forward = false, bool parallel = false) const;

  ki_Mat get_id_mat(const QuadTree * tree,
                    const QuadTreeNode * node) const;
  ki_Mat get_proxy_mat(std::vector<double> center, int num_points,
                       double r, const QuadTree * tree,
                       const std::vector<int> & box_inds) const;


  ki_Mat get_3d(const std::vector<int>& tgt_inds,
                const std::vector<int>& src_inds, bool forward, bool parallel) const;
  ki_Mat get_proxy_mat3d(std::vector<double> center,
                         double r, const QuadTree * tree,
                         const std::vector<int> & box_inds) const;

  ki_Mat forward() const;

  inline double boundary_normals(int i) const;
  inline double boundary_weights(int i) const;
  inline double boundary_curvatures(int i) const;
};  // struct


}  // namespace kern_interp

#endif  // KERN_INTERP_KERNEL_KERNEL_H_
