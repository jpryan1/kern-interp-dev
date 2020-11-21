// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_BOUNDARIES_BOUNDARY_H_
#define KERN_INTERP_BOUNDARIES_BOUNDARY_H_

#include <vector>
#include <memory>
#include "kern_interp/pointvec.h"
#include "kern_interp/ki_mat.h"

#define STOKES_MIXER 0

namespace kern_interp {

struct Hole {
  PointVec center;
  double radius;
  int num_nodes;
};

enum BoundaryCondition {
  ELECTRON_3D,
  LAPLACE_PIPE_HOLES,
  LAPLACE_COW_MESH,
  SINGLE_ELECTRON,
  ALL_ONES,
  ALL_NEG_ONES,
  ALL_ZEROS,
  TANGENT_VEC,
  REVERSE_TANGENT_VEC,
  NORMAL_VEC,
  REVERSE_NORMAL_VEC,
  EX1,
  LEFT_TO_RIGHT_FLOW,
  NO_SLIP,
  HORIZONTAL_VEC,
  EX3A,
  EX3B,
  LAPLACE_CHECK_2D,
  LAPLACE_CHECK_3D,
  STOKES_SPHERES,
  STOKES_PIPE_HOLES,
  STOKES_COW_MESH,
  STOKES_2D_MIX,
  DEFAULT  // This is special, means use BC inherent to experiment/function.
};

class Boundary {
 public:
  int num_outer_nodes = -1;
  std::vector<double> perturbation_parameters;
  std::vector<double> points, normals, curvatures, weights;
  std::vector<Hole> holes;
  ki_Mat boundary_values;

  static bool comp_ang(PointVec a, PointVec b);
  virtual std::unique_ptr<Boundary> clone() const = 0;
  virtual ~Boundary() {}
  void set_boundary_values_size(BoundaryCondition bc);
  void apply_boundary_condition(int start_idx,
                                int end_idx,
                                BoundaryCondition bc);
  virtual void initialize(int n, BoundaryCondition bc) = 0;
  virtual bool is_in_domain(const PointVec& a) const = 0;
};


class CubicBoundary : public Boundary {
 public:
  // Note, the outer nodes must appear first in the point data vectors.

  virtual void get_spline_points(std::vector<double>* outer_x0_points,
                                 std::vector<double>* outer_x1_points) = 0;
  void get_cubics(const std::vector<double>& x0_points,
                  const std::vector<double>& x1_points,
                  std::vector<std::vector<double>>* x0_cubics,
                  std::vector<std::vector<double>>* x1_cubics);
  void interpolate(bool is_interior, int nodes_per_spline,
                   const std::vector<std::vector<double>>& x0_cubics,
                   const std::vector<std::vector<double>>& x1_cubics);
  void find_real_roots_of_cubic(const std::vector<double>& y_cubic,
                                std::vector<double>* t_vals);
  int num_right_intersections(double x, double y, int index);
  bool is_in_domain(const PointVec& a) const override;
};


}  // namespace kern_interp

#endif  // KERN_INTERP_BOUNDARIES_BOUNDARY_H_
