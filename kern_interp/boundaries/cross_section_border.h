// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_BOUNDARIES_CROSS_SECTION_BORDER_H_
#define KERN_INTERP_BOUNDARIES_CROSS_SECTION_BORDER_H_

#include <vector>
#include <memory>
#include "kern_interp/boundaries/boundary.h"

namespace kern_interp {

class CrossSectionBorder : public CubicBoundary {
 public:
  std::vector<double> sorted_2d_outer_knots;
  std::vector<std::vector<double>> sorted_2d_hole_knots;
  void initialize(int N, BoundaryCondition bc) override;

  void get_spline_points(std::vector<double>* outer_x0_spline_points,
                         std::vector<double>* outer_x1_spline_points) override;
  void get_inner_spline_points(
    std::vector<std::vector<double>>* x0_points,
    std::vector<std::vector<double>>* x1_points);

  std::unique_ptr<Boundary> clone() const override {
    return std::make_unique<CrossSectionBorder>(*this);
  }
};

}  // namespace kern_interp

#endif  // KERN_INTERP_BOUNDARIES_CROSS_SECTION_BORDER_H_
