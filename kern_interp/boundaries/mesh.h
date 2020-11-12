// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_BOUNDARIES_MESH_H_
#define KERN_INTERP_BOUNDARIES_MESH_H_

#include <memory>
#include "kern_interp/boundaries/boundary.h"
#include "kern_interp/boundaries/cross_section_border.h"

namespace kern_interp {

class Mesh : public Boundary {
 public:
  double r = 1.0;

  CrossSectionBorder* cross_section_border;
  void create_cross_section_border();
  void initialize(int N, BoundaryCondition bc) override;
  bool is_in_domain(const PointVec& a) const override;
  std::unique_ptr<Boundary> clone() const override {
    return std::make_unique<Mesh>(*this);
  }
};

}  // namespace kern_interp

#endif  // KERN_INTERP_BOUNDARIES_MESH_H_
