// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include "kern_interp/boundaries/donut.h"

namespace kern_interp {

void Donut::initialize(int N, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  if (holes.size() == 0) {
    Hole hole;
    hole.center = PointVec(0.0, 0.0);
    hole.radius = 0.5;
    holes.push_back(hole);
  }

  int hole_nodes = N / 5;
  int num_points = N + hole_nodes * holes.size();
  for (int i = 0; i < N; i++) {
    double ang = i * 2.0 * M_PI / N;
    double x = 0.0 + cos(ang);
    double y = 0.0 + sin(ang);
    points.push_back(x);
    points.push_back(y);
    normals.push_back(cos(ang));
    normals.push_back(sin(ang));
    curvatures.push_back(1);
    weights.push_back(2 * M_PI / N);
  }

  for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
    Hole hole = holes[hole_idx];
    int start_idx = N + hole_nodes * hole_idx;
    int end_idx = N + hole_nodes * (hole_idx + 1);
    for (int i = start_idx; i < end_idx; i++) {
      double ang = (i - start_idx) * 2.0 * M_PI / (end_idx - start_idx);
      double x = hole.center.a[0] + hole.radius * cos(ang);
      double y = hole.center.a[1] + hole.radius * sin(ang);
      points.push_back(x);
      points.push_back(y);
      normals.push_back(-cos(ang));
      normals.push_back(-sin(ang));
      curvatures.push_back(-1.0 / hole.radius);  // 1/r
      weights.push_back((2 * hole.radius * M_PI) / (end_idx - start_idx));
    }
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ki_Mat(num_points * 2, 1);
    for (int i = 0; i < N; i++) {
      boundary_values.set(2 * i, 0, -3 * normals[2 * i + 1]);
      boundary_values.set(2 * i + 1, 0, 3 * normals[2 * i]);
    }
    for (int i = N; i < num_points; i++) {
      boundary_values.set(2 * i, 0, normals[2 * i + 1]);
      boundary_values.set(2 * i + 1, 0, -normals[2 * i]);
    }
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, num_points, bc);
  }
}


bool Donut::is_in_domain(const PointVec& a) const {
  double x = a.a[0] - 0.0;
  double y = a.a[1] - 0.0;
  double eps = 1e-2;

  double dist = sqrt(pow(x, 2) + pow(y, 2));
  if (dist + eps > 1) return false;
  for (Hole hole : holes) {
    PointVec r = a - hole.center;
    if (r.norm() - eps < hole.radius) return false;
  }
  return true;
}

}  // namespace kern_interp
