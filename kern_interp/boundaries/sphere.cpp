// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "kern_interp/boundaries/sphere.h"
#include "kern_interp/legendre.h"

namespace kern_interp {

void Sphere::initialize(int sz_param, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  std::vector<double> file_points, file_weights;
  string line;
  ifstream myfile("kern_interp/boundaries/tri11146.txt");
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      stringstream s_stream(line);
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_points.push_back(std::stod(pt));
      }
      string wt;
      getline(s_stream, wt, ',');
      file_weights.push_back(std::stod(wt));
    }
    myfile.close();
  }
  num_outer_nodes = file_points.size() / 3;

  std::vector<double> file_hole_points, file_hole_weights;
  string hole_line;
  ifstream myholefile("kern_interp/boundaries/tri570.txt");
  if (myholefile.is_open()) {
    while (getline(myholefile, hole_line)) {
      stringstream s_stream(hole_line);
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_hole_points.push_back(std::stod(pt));
      }
      string wt;
      getline(s_stream, wt, ',');
      file_hole_weights.push_back(std::stod(wt));
    }
    myholefile.close();
  }


  int num_hole_points = file_hole_points.size() / 3;

  if (holes.size() == 0) {
    Hole hole;
    hole.center = PointVec(0.5, 0.5, 0.5);
    hole.radius = 0.1;
    hole.num_nodes = num_hole_points;
    holes.push_back(hole);
  }

  int num_holes = holes.size();

  int total_points = num_outer_nodes +
                     (num_holes * num_hole_points);

  for (int i = 0; i < num_outer_nodes; i++) {
    points.push_back(0.5 + file_points[3 * i]);
    points.push_back(0.5 + file_points[3 * i + 1]);
    points.push_back(0.5 + file_points[3 * i + 2]);
    normals.push_back(file_points[3 * i]);
    normals.push_back(file_points[3 * i + 1]);
    normals.push_back(file_points[3 * i + 2]);
    weights.push_back(file_weights[i]);
  }

  for (Hole hole : holes) {
    for (int i = 0; i < num_hole_points; i++) {
      points.push_back(hole.center.a[0] + hole.radius * file_hole_points[3 * i]);
      points.push_back(hole.center.a[1] + hole.radius * file_hole_points[3 * i + 1]);
      points.push_back(hole.center.a[2] + hole.radius * file_hole_points[3 * i + 2]);
      normals.push_back(file_hole_points[3 * i]);
      normals.push_back(file_hole_points[3 * i + 1]);
      normals.push_back(file_hole_points[3 * i + 2]);
      weights.push_back(hole.radius * hole.radius * file_hole_weights[i]);
    }
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ki_Mat(total_points, 1);
    apply_boundary_condition(0, total_points, ELECTRON_3D);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, total_points, bc);
  }
}

bool Sphere::is_in_domain(const PointVec& a) const {
  PointVec center(0.5, 0.5, 0.5);
  double eps = 0.1;
  double dist = (center - a).norm();
  if (dist + eps > r) return false;
  for (Hole hole : holes) {
    dist = (hole.center - a).norm();
    if (dist - eps < hole.radius) return false;
  }
  return true;
}

}  // namespace kern_interp
