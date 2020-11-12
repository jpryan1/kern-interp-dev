// Copyright 2019 John Paul Ryan
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "kern_interp/boundaries/mesh.h"
#include "kern_interp/legendre.h"

namespace kern_interp {

void Mesh::initialize(int sz_param, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  std::vector<double> file_points, file_normals, file_weights;
  string line, outer_filename;
  // ifstream myfile("kern_interp/boundaries/meshes/spheres/sphere_4192_faces.txt");
  if (sz_param > 0) {
    outer_filename = "kern_interp/boundaries/meshes/cow_mesh.txt";

  } else {
    outer_filename = "kern_interp/boundaries/meshes/pipe_mesh.txt";
  }

  ifstream myfile(outer_filename);
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      stringstream s_stream(line);
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_points.push_back(std::stod(pt));
      }
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_normals.push_back(std::stod(pt));
      }
      string wt;
      getline(s_stream, wt, ',');
      file_weights.push_back(std::stod(wt));
    }
    myfile.close();
  } else {
    std::cout << "N0 File" << std::endl;
    exit(0);
  }
  num_outer_nodes = file_points.size() / 3;

  // Now we read in the data for the interior hole
  std::vector<double> file_hole_points, file_hole_normals, file_hole_weights;
  string hole_line;
  ifstream
  myholefile("kern_interp/boundaries/meshes/ellipsoid.txt");
  if (myholefile.is_open()) {
    while (getline(myholefile, hole_line)) {
      stringstream s_stream(hole_line);
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_hole_points.push_back(std::stod(pt));
      }
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_hole_normals.push_back(std::stod(pt));
      }
      string wt;
      getline(s_stream, wt, ',');
      file_hole_weights.push_back(std::stod(wt));
    }
    myholefile.close();
  }
  int num_hole_points = file_hole_points.size() / 3;

  if (perturbation_parameters.size() == 0) {
    perturbation_parameters.push_back(0.0);
  }

  if (holes.size() == 0) {
    Hole hole;
    hole.center = PointVec(0.0, 0.0, 5.0);
    hole.radius = 0.5*sqrt(2);
    hole.num_nodes = num_hole_points;
    holes.push_back(hole);
    hole.center = PointVec(0.0, 0.0, 2.0);
    holes.push_back(hole);
    hole.center = PointVec(0.0, 0.0, 8.0);
    holes.push_back(hole);
  }

  int num_holes = holes.size();

  int total_points = num_outer_nodes +
                     (num_holes * num_hole_points);
  for (int i = 0; i < num_outer_nodes; i++) {
    points.push_back(file_points[3 * i]);
    points.push_back(file_points[3 * i + 1]);
    points.push_back(file_points[3 * i + 2]);
    normals.push_back(file_normals[3 * i ]);
    normals.push_back(file_normals[3 * i + 1]);
    normals.push_back(file_normals[3 * i + 2]);
    weights.push_back(file_weights[i]);
  }
  // For the interior hole, we dilate to 0.1x the original size.
  for (int hole_idx = 0; hole_idx < holes.size(); hole_idx++) {
    Hole hole = holes[hole_idx];
    double ang;
    if (hole_idx == 1) {
      ang = M_PI / 3.0;
    } else if (hole_idx == 2) {
      ang = -M_PI / 3.0;
    } else {
      ang = 0;
    }
    for (int i = 0; i < num_hole_points; i++) {
      double x = file_hole_points[3 * i];
      double y = cos(ang) * file_hole_points[3 * i + 1] + sin(
                   ang) * file_hole_points[3 * i + 2];
      double z = -sin(ang) * file_hole_points[3 * i + 1] + cos(
                   ang) * file_hole_points[3 * i + 2];
      double nx = file_hole_normals[3 * i];
      double ny = cos(ang) * file_hole_normals[3 * i + 1] + sin(
                    ang) * file_hole_normals[3 * i + 2];
      double nz = -sin(ang) * file_hole_normals[3 * i + 1] + cos(
                    ang) * file_hole_normals[3 * i + 2];
      points.push_back(hole.center.a[0]
                       + hole.radius * x);
      points.push_back(hole.center.a[1]
                       +  hole.radius * y);
      points.push_back(hole.center.a[2]
                       +  hole.radius * z);
      normals.push_back(nx);
      normals.push_back(ny);
      normals.push_back(nz);
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

  create_cross_section_border();
  cross_section_border->initialize(sz_param, bc);
}


void Mesh::create_cross_section_border() {
  // get faces near plane for outer
  cross_section_border = new CrossSectionBorder();
  double height_thresh = 0.1;
  std::vector<PointVec> outer_2d_points;
  for (int i = 0; i < num_outer_nodes; i++) {
    if (abs(points[3 * i]) < height_thresh) {
      outer_2d_points.push_back(PointVec(points[3 * i + 1], points[3 * i + 2]));
    }
  }

  // sort points by angle
  std::sort(outer_2d_points.begin(), outer_2d_points.end(), comp_ang);
  for (int i = 0; i < outer_2d_points.size(); i++) {
    PointVec pv = outer_2d_points[i];
    cross_section_border->sorted_2d_outer_knots.push_back(pv.a[0]);
    cross_section_border->sorted_2d_outer_knots.push_back(pv.a[1]);
  }
  // same with holes, sort by angle about center
  double point_idx = num_outer_nodes;
  for (int i = 0; i < holes.size(); i++) {
    Hole hole = holes[i];
    std::vector<PointVec> hole_knots;
    for (int j = 0; j < hole.num_nodes; j++) {
      PointVec holept3d(points[3 * (point_idx + j)], points[3 * (point_idx + j) + 1],
                        points[3 * (point_idx + j) + 2]);
      if (abs(holept3d.a[0]) < height_thresh) {
        hole_knots.push_back(PointVec(holept3d.a[1], holept3d.a[2]));
      }
    }
    std::sort(hole_knots.begin(), hole_knots.end(), comp_ang);
    std::vector<double> sorted_hole_knots;
    for (int j = 0; j < hole_knots.size(); j++) {
      PointVec pv = hole_knots[j];
      sorted_hole_knots.push_back(pv.a[0]);
      sorted_hole_knots.push_back(pv.a[1]);
    }
    cross_section_border->sorted_2d_hole_knots.push_back(sorted_hole_knots);
    point_idx += hole.num_nodes;
  }

  // populate and initialize cross section border
}


bool Mesh::is_in_domain(const PointVec& a) const {
  return cross_section_border->is_in_domain(PointVec(a.a[1], a.a[2]));
}

}  // namespace kern_interp
