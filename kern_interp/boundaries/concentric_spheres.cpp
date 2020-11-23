// Copyright 2019 John Paul Ryan
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "kern_interp/boundaries/concentric_spheres.h"
#include "kern_interp/legendre.h"

namespace kern_interp {

void Sphere::initialize(int sz_param, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  string big_sphere, little_sphere;
  switch (sz_param) {
    case 0: {
      big_sphere = "sphere_19126_faces.txt";
      little_sphere = "sphere_570_faces.txt";
      break;
    } case 1: {
      big_sphere = "sphere_4192_faces.txt";
      little_sphere = "sphere_1032_faces.txt";
      break;
    } case 2: {
      big_sphere = "sphere_7970_faces.txt";
      little_sphere = "sphere_2344_faces.txt";
      break;
    } case 3: {
      big_sphere = "sphere_14870_faces.txt";
      little_sphere = "sphere_4192_faces.txt";
      break;
    } case 4: {
      big_sphere = "sphere_29214_faces.txt";
      little_sphere = "sphere_7970_faces.txt";
      break;
    } case 5: {
      big_sphere = "sphere_59744_faces.txt";
      little_sphere = "sphere_14870_faces.txt";
      break;
    } default: {
      big_sphere = "sphere_2344_faces.txt";
      little_sphere = "sphere_570_faces.txt";
      break;
    }
  }
  string data_dir = "kern_interp/boundaries/meshes/spheres/";
  std::vector<double> file_points, file_weights;
  string line;
  ifstream myfile(data_dir + big_sphere);
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
  } else {
    std::cout << "N0 File" << std::endl;
    exit(0);
  }
  num_outer_nodes = file_points.size() / 3;

  // Now we read in the data for the interior hole
  std::vector<double> file_hole_points, file_hole_weights;
  string hole_line;
  ifstream myholefile(data_dir + little_sphere);
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
  if (perturbation_parameters.size() == 0) {
    perturbation_parameters.push_back(0.0);
  }

  if (holes.size() == 0) {
    Hole hole;
    hole.center = PointVec(perturbation_parameters[0], perturbation_parameters[0],
                           perturbation_parameters[0]);
    hole.radius = 0.15;
    hole.num_nodes = num_hole_points;
    holes.push_back(hole);
  } else {
    holes[0].center = PointVec(perturbation_parameters[0],
                               perturbation_parameters[0],
                               perturbation_parameters[0]);
  }

  int num_holes = holes.size();
  int total_points = num_outer_nodes +
                     (num_holes * num_hole_points);
  for (int i = 0; i < num_outer_nodes; i++) {
    points.push_back(file_points[3 * i]);
    points.push_back(file_points[3 * i + 1]);
    points.push_back(file_points[3 * i + 2]);
    normals.push_back(file_points[3 * i]);
    normals.push_back(file_points[3 * i + 1]);
    normals.push_back(file_points[3 * i + 2]);
    weights.push_back(file_weights[i]);
  }

  // // For the interior hole, we dilate to 0.1x the original size.
  for (Hole hole : holes) {
    for (int i = 0; i < num_hole_points; i++) {
      // note the normals assume the unit norm of the file points
      points.push_back(hole.center.a[0]
                       + hole.radius * file_hole_points[3 * i]);
      points.push_back(hole.center.a[1]
                       + hole.radius * file_hole_points[3 * i + 1]);
      points.push_back(hole.center.a[2]
                       + hole.radius * file_hole_points[3 * i + 2]);
      normals.push_back(-file_hole_points[3 * i]);
      normals.push_back(-file_hole_points[3 * i + 1]);
      normals.push_back(-file_hole_points[3 * i + 2]);
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


void Sphere::create_cross_section_border() {
  // get faces near plane for outer
  if (cross_section_border != nullptr) {
    delete cross_section_border;
  }
  cross_section_border = new CrossSectionBorder();
  for (int i = 0; i < 100; i++) {
    double ang = (i / 100.0) * 2 * (M_PI);
    cross_section_border->sorted_2d_outer_knots.push_back(cos(ang));
    cross_section_border->sorted_2d_outer_knots.push_back(sin(ang));
  }

  Hole hole = holes[0];
  std::vector<double> hole_knots;
  for (int i = 0; i < 100; i++) {
    double ang = (i / 100.0) * 2 * (M_PI);
    hole_knots.push_back(hole.center.a[1] + hole.radius * cos(ang));
    hole_knots.push_back(hole.center.a[2] + hole.radius * sin(ang));
  }
  cross_section_border->sorted_2d_hole_knots.push_back(hole_knots);

  // populate and initialize cross section border
}


bool Sphere::is_in_domain(const PointVec& a) const {
  PointVec center(0.0, 0.0, 0.0);
  double eps = 2e-2;
  double dist = (center - a).norm();
  if (dist + eps > r) return false;
  for (Hole hole : holes) {
    dist = (hole.center - a).norm();
    if (dist - eps < hole.radius) return false;
  }
  return true;
}

}  // namespace kern_interp
