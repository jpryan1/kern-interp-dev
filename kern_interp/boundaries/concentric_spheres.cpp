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
// 100 = 6 + x + 29x
  string big_sphere, little_sphere;
  switch (sz_param) {
    case 0: {
      big_sphere = "sphere_7680_faces.txt";
      little_sphere = "sphere_240_faces.txt";
      // big_sphere = "sphere_29214_faces.txt";
      // little_sphere = "sphere_1032_faces.txt";
      break;
    } case 1: {
      big_sphere = "sphere_20496_faces.txt";
      little_sphere = "sphere_648_faces.txt";
      break;
    } case 2: {
      big_sphere = "sphere_47632_faces.txt";
      little_sphere = "sphere_1488_faces.txt";
      break;
    } case 3: {
      big_sphere = "sphere_101070_faces.txt";
      little_sphere = "sphere_3206_faces.txt";
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


  std::vector<double> file_hole_cow_points, file_hole_cow_normals,
      file_hole_cow_weights;
  ifstream mycowfile("kern_interp/boundaries/meshes/cow_mesh.txt");
  string cow_line;
  if (mycowfile.is_open()) {
    while (getline(mycowfile, cow_line)) {
      stringstream s_stream(cow_line);
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_hole_cow_points.push_back(std::stod(pt));
      }
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_hole_cow_normals.push_back(std::stod(pt));
      }
      string wt;
      getline(s_stream, wt, ',');
      file_hole_cow_weights.push_back(std::stod(wt));
    }
    mycowfile.close();
  }
  int num_cow_points = file_hole_cow_points.size() / 3;

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
    hole.radius = 0.1;
    hole.num_nodes = num_hole_points;

    hole.center = PointVec(0, -0.45, 0.1);

    holes.push_back(hole);

    hole.center = PointVec(0, 0.3, 0);
    hole.num_nodes = num_cow_points;
    hole.radius = 0.4;
    holes.push_back(hole);

  } else {
    holes[0].center = PointVec(0, -0.45, perturbation_parameters[0]);
  }

  int num_holes = holes.size();
  int total_points = num_outer_nodes +
                     (num_hole_points + num_cow_points);
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
  Hole hole = holes[0];
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
  hole = holes[1];
  for (int i = 0; i < hole.num_nodes; i++) {
    double x = file_hole_cow_points[3 * i];
    double y =  file_hole_cow_points[3 * i + 1];
    double z = file_hole_cow_points[3 * i + 2];
    double nx = file_hole_cow_normals[3 * i];
    double ny =  file_hole_cow_normals[3 * i + 1] ;
    double nz = file_hole_cow_normals[3 * i + 2];
    points.push_back(hole.center.a[0] + hole.radius * x);
    points.push_back(hole.center.a[1] +  hole.radius * y);
    points.push_back(hole.center.a[2] +  hole.radius * z);
    normals.push_back(-nx);
    normals.push_back(-ny);
    normals.push_back(-nz);
    weights.push_back(hole.radius * hole.radius * file_hole_cow_weights[i]);
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


void Sphere::order_mesh_nodes(std::vector<double>* points, double cx,
                              double cy) {
  std::vector<double> ordered;
  for (int i = 0; i < 2; i++) {
    ordered.push_back((*points)[i]);
  }
  double desired_size = points->size();
  while (ordered.size() != desired_size) {
    double curr_x = ordered[ordered.size() - 2];
    double curr_y = ordered[ordered.size() - 1];

    double closest[2];
    double closest_dist = 100;
    int point_idx = 0;
    for (int i = 0; i < points->size(); i += 2) {
      double point_x = (*points)[i];
      double point_y = (*points)[i + 1];
      if (ordered.size() > 2) {
        if (abs(point_x - ordered[ordered.size() - 4]) < 1e-12
            && abs(point_y - ordered[ordered.size() - 3]) < 1e-12) {
          continue;
        }
      }
      double dist = sqrt(pow(point_x - curr_x, 2) + pow(point_y - curr_y, 2));
      if (abs(dist) < 1e-12) {
        continue;
      }
      if (dist < closest_dist) {
        closest_dist = dist;
        closest[0] = point_x;
        closest[1] = point_y;
        point_idx = i;
      }
    }
    points->erase(points->begin() + point_idx, points->begin() + point_idx + 2);
    ordered.push_back(closest[0]);
    ordered.push_back(closest[1]);
  }
  double tot_ang = 0;
  for (int i = 0; i < ordered.size() - 2; i += 2) {
    tot_ang += ((ordered[i + 1] - cy) - (ordered[i] - cx)) * ((
                 ordered[i + 2] - cx) + (ordered[i + 3] - cy));
  }
  tot_ang += (ordered[ordered.size() - 1] - ordered[ordered.size() - 2]) *
             (ordered[0] + ordered[1]);

  if (tot_ang > 0) {
    points->clear();
    for (int i = ordered.size() - 2; i >= 0; i -= 2) {
      points->push_back(ordered[i]);
      points->push_back(ordered[i + 1]);
    }
  } else {
    points->clear();

    for (int i = 0; i < ordered.size(); i += 2) {

      points->push_back(ordered[i]);
      points->push_back(ordered[i + 1]);
    }
  }
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
  // for (Hole hole : holes) {
  Hole hole = holes[0];
  std::vector<double> hole_knots;
  for (int i = 0; i < 100; i++) {
    double ang = (i / 100.0) * 2 * (M_PI);
    hole_knots.push_back(hole.center.a[1] + hole.radius * cos(ang));
    hole_knots.push_back(hole.center.a[2] + hole.radius * sin(ang));
    // }
  }
  cross_section_border->sorted_2d_hole_knots.push_back(hole_knots);
  hole_knots.clear();

  // populate and initialize cross section border
  string hole_cow_filename = "kern_interp/boundaries/meshes/cow_mesh_cross.txt";
  string line;
  hole = holes[1];
  ifstream myholefile;
  // if (i == 0) {
  myholefile = ifstream(hole_cow_filename);
  // } else {
  //   myholefile = ifstream(hole_ellipsoid_filename);
  // }
  if (myholefile.is_open()) {
    while (getline(myholefile, line)) {
      stringstream s_stream(line);
      std::vector<double> pt_vec;
      for (int d = 0; d < 2; d++) {
        string pt;
        getline(s_stream, pt, ',');
        pt_vec.push_back(std::stod(pt));
      }
      PointVec rotated(pt_vec[0], pt_vec[1]);
      hole_knots.push_back(rotated.a[0]*hole.radius);
      hole_knots.push_back(rotated.a[1]*hole.radius);
    }
    myholefile.close();
  } else {
    std::cout << "N0 File" << std::endl;
    exit(0);
  }
  // the args are called x and y although we are passing y and z - bad code
  order_mesh_nodes(&hole_knots, hole.center.a[1], hole.center.a[2]);
  std::vector<double> sorted_hole_knots;
  for (int j = 0; j < hole_knots.size(); j += 2) {
    sorted_hole_knots.push_back(hole.center.a[1] + hole_knots[j]);
    sorted_hole_knots.push_back(hole.center.a[2] + hole_knots[j + 1]);
  }
  cross_section_border->sorted_2d_hole_knots.push_back(sorted_hole_knots);
}


bool Sphere::is_in_domain(const PointVec& a) const {
  return cross_section_border->is_in_domain(PointVec(a.a[1], a.a[2]));

  // PointVec center(0.0, 0.0, 0.0);
  // double eps = 2e-2;
  // double dist = (center - a).norm();
  // if (dist + eps > r) return false;
  // for (Hole hole : holes) {
  //   dist = (hole.center - a).norm();
  //   if (dist - eps < hole.radius) return false;
  // }
  // return true;
}

}  // namespace kern_interp
