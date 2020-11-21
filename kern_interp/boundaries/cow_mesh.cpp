// Copyright 2019 John Paul Ryan
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "kern_interp/boundaries/cow_mesh.h"
#include "kern_interp/legendre.h"

#define CIRCLE_RAD 1

namespace kern_interp {

void CowMesh::initialize(int sz_param, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();

  std::vector<double> file_points, file_normals, file_weights;
  string line, outer_filename;
  outer_filename = "kern_interp/boundaries/meshes/spheres/sphere_9496_faces.txt";


  ifstream myfile(outer_filename);
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      stringstream s_stream(line);
      for (int d = 0; d < 3; d++) {
        string pt;
        getline(s_stream, pt, ',');
        file_points.push_back(std::stod(pt));
      }
      // for (int d = 0; d < 3; d++) {
      //   string pt;
      //   getline(s_stream, pt, ',');
      //   file_normals.push_back(std::stod(pt));
      // }
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
  std::vector<double> file_hole_cow_points, file_hole_cow_normals,
      file_hole_cow_weights;
  string hole_line;
  ifstream myholefile("kern_interp/boundaries/meshes/cow_mesh.txt");
  if (myholefile.is_open()) {
    while (getline(myholefile, hole_line)) {
      stringstream s_stream(hole_line);
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
    myholefile.close();
  }
  int num_hole_cow_points = file_hole_cow_points.size() / 3;

  // std::vector<double> file_hole_ellipsoid_points, file_hole_ellipsoid_normals,
  //     file_hole_ellipsoid_weights;
  // string ellipsoid_hole_line;
  // myholefile = ifstream("kern_interp/boundaries/meshes/ellipsoid.txt");
  // if (myholefile.is_open()) {
  //   while (getline(myholefile, ellipsoid_hole_line)) {
  //     stringstream s_stream(ellipsoid_hole_line);
  //     for (int d = 0; d < 3; d++) {
  //       string pt;
  //       getline(s_stream, pt, ',');
  //       file_hole_ellipsoid_points.push_back(std::stod(pt));
  //     }
  //     for (int d = 0; d < 3; d++) {
  //       string pt;
  //       getline(s_stream, pt, ',');
  //       file_hole_ellipsoid_normals.push_back(std::stod(pt));
  //     }
  //     string wt;
  //     getline(s_stream, wt, ',');
  //     file_hole_ellipsoid_weights.push_back(std::stod(wt));
  //   }
  //   myholefile.close();
  // }
  // int num_hole_ellipsoid_points = file_hole_ellipsoid_points.size() / 3;

  if (perturbation_parameters.size() == 0) {
    perturbation_parameters.push_back(0.0);
  }

  if (holes.size() == 0) {
    Hole hole;
    hole.center = PointVec(0.0, 0.0, 0.0);
    hole.radius = 0.5;
    hole.num_nodes = num_hole_cow_points;
    holes.push_back(hole);

    // hole.center = PointVec(0.0, -0.4, 0.3);
    // hole.radius = 0.2;
    // hole.num_nodes = num_hole_ellipsoid_points;
    // holes.push_back(hole);
    // hole.center = PointVec(0.0, -0.4, -0.3);
    // holes.push_back(hole);
  }


  int num_holes = holes.size();

  int total_points = num_outer_nodes + num_hole_cow_points;// +
                     // ((num_holes - 1) * num_hole_ellipsoid_points);
  for (int i = 0; i < num_outer_nodes; i++) {
    points.push_back(CIRCLE_RAD * file_points[3 * i]);
    points.push_back(CIRCLE_RAD * file_points[3 * i + 1]);
    points.push_back(CIRCLE_RAD * file_points[3 * i + 2]);
    // NOTE normals are just unscaled points for sphere
    normals.push_back(file_points[3 * i ]);
    normals.push_back(file_points[3 * i + 1]);
    normals.push_back(file_points[3 * i + 2]);
    weights.push_back(CIRCLE_RAD * CIRCLE_RAD * file_weights[i]);
  }

// Now the hole cow
  Hole hole = holes[0];

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


  // for (int hole_idx = 1; hole_idx < holes.size(); hole_idx++) {
  //   Hole hole = holes[hole_idx];
  //   double ang = 0;
  //   // if (hole_idx == 1) {
  //   //   ang = M_PI / 3.0;
  //   // } else if (hole_idx == 2) {
  //   //   ang = -M_PI / 3.0;
  //   // } else {
  //   //   ang = 0;
  //   // }
  //   for (int i = 0; i < hole.num_nodes; i++) {
  //     double x = file_hole_ellipsoid_points[3 * i];
  //     double y = cos(ang) * file_hole_ellipsoid_points[3 * i + 1] + sin(
  //                  ang) * file_hole_ellipsoid_points[3 * i + 2];
  //     double z = -sin(ang) * file_hole_ellipsoid_points[3 * i + 1] + cos(
  //                  ang) * file_hole_ellipsoid_points[3 * i + 2];
  //     double nx = file_hole_ellipsoid_normals[3 * i];
  //     double ny = cos(ang) * file_hole_ellipsoid_normals[3 * i + 1] + sin(
  //                   ang) * file_hole_ellipsoid_normals[3 * i + 2];
  //     double nz = -sin(ang) * file_hole_ellipsoid_normals[3 * i + 1] + cos(
  //                   ang) * file_hole_ellipsoid_normals[3 * i + 2];
  //     points.push_back(hole.center.a[0]
  //                      + hole.radius * x);
  //     points.push_back(hole.center.a[1]
  //                      +  hole.radius * y);
  //     points.push_back(hole.center.a[2]
  //                      +  hole.radius * z);
  //     normals.push_back(-nx);
  //     normals.push_back(-ny);
  //     normals.push_back(-nz);
  //     weights.push_back(hole.radius * hole.radius * file_hole_ellipsoid_weights[i]);
  //   }
  // }

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
 

void CowMesh::order_mesh_nodes(std::vector<double>* points, double cx,
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

void CowMesh::create_cross_section_border() {
  // get faces near plane for outer

  string line;

  string hole_cow_filename = "kern_interp/boundaries/meshes/cow_mesh_cross.txt";

  // string hole_ellipsoid_filename =
  //   "kern_interp/boundaries/meshes/ellipsoid_cross.txt";

  cross_section_border = new CrossSectionBorder();

  for (int i = 0; i < 100; i++) {
    double ang = (i / 100.0) * 2 * (M_PI);
    cross_section_border->sorted_2d_outer_knots.push_back(CIRCLE_RAD * cos(ang));
    cross_section_border->sorted_2d_outer_knots.push_back(CIRCLE_RAD * sin(ang));
  }

  for (int i = 0; i < holes.size(); i++) {
    Hole hole = holes[i];
    std::vector<double> hole_knot_points;
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

        double ang = 0;
        // if (i > 0) {
        //   if (i == 1) {
        //     ang = M_PI / 3.0;
        //   } else if (i == 2) {
        //     ang = -M_PI / 3.0;
        //   } else {
        //     ang = 0;
        //   }
        // }
        PointVec rotated(cos(ang)*pt_vec[0] + sin(ang)*pt_vec[1],
                         -sin(ang)*pt_vec[0] + cos(ang)*pt_vec[1]);
        hole_knot_points.push_back(rotated.a[0]*hole.radius);
        hole_knot_points.push_back(rotated.a[1]*hole.radius);
      }
      myholefile.close();
    } else {
      std::cout << "N0 File" << std::endl;
      exit(0);
    }
    // the args are called x and y although we are passing y and z - bad code
    order_mesh_nodes(&hole_knot_points, hole.center.a[1], hole.center.a[2]);
    std::vector<double> sorted_hole_knots;
    for (int j = 0; j < hole_knot_points.size(); j += 2) {
      sorted_hole_knots.push_back(hole.center.a[1] + hole_knot_points[j]);
      sorted_hole_knots.push_back(hole.center.a[2] + hole_knot_points[j + 1]);
    }
    cross_section_border->sorted_2d_hole_knots.push_back(sorted_hole_knots);
  }

  // populate and initialize cross section border
}


bool CowMesh::is_in_domain(const PointVec& a) const {
  return cross_section_border->is_in_domain(PointVec(a.a[1], a.a[2]));
}

}  // namespace kern_interp



// old code, delete when fixed
// void CowMesh::create_cross_section_border(int sz_param) {
//   // get faces near plane for outer

//   string line, outer_filename;
//   outer_filename = "kern_interp/boundaries/meshes/cow_mesh_cross.txt";

//   string hole_filename = "kern_interp/boundaries/meshes/ellipsoid_cross.txt";

//   std::vector<double> outer_2d_points;
//   cross_section_border = new CrossSectionBorder();

//   ifstream myfile(outer_filename);
//   if (myfile.is_open()) {
//     while (getline(myfile, line)) {
//       stringstream s_stream(line);
//       std::vector<double> pt_vec;
//       for (int d = 0; d < 2; d++) {
//         string pt;
//         getline(s_stream, pt, ',');
//         pt_vec.push_back(std::stod(pt));
//       }
//       outer_2d_points.push_back(pt_vec[0]);
//       outer_2d_points.push_back(pt_vec[1]);
//     }
//     myfile.close();
//   } else {
//     std::cout << "N0 File" << std::endl;
//     exit(0);
//   }


//   // sort points by angle
//   order_mesh_nodes(&outer_2d_points);

//   // std::sort(outer_2d_points.begin(), outer_2d_points.end(), comp_ang);
//   for (int i = 0; i < outer_2d_points.size(); i++) {
//     cross_section_border->sorted_2d_outer_knots.push_back(outer_2d_points[i]);
//   }
//   // same with holes, sort by angle about center
//   double point_idx = num_outer_nodes;
//   for (int i = 0; i < holes.size(); i++) {
//     Hole hole = holes[i];
//     std::vector<PointVec> hole_knots;

//     ifstream myholefile(hole_filename);
//     if (myholefile.is_open()) {
//       while (getline(myholefile, line)) {
//         stringstream s_stream(line);
//         std::vector<double> pt_vec;
//         for (int d = 0; d < 2; d++) {
//           string pt;
//           getline(s_stream, pt, ',');
//           pt_vec.push_back(std::stod(pt));
//         }

//         double ang = 0;
//         // if (i == 1) {
//         //   ang = M_PI / 3.0;
//         // } else if (i == 2) {
//         //   ang = -M_PI / 3.0;
//         // } else {
//         //   ang = 0;
//         // }
//         PointVec rotated(cos(ang)*pt_vec[0] + sin(ang)*pt_vec[1],
//                          -sin(ang)*pt_vec[0] + cos(ang)*pt_vec[1]);
//         hole_knots.push_back((rotated * hole.radius));
//       }
//       myholefile.close();
//     } else {
//       std::cout << "N0 File" << std::endl;
//       exit(0);
//     }
//     std::sort(hole_knots.begin(), hole_knots.end(), comp_ang);
//     std::vector<double> sorted_hole_knots;
//     for (int j = 0; j < hole_knots.size(); j++) {
//       PointVec pv = hole_knots[j];
//       sorted_hole_knots.push_back(hole.center.a[1] + pv.a[0]);
//       sorted_hole_knots.push_back(hole.center.a[2] + pv.a[1]);
//     }
//     cross_section_border->sorted_2d_hole_knots.push_back(sorted_hole_knots);
//     point_idx += hole.num_nodes;
//   }

//   // populate and initialize cross section border
// }
