// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include "kern_interp/boundaries/cross_section_border.h"

namespace kern_interp {

void CrossSectionBorder::get_spline_points(std::vector<double>* x0_points,
    std::vector<double>* x1_points) {
  for (int i = 0; i < sorted_2d_outer_knots.size(); i += 2) {
    x0_points->push_back(sorted_2d_outer_knots[i]);
    x1_points->push_back(sorted_2d_outer_knots[i + 1]);
  }
}

void CrossSectionBorder::get_inner_spline_points(
  std::vector<std::vector<double>>* x0_points,
  std::vector<std::vector<double>>* x1_points) {
  for (int hole_idx = 0; hole_idx < sorted_2d_hole_knots.size(); hole_idx++) {
    std::vector<double> hole_spline_knots = sorted_2d_hole_knots[hole_idx];
    std::vector<double> x0_inner_knots, x1_inner_knots;
    for (int i = 0; i < hole_spline_knots.size(); i += 2) {
      x0_inner_knots.push_back(hole_spline_knots[i]);
      x1_inner_knots.push_back(hole_spline_knots[i + 1]);
    }
    x0_points->push_back(x0_inner_knots);
    x1_points->push_back(x1_inner_knots);
  }
}


PointVec get_hole_center(const std::vector<double>& x0,
                         const std::vector<double>& x1) {
  double avg0 = 0.;
  double avg1 = 0.;
  for (int i = 0; i < x0.size(); i++) {
    avg0 += x0[i];
    avg1 += x1[i];
  }
  return PointVec(avg0 / x0.size(), avg1 / x1.size());
}


double get_hole_radius(const std::vector<double>& x0,
                       const std::vector<double>& x1,
                       const PointVec& center) {
  double rad = 0.;
  for (int i = 0; i < x0.size(); i++) {
    rad = std::max(rad, sqrt(pow(center.a[0] - x0[i], 2) + pow(center.a[1] - x1[i],
                             2)));
  }
  return rad + 1e-3;
}


void CrossSectionBorder::initialize(int N, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  holes.clear();
  int OUTER_NODES_PER_SPLINE = 3;
  int INNER_NODES_PER_SPLINE = 3;


  // Populate num_outer_nodes...done
  // Populate num_nodes for all interior holes...done
  // interior holes center set to centroid, radius set to max dist from centroid...done
  // call get_cubics and interpolate for outer and for all holes...done

  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;

  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);
  interpolate(false, OUTER_NODES_PER_SPLINE,
              outer_x0_cubics, outer_x1_cubics);
  num_outer_nodes = OUTER_NODES_PER_SPLINE *  outer_x0_spline_points.size();

  std::vector<std::vector<double>> inner_x0_spline_points, inner_x1_spline_points;
  get_inner_spline_points(&inner_x0_spline_points, &inner_x1_spline_points);
  for (int hole_idx = 0; hole_idx < inner_x1_spline_points.size(); hole_idx++) {
    std::vector<double> x0_inner_knots = inner_x0_spline_points[hole_idx];
    std::vector<double> x1_inner_knots = inner_x1_spline_points[hole_idx];
    std::vector<std::vector<double>> inner_x0_cubics, inner_x1_cubics;
    get_cubics(x0_inner_knots, x1_inner_knots,
               &inner_x0_cubics, &inner_x1_cubics);
    interpolate(true, INNER_NODES_PER_SPLINE, inner_x0_cubics, inner_x1_cubics);
    Hole hole;
    hole.center = get_hole_center(x0_inner_knots, x1_inner_knots);
    hole.radius = get_hole_radius(x0_inner_knots, x1_inner_knots,
                                  hole.center);
    hole.num_nodes = INNER_NODES_PER_SPLINE * x0_inner_knots.size();
    holes.push_back(hole);
  }

  std::ofstream bound_out;
  bound_out.open("output/data/cross_section_bound.txt");
  for (int i = 0; i < num_outer_nodes*2; i += 2) {
    bound_out << points[i] << "," << points[i + 1]
              << std::endl;
  }
  bound_out.close();
}


}  // namespace kern_interp
