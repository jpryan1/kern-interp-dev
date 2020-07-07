// Copyright 2019 John Paul Ryan
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include "kern_interp/boundaries/ex2boundary.h"

namespace kern_interp {

void Ex2Boundary::get_spline_points(std::vector<double>* x0_points,
                                    std::vector<double>* x1_points) {


  for (int i = 0; i <= 3; i++) {
    x0_points->push_back(0);
    x1_points->push_back(3 - 0.5 * (i / 3.0));
  }


  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25 * (i / 3.0));
    x1_points->push_back(2.5);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25);
    x1_points->push_back(2.5 + 0.25 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25 - 0.25 * (i / 3.0));
    x1_points->push_back(2.75);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5);
    x1_points->push_back(2.75 - 0.5 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5 + 0.5 * (i / 3.0));
    x1_points->push_back(2.25);
  }



  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0);
    x1_points->push_back(2.25 - 0.5 * (i / 3.0));
  }


  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5 * (i / 3.0));
    x1_points->push_back(1.75);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5);
    x1_points->push_back(1.75 - 0.5 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5 + 0.25 * (i / 3.0));
    x1_points->push_back(1.25);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25);
    x1_points->push_back(1.25 + 0.25 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25 + 0.25 * (i / 3.0));
    x1_points->push_back(1.5);
  }


  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0);
    x1_points->push_back(1.5 - 0.5 * (i / 3.0));
  }



  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5 * (i / 3.0));
    x1_points->push_back(1.0);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5);
    x1_points->push_back(1.0 - 0.5 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.5 + 0.25 * (i / 3.0));
    x1_points->push_back(0.5);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25);
    x1_points->push_back(0.5 + 0.25 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(-0.25 + 0.25 * (i / 3.0));
    x1_points->push_back(0.75);
  }

  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0);
    x1_points->push_back(0.75 - 0.75 * (i / 3.0));
  }

// bottom wall
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5 * (i / 3.0));
    x1_points->push_back(0);
  }




  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5);
    x1_points->push_back(0.5 * (i / 3.0));
  }


  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5 + 0.25 * (i / 3.0));
    x1_points->push_back(0.5);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.75);
    x1_points->push_back(0.5- 0.25 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.75+ 0.25 * (i / 3.0));
    x1_points->push_back(0.25);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(1);
    x1_points->push_back(0.25+0.5 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(1- 0.5 * (i / 3.0));
    x1_points->push_back(0.75);
  }



  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5);
    x1_points->push_back(0.75 + 1.0 * (i / 3.0));
  }


  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5+0.25 * (i / 3.0));
    x1_points->push_back(1.75);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.75);
    x1_points->push_back(1.75 - 0.25 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.75 + 0.25 * (i / 3.0));
    x1_points->push_back(1.5);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(1.0);
    x1_points->push_back(1.5 + 0.5 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(1.0-0.5 * (i / 3.0));
    x1_points->push_back(2.0);
  }


  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5);
    x1_points->push_back(2.0 + 0.25 * (i / 3.0));
  }



  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5+0.5 * (i / 3.0));
    x1_points->push_back(2.25);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(1.0);
    x1_points->push_back(2.25 + 0.5 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(1.0- 0.25 * (i / 3.0));
    x1_points->push_back(2.75);
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.75);
    x1_points->push_back(2.75 - 0.25 * (i / 3.0));
  }
  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.75 - 0.25 * (i / 3.0));
    x1_points->push_back(2.5);
  }

  for (int i = 1; i <= 3; i++) {
    x0_points->push_back(0.5);
    x1_points->push_back(2.5 + 0.5 * (i / 3.0));
  }


// top wall
  for (int i = 1; i < 3; i++) {
    x0_points->push_back(0.5-0.5 * (i / 3.0));
    x1_points->push_back(3.0);
  }

// // top wall
//   for (int i = 0; i < 12; i++) {
//     x0_points->push_back(-1 + 3 * (i / 12.0));
//     x1_points->push_back(0.25);
//   }
//   for (int i = 0; i < 2; i++) {
//     x0_points->push_back(2);
//     x1_points->push_back(0.25 + 0.5 * i / 2.0);
//   }

//   for (int i = 0; i < 12; i++) {
//     x0_points->push_back(2 - 3 * (i / 12.0));
//     x1_points->push_back(0.75);
//   }
//   for (int i = 0; i < 2; i++) {
//     x0_points->push_back(-1);
//     x1_points->push_back(0.75 - 0.5 * i / 2.0);
//   }
}


void Ex2Boundary::initialize(int N, BoundaryCondition bc) {
  points.clear();
  normals.clear();
  weights.clear();
  curvatures.clear();
  holes.clear();

  // if (perturbation_parameters.size() == 0) {
  //   for (int i = 0; i <= 7; i++) perturbation_parameters.push_back(0.5);
  // }

  int NUM_CIRCLE_POINTS = (N / 4) / 8;

  Hole circle;
  circle.radius = 0.05;
  circle.num_nodes =  NUM_CIRCLE_POINTS;
  // for (int i = 2; i <= 3; i++) {
    circle.center = PointVec(0.25, 0.5);
    // holes.push_back(circle);
    // circle.center = PointVec(0.25, 2.0);
    holes.push_back(circle);
  // }


  std::vector<double> outer_x0_spline_points, outer_x1_spline_points;
  get_spline_points(&outer_x0_spline_points, &outer_x1_spline_points);


  int OUTER_NUM_SPLINE_POINTS = outer_x0_spline_points.size();
  int OUTER_NODES_PER_SPLINE = (3 * N / 4) / OUTER_NUM_SPLINE_POINTS;
  num_outer_nodes = OUTER_NODES_PER_SPLINE * OUTER_NUM_SPLINE_POINTS;

  std::vector<std::vector<double>> outer_x0_cubics, outer_x1_cubics;
  get_cubics(outer_x0_spline_points, outer_x1_spline_points,
             &outer_x0_cubics, &outer_x1_cubics);

  interpolate(false, OUTER_NODES_PER_SPLINE,
              outer_x0_cubics, outer_x1_cubics);

  for (int i = 0; i < holes.size(); i++) {
    Hole circle = holes[i];
    for (int i = 0; i < NUM_CIRCLE_POINTS; i++) {
      double ang = (2.0 * M_PI * i) / NUM_CIRCLE_POINTS;
      points.push_back(circle.center.a[0] + circle.radius * cos(ang));
      points.push_back(circle.center.a[1] + circle.radius * sin(ang));
      normals.push_back(-cos(ang));
      normals.push_back(-sin(ang));
      curvatures.push_back(-1.0 / circle.radius);
      weights.push_back(2.0 * M_PI * circle.radius / NUM_CIRCLE_POINTS);
    }
  }

  if (bc == BoundaryCondition::DEFAULT) {
    boundary_values = ki_Mat(weights.size() * 2, 1);
    apply_boundary_condition(0, weights.size(), LEFT_TO_RIGHT_FLOW);
  } else {
    set_boundary_values_size(bc);
    apply_boundary_condition(0, weights.size(), bc);
  }
}


}  // namespace kern_interp
